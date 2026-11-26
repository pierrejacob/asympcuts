"""Run datasets in parallel in a NEW folder. A fitting error stops the study."""
import os
for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"):
    os.environ[name] = "1" # sets several numerical libraries to use one computational thread each, prevents thread oversubscription

from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import get_context
import csv
import importlib.metadata
import json
import platform
import subprocess
import time
import numpy as np
from threadpoolctl import threadpool_limits
from data import load_population, make_data
from baselines import fit_baselines
from summarize import summarize

# All simulation controls are here.
DATASET_IDS = range(2001, 2002)  # One ID; full study: range(2001, 3001). These are the identifiers that enter the random seed calculation in file data.py. Mode full introduced below will change this range to range(2001,3001).
SCENARIOS = ["power", "size"]
N = 2675  # Default sample size; use --n 400 2675 to run both sizes.
ROOT_SEED = 2026091801
BASELINE_DRAWS = 1000 # Requests 1000 PBMI draws and 1000 pairs-bootstrap draws per dataset.
YIU_DRAWS, BURNIN, CHAINS, TREES = 16000, 4000, 4, 200
WORKERS = 8  # Dataset workers; each uses four BART chains during its outcome fit.
OUTPUT = Path(__file__).resolve().parent / "results_simple"
RSCRIPT = os.environ.get("RSCRIPT", "Rscript") # Selects the command used to launch R.
COLUMNS = ["n", "scenario", "dataset_id", "method", "truth", "estimate", "lower", "upper"]

def initialize_worker(settings):
    """Each worker loads its own population once, then reuses it across datasets."""
    global POPULATION, SETTINGS
    threadpool_limits(limits=1) # limits the supported numerical libraries to one thread and reinforces the settings from the top of this file
    SETTINGS = settings
    POPULATION = load_population()


def run_one_dataset(job):
    scenario, dataset_id = job
    s = SETTINGS
    started = time.perf_counter()
    package = Path(__file__).resolve().parent
    folder = Path(s["output"]) / scenario / str(dataset_id)
    folder.mkdir(parents=True)  # Every job owns a different folder.
    print("Starting:", scenario, dataset_id, flush=True)
    try:
        d = make_data(POPULATION, dataset_id, scenario, s["n"], s["root_seed"])
        observed = folder / "observed.csv"
        np.savetxt(observed, np.column_stack([d["x"], d["a"], d["y"]]), fmt="%.17g",
                   delimiter=",", comments="", header=",".join(POPULATION["names"] + ["A", "Y"]))  # Saves the observed dataset as CSV: n rows and 11 columns (nine covariates, A and Y).
        command = [s["rscript"], "--vanilla", str(package / "yiu.R"), str(observed), str(folder)]
        command += [str(v) for v in (d["seed"], s["yiu_draws"], s["burnin"], s["chains"], s["trees"])]
        with (folder / "yiu.log").open("w") as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
        baseline = fit_baselines(d["x"], d["a"], d["y"], dataset_id, s["root_seed"], s["baseline_draws"])
        draws = baseline["draws"]
        np.savetxt(folder / "baseline_draws.csv", np.column_stack([draws["pbmi"], draws["bootstrap"]]),
                   delimiter=",", header="pbmi,bootstrap", comments="", fmt="%.17g") 
        with (folder / "fit_info.csv").open("a") as info:
            info.write("baseline_max_scaled_gradient," + str(baseline["max_scaled_gradient"]) + "\n") # adds an optimisation diagnostic - records the largest absolute gradient component, after the objective's scaling, encountered at the final solutions across the baseline optimisations
        draws["yiu"] = np.loadtxt(folder / "yiu_draws.csv", delimiter=",", skiprows=1)
        rows = []
        for method, values in draws.items():
            if not np.isfinite(values).all():
                raise ValueError("Nonfinite draws for " + method)
            lower, upper = np.quantile(values, [0.025, 0.975]) # calculates interval endpoints as in Sec. 5.2.3
            point = values.mean() if method == "yiu" else baseline["point_estimate"] # chooses the point estimate: the mean of corrected ATT draws for Yiu et al. and the MAP estimate for the baseline methods
            rows.append(dict(n=s["n"], scenario=scenario, dataset_id=dataset_id, method=method,
                             truth=d["truth"], estimate=point, lower=lower, upper=upper))
        with (folder / "summary.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=COLUMNS)
            writer.writeheader()
            writer.writerows(rows)
        elapsed = time.perf_counter() - started
        (folder / "elapsed_seconds.txt").write_text(str(elapsed) + "\n")
        return rows, elapsed
    except Exception as error:
        raise RuntimeError(f"Fit failed at {scenario}, dataset {dataset_id}; inspect {folder}") from error


def main():
    started = time.perf_counter()
    jobs = [(scenario, dataset_id) for scenario in SCENARIOS for dataset_id in DATASET_IDS]
    if WORKERS < 1 or not jobs:
        raise ValueError("Use at least one worker and one dataset.")
    OUTPUT.mkdir(exist_ok=False)  # Creates the output directory. No automatic resume and no overwriting an earlier study.
    settings = dict(dataset_ids=list(DATASET_IDS), scenarios=SCENARIOS, n=N,
                    root_seed=ROOT_SEED, baseline_draws=BASELINE_DRAWS,
                    yiu_draws=YIU_DRAWS, burnin=BURNIN, chains=CHAINS, trees=TREES,
                    workers=WORKERS, output=str(OUTPUT), rscript=RSCRIPT,
                    python=platform.python_version(), versions={name: importlib.metadata.version(name)
                    for name in ("numpy", "scipy", "threadpoolctl")})
    (OUTPUT / "settings.json").write_text(json.dumps(settings, indent=2))
    check = 'stopifnot(requireNamespace("BART"), requireNamespace("posterior"), requireNamespace("rstanarm")); sessionInfo()'
    with (OUTPUT / "R_versions.txt").open("w") as handle:
        subprocess.run([RSCRIPT, "--vanilla", "-e", check], check=True, stdout=handle)
    print(f"Running {len(jobs)} datasets with up to {WORKERS} workers", flush=True)
    with (OUTPUT / "per_dataset.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        with ProcessPoolExecutor(max_workers=min(WORKERS, len(jobs)), mp_context=get_context("spawn"),
                                 initializer=initialize_worker, initargs=(settings,)) as pool:  #spawn = fresh interpreter per worker
            completed = 0
            # Submit only one small batch: after an error, no later batch starts.
            for start in range(0, len(jobs), WORKERS):
                futures = [pool.submit(run_one_dataset, job) for job in jobs[start:start + WORKERS]]
                for future in as_completed(futures):
                    try:
                        rows, seconds = future.result()
                    except Exception:
                        print("A fit failed. No further batches will start. Waiting for the other "
                              "datasets in this batch to finish, then exiting with an error.", flush=True)
                        raise
                    writer.writerows(rows)  #Only the parent writes the combined table, avoiding several workers working on it simultaneously.
                    handle.flush()
                    completed += 1
                    print(f"Completed {completed}/{len(jobs)}: {rows[0]['scenario']} "
                          f"{rows[0]['dataset_id']} ({seconds:.1f} seconds)", flush=True)
    summarize(OUTPUT) #Reached only if every requested fit succeeded. Writes aggregate.csv and the diagnostic summaries.
    elapsed = time.perf_counter() - started
    (OUTPUT / "elapsed_seconds.txt").write_text(str(elapsed) + "\n")
    print("Finished:", OUTPUT, "\nTotal wall-clock hours:", elapsed / 3600)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["example", "pilot", "full"])
    parser.add_argument("--n", type=int, nargs="+", choices=[400, 2675], default=[N],
                        help="Sample sizes to run, in order (default: %(default)s)")
    args = parser.parse_args()
    if len(set(args.n)) != len(args.n):
        parser.error("Specify each sample size only once.")

    if args.mode is not None:  # Otherwise keep DATASET_IDS and OUTPUT set above.
        stop, folder = {"example": (2002, "results_simple"),
                        "pilot": (2021, "results_pilot20"),
                        "full": (3001, "results_full1000")}[args.mode]
        DATASET_IDS = range(2001, stop) 
        OUTPUT = Path(__file__).resolve().parent / folder

    # Separate result folders prevent mixing the two sample sizes.
    base_output = OUTPUT
    output_folders = [base_output.with_name(f"{base_output.name}_n{n}")
                      for n in args.n]
    for folder in output_folders:
        if folder.exists():
            parser.error(f"Output already exists: {folder}. No runs were started.")

    # Each main() call finishes and closes its worker pool before the next size.
    for N, OUTPUT in zip(args.n, output_folders):
        print(f"Starting sample size n={N}", flush=True)
        main()
