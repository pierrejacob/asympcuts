"""Coverage, size/power and interval length; no fits are filtered on diagnostics."""
import csv
import json
from pathlib import Path
import numpy as np


def summarize(output):
    output = Path(output)
    settings = json.loads((output / "settings.json").read_text())  # reads the settings file created by run.py
    with (output / "per_dataset.csv").open() as handle:
        rows = list(csv.DictReader(handle))
        # Guard against mixing results from different sample sizes in one folder.
    if any("n" in row and int(row["n"]) != settings["n"] for row in rows):
        raise ValueError("Result rows have a different sample size from settings.json")
    summaries = []
    for scenario in settings["scenarios"]:
        for method in ("pbmi", "bootstrap", "yiu"):
            group = [r for r in rows if r["scenario"] == scenario and r["method"] == method]
            if sorted(int(r["dataset_id"]) for r in group) != sorted(settings["dataset_ids"]):
                raise ValueError("Missing or duplicate results: " + scenario + " / " + method)
            values = np.array([[float(r[k]) for k in ("truth", "estimate", "lower", "upper")] for r in group])
            if not np.isfinite(values).all() or np.any(values[:, 2] > values[:, 3]):
                raise ValueError("Nonfinite estimate or invalid interval")
            truth, estimate, lower, upper = values.T
            number = len(group)
            coverage = np.mean((lower <= truth) & (truth <= upper))
            rejection = np.mean((lower > 0) | (upper < 0))
            widths = upper - lower
            summaries.append(dict(n=settings["n"], scenario=scenario, method=method, repetitions=number,
                coverage=coverage, rejection_rate=rejection,
                interpretation="size" if scenario == "size" else "power",
                mean_interval_length=widths.mean(),
                coverage_mcse=np.sqrt(coverage * (1 - coverage) / number),
                rejection_mcse=np.sqrt(rejection * (1 - rejection) / number),
                length_mcse=widths.std(ddof=1) / np.sqrt(number) if number > 1 else None))
    with (output / "aggregate.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    summarize_diagnostics(output, settings)
    summarize_stan_diagnostics(output, settings)

def summarize_diagnostics(output, settings):
    """Count flagged DATASETS, keeping outcome and propensity results separate."""
    rows = []
    for scenario in settings["scenarios"]:
        for dataset_id in settings["dataset_ids"]:
            folder = output / scenario / str(dataset_id)
            with (folder / "diagnostics.csv").open() as handle:
                diagnostics = list(csv.DictReader(handle))
            # Missing diagnostics must not be counted as a passed check.
            expected = ["m0_mean_treated", "m0_mean_control",
                        "e_mean_treated", "e_mean_control",
                        "max_control_odds", "att_without_bb", "corrected_att"]
            expected += [f"alpha_{j}" for j in range(10)]
            if sorted(r["variable"] for r in diagnostics) != sorted(expected):
                raise ValueError("Missing or duplicate diagnostics in " + str(folder))
            rows += [dict(scenario=scenario, dataset_id=dataset_id, **r) for r in diagnostics]
    flags = ["rhat_flag", "ess_flag", "endpoint_mcse_flag", "unavailable_flag", "diagnostic_flag"] # five kinds of flags recorded by diagnostics.R
    summaries = []
    for scenario in settings["scenarios"]:
        components = sorted({r["component"] for r in rows if r["scenario"] == scenario})
        for component in components:
            group = [r for r in rows if r["scenario"] == scenario and r["component"] == component]
            record = dict(n=settings["n"], scenario=scenario, component=component, datasets=len(settings["dataset_ids"]))
            for flag in flags:
                if any(r[flag] not in ("TRUE", "FALSE") for r in group):
                    raise ValueError("Missing diagnostic flag: " + flag)
                record[flag + "_datasets"] = len({r["dataset_id"] for r in group if r[flag] == "TRUE"}) # count flagged datasets for this particular type of flag
            summaries.append(record)
    with (output / "diagnostic_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)

def summarize_stan_diagnostics(output, settings):
    """Count datasets with Stan warnings or unavailable diagnostics."""
    output = Path(output)
    threshold = 0.3  # Flag E-BFMI below this value in any chain.
    summaries = []

    for scenario in settings["scenarios"]:
        flags_by_dataset = []

        for dataset_id in settings["dataset_ids"]:
            path = output / scenario / str(dataset_id) / "fit_info.csv"
            info = {}

            # Missing files or values remain unavailable, rather than passing.
            if path.exists():
                with path.open() as handle:
                    for row in csv.DictReader(handle):
                        name = row["item"]
                        if name in info:
                            raise ValueError("Duplicate item in " + str(path))
                        try:
                            info[name] = float(row["value"])
                        except (TypeError, ValueError):
                            info[name] = np.nan  # Includes R's "NA".

            divergences = info.get("stan_divergences", np.nan)
            treedepth = info.get("stan_max_treedepth_hits", np.nan)
            ebfmi = np.array([
                info.get(f"stan_ebfmi_{chain}", np.nan)
                for chain in range(1, settings["chains"] + 1)
            ])

            # Counts must be finite, nonnegative integers.
            div_ok = (np.isfinite(divergences) and divergences >= 0
                      and divergences == np.floor(divergences))
            depth_ok = (np.isfinite(treedepth) and treedepth >= 0
                        and treedepth == np.floor(treedepth))
            ebfmi_ok = np.isfinite(ebfmi) & (ebfmi >= 0)

            flags = dict(
                divergence=div_ok and divergences > 0,
                treedepth=depth_ok and treedepth > 0,
                low_ebfmi=np.any(ebfmi[ebfmi_ok] < threshold),
                unavailable_divergence=not div_ok,
                unavailable_treedepth=not depth_ok,
                unavailable_ebfmi=not np.all(ebfmi_ok),
            )
            # A dataset can have both a warning and unavailable diagnostics.
            flags["any_stan_flag"] = any(flags.values())
            flags_by_dataset.append(flags)

        record = dict(n=settings["n"], scenario=scenario,
                      datasets=len(settings["dataset_ids"]),
                      ebfmi_threshold=threshold)
        for name in flags_by_dataset[0]:
            record[name + "_datasets"] = sum(
                int(flags[name]) for flags in flags_by_dataset
            )
        summaries.append(record)

    with (output / "stan_diagnostic_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
            
if __name__ == "__main__":
    import sys
    summarize(sys.argv[1])  # Gives that argument to summarize. Example: python summarize.py results_simple
