"""The frozen population is used only to generate observed data and the truth."""
from pathlib import Path
import json
import numpy as np
from scipy.special import expit


def load_population():
    folder = Path(__file__).resolve().parent / "data"
    config = json.loads((folder / "population_config.json").read_text())
    with np.load(folder / "population.npz", allow_pickle=False) as bank:
        x = (bank["X"] - config["center"]) / config["scale"] # loads and standardizes the covariates
        y0 = bank["y0"].astype(np.int8) # loads potential outcomes under control
        y1 = bank["y1"].astype(np.int8) # loads poential outcomes under treatment
    e0 = np.load(folder / "true_propensity.npy", allow_pickle=False) # loads the true propensities
    alpha = np.array(config["alpha0"])
    if not np.allclose(e0, expit(alpha[0] + x @ alpha[1:]), atol=2e-14, rtol=0):
        raise ValueError("The saved treatment probabilities do not match the model")
    truth = float(np.dot(e0, y1.astype(float) - y0.astype(float)) / e0.sum()) # calculates the true population ATT
    return dict(x=x, y0=y0, y1=y1, e0=e0, truth=truth, names=config["covariates"])


def make_data(population, dataset_id, scenario, n, root_seed):
    p = population
    # Stage 0: data; stage 17: the seed passed to R.
    rng = np.random.default_rng(np.random.SeedSequence([root_seed, dataset_id, 0]))
    rows = rng.integers(len(p["e0"]), size=n) # selects the rows to appear in the dataset
    a = (rng.random(n) < p["e0"][rows]).astype(np.int8) # generates treatment assignments by comparing i.i.d. uniform(0,1) draws to propensities
    if scenario == "power":
        y = np.where(a, p["y1"][rows], p["y0"][rows]) # constructs the observed outcome
        truth = p["truth"] # the true population ATT in the power scenario
    elif scenario == "size":
        y = p["y0"][rows]  # Both potential outcomes equal y0: a sharp null.
        truth = 0.0 # the true population ATT in the size scenario
    else:
        raise ValueError("Scenario must be power or size")
    seed_value = np.random.SeedSequence([root_seed, dataset_id, 17]).generate_state(1)[0] # generates one reproducible integer for the R calculations
    seed = 1 + int(seed_value) % 2147483643 # maps that integer into a range accepted by the R code
    return dict(x=p["x"][rows], a=a, y=y, truth=truth, seed=seed) 
