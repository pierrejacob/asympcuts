"""Place beside baselines.py and run: python make_k1_data.py

One-time setup: python -m pip install "pyarrow<25"
Downloads Breunig, Liu and Yu's original bank and creates data/ from scratch.
Source: https://doi.org/10.5281/zenodo.14015435 (archive license: CC BY 4.0).
"""
from pathlib import Path
from io import BytesIO
from urllib.request import urlopen
from zipfile import ZipFile
import hashlib
import json
import numpy as np
from pyarrow.feather import read_table
from scipy.optimize import brentq
from scipy.special import expit
from baselines import logistic_map

folder = Path(__file__).resolve().parent
archive = folder / "Replication.zip"
target = folder / "data"
if target.exists():
    raise FileExistsError(f"{target} already exists; move it before running again.")

# 1. Download once, verify the published archive, and read the original Feather bank.
url = "https://zenodo.org/api/records/14015435/files/Replication.zip/content"
archive_md5 = "8013e263ed18da795c5ccb430be489ce" # a fingerprint that will help us check that the file we downloaded matches the expected file
if not archive.exists():
    print("Downloading the replication archive (about 19 MB)...", flush=True)
    with urlopen(url, timeout=60) as response:
        archive.write_bytes(response.read())
if hashlib.md5(archive.read_bytes()).hexdigest() != archive_md5:
    raise ValueError("Replication.zip differs from the published archive.")
member = "Replication/Data/psid_generated.feather"
with ZipFile(archive) as zipped:
    table = read_table(BytesIO(zipped.read(member)))
data = {name: table[name].to_numpy() for name in table.column_names}

# 2. Prepare the nine covariates and binary potential outcomes; keep every row.
# We omit nodegree, as in the experiment desciption.
data["emp74"], data["emp75"] = data["re74"] > 0, data["re75"] > 0 # calculating employment indicators
names = ["age", "education", "black", "hispanic", "married",
         "re74", "re75", "emp74", "emp75"]
raw_x = np.array([data[name] for name in names], dtype=float).T
original_a = data["t"].astype(np.int32) # save the original synthetic treatment labels as integers
observed, counterfactual = data["re78"] > 0, data["re78_cf"] > 0 # calculating employment indicators in 1978.
# re78 is factual under the ORIGINAL treatment label; re78_cf is counterfactual.
y0 = np.where(original_a == 0, observed, counterfactual).astype(np.int32) # Construct Y(0) for everyone: factual if originally untreated, counterfactual otherwise.
y1 = np.where(original_a == 1, observed, counterfactual).astype(np.int32) # Construct Y(1) for everyone: factual if originally treated, counterfactual otherwise.


# 3. Standardize using the full bank, then fit the original treatment labels.
center, scale = raw_x.mean(axis=0), raw_x.std(axis=0)
x = (raw_x - center) / scale
z = np.column_stack((np.ones(len(x)), x)) # add a column of ones for the intercept
prior_sd = [10] + [5] * 9  # Independent zero-mean Normal priors.
fitted, _ = logistic_map(z, original_a, np.ones(len(x)), prior_sd, np.zeros(10)) # fits logistic regression to all original synthetic treatment labels

# 4. Rescale the fitted slopes to give the linear predictor population SD 1.
slopes = fitted[1:] / (x @ fitted[1:]).std()
score = x @ slopes

# 5. Calibrate the intercept to mean treatment probability 0.30.
intercept = brentq(lambda a: expit(a + score).mean() - 0.30, -50, 50)
alpha = np.r_[intercept, slopes]
e0 = expit(intercept + score) # calculates new treatment probabilities
truth = float(e0 @ (y1 - y0) / e0.sum()) # calculates the true population ATT

# 6. Save the three files read by data.py.
target.mkdir()
bank_path = target / "population.npz"
np.savez_compressed(bank_path, X=raw_x, y0=y0, y1=y1, original_treatment=original_a)
config = dict(
    covariates=names, center=center.tolist(), scale=scale.tolist(),
    alpha0=alpha.tolist(), reference_coefficients=fitted.tolist(),
    reference_prior_sd=prior_sd, propensity_mean=float(e0.mean()),
    population_att=truth, population_npz_sha256=hashlib.sha256(bank_path.read_bytes()).hexdigest(),
    source_url=url, source_archive_md5=archive_md5, source_member=member,
    source_attribution="Breunig, Liu and Yu replication archive; WGAN bank from Athey et al.",
    source_archive_license="CC BY 4.0",
    calibration_rule="Logistic MAP on original labels; predictor SD 1; mean propensity 0.30",
)
np.save(target / "true_propensity.npy", e0) # saves the million true propensity scores as an array
(target / "population_config.json").write_text(json.dumps(config, indent=2) + "\n")
print(f"Created {target}\nMean treatment probability: {e0.mean():.12f}\nPower ATT: {truth:.12f}")
