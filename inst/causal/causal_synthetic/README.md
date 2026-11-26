This package reproduces the ATT experiment in Section 5.2.2 and Appendix E. 

## 1. Setup

Use a Unix system with at least four CPU cores, Python 3.12, and R with its
recommended packages and source-package compilation tools. Install

- **Python:** NumPy 2.3.5, SciPy 1.17.0, and threadpoolctl 3.6.0, as listed in
  `requirements.txt`.
- **R:** Rcpp, posterior, rstanarm, rstan, and BART 2.9.10, including their dependencies.
The experiment requires BART version 2.9.10.
  

The paper used R 4.3.3, posterior 1.5.0, rstanarm 2.32.1, rstan 2.32.5. and BART 2.9.10.
Use these versions for the closest numerical reproduction. Installing R packages
from source requires the appropriate compiler tools.

Run the commands below from the extracted `yiu_att_replication` directory. If your Python
3.12 interpreter is called `python3.12`, use that in place of `python`.

## 2. Run the experiment

Before the first run, install the Python dependencies and build the simulation
population. Run these commands once. The script creates `data/` so that directory
should not already exist. 

```bash
python -m pip install -r requirements.txt "pyarrow<25"
python make_k1_data.py
```

The script downloads `Replication.zip` from the Breunig, Liu and Yu replication archive on Zenodo
(https://zenodo.org/records/14015435), checks its checksum, and creates `data/population.npz`, 
`data/true_propensity.npy` and `data/population_config.json`. If you already have `Replication.zip`,
or prefer to download it manually, place it in `yiu_att_replication/` without unzipping it. 
Then run `python make_k1_data.py` to use that copy.

Then check the installation using one dataset per scenario and sample size:

```bash
python run.py --mode example --n 400 2675
```

Then reproduce the full experiment:

```bash
python run.py --mode full --n 400 2675
```

Results go to `results_full1000_n400/` and `results_full1000_n2675/`, running the
sample sizes sequentially. Each has 1,000 datasets per scenario (IDs 2001–3000).
The supplied settings match the paper:

- 1,000 draws for PBMI and 1,000 pairs-bootstrap replicates per dataset;
- four chains for each Yiu nuisance model, each with 4,000 warmup iterations
  and 4,000 retained draws, yielding 16,000 corrected ATT draws;
- 200 BART trees and empirical 2.5%–97.5% interval endpoints.

Keep the seed and statistical settings unchanged. Reduce `WORKERS = 8` in
`run.py` for smaller machines; each worker can use four BART processes. Allow
several tens of GB for outputs. An optional timing pilot uses
`python run.py --mode pilot --n 400 2675` (20 datasets per scenario and size).

Each run requires a new output folder. Before repeating a run, rename or move 
its previous output folder. Interrupted runs are not resumed automatically.

## 3. Obtain the paper's tables

Summaries are written automatically after all datasets for a sample size finish.
In each full-run folder:

| File | Use |
| --- | --- |
| `aggregate.csv` | Main table: `coverage`, `rejection_rate`, and `mean_interval_length`, plus simulation Monte Carlo standard errors. |
| `diagnostic_summary.csv` | Appendix diagnostic counts: use `diagnostic_flag_datasets` for the number of datasets with any warning in each component. |
| `stan_diagnostic_summary.csv` | Additional Stan divergence, maximum-tree-depth, E-BFMI, and unavailable-diagnostic counts. |
| `per_dataset.csv` | Estimates and interval endpoints for each dataset and method. |

`power` denotes the nonzero-effect scenario; `size` denotes the sharp null.
`rejection_rate` is respectively power or size. Methods are `pbmi`, `bootstrap`
(pairs bootstrap), and `yiu`. 

Flags do not exclude datasets or draws.
Detailed results and traces are under `scenario/dataset_id/`.
`settings.json` and `R_versions.txt` record the run configuration and software.


## 4. Data provenance

The population comes from the Athey et al. synthetic bank distributed in the
[Breunig, Liu and Yu replication archive](https://doi.org/10.5281/zenodo.14015435)
(CC BY 4.0). `make_k1_data.py` reads `Replication/Data/psid_generated.feather` from that archive,
defines employment as positive 1978 earnings and calibrates the
treatment-assignment model.

Original Yiu et al. routines are in `original_yiu/`; adaptations are in `yiu.R`
and `helpers/`.
