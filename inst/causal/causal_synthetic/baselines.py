"""PBMI and pairs bootstrap for ATT, using ordinary SciPy logistic MAP fits."""
import numpy as np
from scipy.optimize import minimize, root
from scipy.special import expit


def logistic_map(x, y, weights, prior_sd, start): # takes the design matrix, binary response, observation weights, prior standard deviation for the coefficients, initial coefficient vector for the numerical optimiser
    """Minimize the negative weighted log posterior; check its final gradient."""
    precision = 1 / np.asarray(prior_sd) ** 2
    scale = 1 + weights.sum()
    sign = 1 - 2 * y

    def objective(beta):
        # logaddexp(0, t) computes log(1 + exp(t)) without overflow.
        loss = weights @ np.logaddexp(0, sign * (x @ beta)) # the negative weighted Bernoulli logistic log-likelihood
        return (loss + 0.5 * np.dot(precision * beta, beta)) / scale

    def gradient(beta):
        # This form avoids subtracting a probability very close to one.
        residual = sign * expit(sign * (x @ beta))
        return (x.T @ (weights * residual) + precision * beta) / scale

    def hessian(beta):
        eta = x @ beta
        curvature = weights * expit(eta) * expit(-eta)
        return (x.T @ (curvature[:, None] * x) + np.diag(precision)) / scale

    result = minimize(objective, start, method="trust-exact", jac=gradient,
                      hess=hessian, options={"gtol": 1e-9, "maxiter": 100})
    beta = result.x
    if np.max(np.abs(gradient(beta))) > 1e-9:
        # Near the maximum, objective differences may round to zero. Solve the
        # score equation with SciPy instead; strict convexity gives one root.
        beta = root(gradient, beta, jac=hessian, options={"xtol": 1e-10}).x
    error = float(np.max(np.abs(gradient(beta))))
    if not np.isfinite(beta).all() or not np.isfinite(error) or error > 1e-9:
        raise RuntimeError(f"MAP gradient {error:.3g}: {result.message}")
    return beta, error


def fit_att(x, a, y, weights, start):
    """Fit propensity, then the control outcome, then the weighted ATT ratio."""
    control, treated = a == 0, a == 1
    if weights[control].sum() <= 0 or weights[treated].sum() <= 0:
        raise RuntimeError("Both treatment arms need positive total weight")
    try:
        alpha, error1 = logistic_map(x, a, weights, [10] + [5] * 9, start[0]) # fits the logistic treatment-assignment model
    except RuntimeError as error:
        raise RuntimeError(f"propensity module: {error}") from error
    eta = x @ alpha # calculates the fitted propensity linear predictors
    e = expit(eta) # converts them into estimated propensity scores
    basis = np.column_stack((np.ones(len(a)), e, e**2, e**3)) # constructs the outcome-model design matrix
    try:
        beta, error2 = logistic_map(basis[control], y[control], weights[control], # outcome model is fitted using controls only
                                    np.full(4, 2.5), start[1])
    except RuntimeError as error:
        raise RuntimeError(f"control outcome module: {error}") from error
    residual = y - expit(basis @ beta) # calculates the residuals for every observation, including treated observations
    # exp(eta) equals e/(1-e), avoiding division by an e rounded to one.
    active = control & (weights > 0)
    with np.errstate(over="ignore", invalid="ignore"): # temporarily suppress numpy warnings - we'll check for them later so that we can raise a clear error
        odds = np.exp(eta[active]) # calculates e_i/(1-e_i) for the active controls
        numerator = weights[treated] @ residual[treated] # first part of the numerator, i.e. weighted sum of treated residuals
        numerator -= (weights[active] * odds) @ residual[active] # second part of the numerator, i.e. negative propensity-odds-weighted control residuals
        estimate = numerator / weights[treated].sum() # divides by the total treated weight
    if not np.isfinite(odds).all() or not np.isfinite(estimate):
        raise RuntimeError("Nonfinite ATT augmentation; no clipping is applied")
    return float(estimate), (alpha, beta), max(error1, error2)


def fit_baselines(x, a, y, dataset_id, root_seed, draws=1000):
    """Return the shared MAP point, PBMI/bootstrap draws and the worst score. This function is called by run.py"""
    x, a, y = np.asarray(x, float), np.asarray(a, float), np.asarray(y, float)
    if a.ndim != 1 or x.shape != (len(a), 9) or y.shape != a.shape:
        raise ValueError("Use nine covariates and matching treatment/outcome vectors")
    if not np.isfinite(x).all() or not np.isin(a, [0, 1]).all() or not np.isin(y, [0, 1]).all():
        raise ValueError("Covariates must be finite; treatment and outcome must be binary")
    if isinstance(draws, bool) or int(draws) != draws or draws < 2:
        raise ValueError("draws must be an integer >= 2")
    draws = int(draws)
    n = len(a)
    x = np.column_stack((np.ones(n), x)) # adds an intercept column
    try:
        point, start, worst = fit_att(x, a, y, np.ones(n), (np.zeros(10), np.zeros(4))) # obtain a two-step MAP estimate of the model coefficients, then insert those into ATT formula to estimate ATT
    except RuntimeError as error:
        raise RuntimeError(f"Dataset {dataset_id}, original fit: {error}") from error
    samples = {}
    for method, stage in (("pbmi", 1), ("bootstrap", 2)): # stage distinguishes random number streams
        rng = np.random.default_rng(np.random.SeedSequence([root_seed, dataset_id, stage]))
        samples[method] = np.empty(draws)
        for b in range(draws):
            # RAW exponential weights preserve the prior/likelihood balance. We should not rescale the weights.
            weights = (rng.exponential(size=n) if method == "pbmi" else
                       np.bincount(rng.integers(n, size=n), minlength=n).astype(float)) # Selects exponential weights for PBMI. For bootstrap samples n row indices uniformly with replacement from 0 to n-1 and counts how often each observation was selected
            try:
                # Always start at the original MAP; reuse the same weights in both fits and ATT.
                samples[method][b], _, error = fit_att(x, a, y, weights, start) 
            except RuntimeError as error:
                raise RuntimeError(f"Dataset {dataset_id}, {method} draw {b+1}: {error}") from error
            worst = max(worst, error) # largest optimisation gradient diagnostic across the original fits and all PBMI and bootstrap fits.
    return {"point_estimate": point, "draws": samples, "max_scaled_gradient": worst}
