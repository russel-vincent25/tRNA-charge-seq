"""Does fitting a TRUNCATED beta-binomial sample with a NAIVE (untruncated)
likelihood materially bias (mu, rho)?

Mirrors trnaseq.modifications.modification_caller.estimate_channel_backgrounds:
  draw k_i ~ BetaBinom(n_i, a, b);  retain k_i/n_i < 0.5 * threshold (= 0.05);
  fit naive BB MLE vs truncated BB MLE on the retained sample.
Uncertainty = 95% CI of the mean over independent replicate datasets.
"""
import sys, time
import numpy as np
from scipy import stats
from scipy.optimize import minimize
from scipy.special import expit, logit

THRESH = 0.10               # DERIVATION_THRESHOLDS (all three channels)
TRUNC_FRAC = 0.5 * THRESH   # 0.05, exactly the pipeline's "null bulk" cut
MAXITER = 400               # production setting

def ab(mu, rho):
    return mu * (1 - rho) / rho, (1 - mu) * (1 - rho) / rho

def _logcdf_trunc(c, n, a, b, uniq):
    """log P(X <= c_i) with the cdf evaluated once per distinct coverage."""
    un, inv = uniq
    vals = stats.betabinom.cdf(np.floor(TRUNC_FRAC * un).astype(np.int64), un, a, b)
    return np.log(np.clip(vals, 1e-300, 1.0))[inv]

def fit(k, n, truncated, maxiter=MAXITER):
    total_k, total_n = float(k.sum()), float(n.sum())
    if total_k == 0:
        return 0.5 / total_n, 0.0
    mu0 = total_k / total_n
    uniq = None
    if truncated:
        un, inv = np.unique(n, return_inverse=True)
        uniq = (un, inv)
        c = np.floor(TRUNC_FRAC * n).astype(np.int64)
    def nll(t):
        a, b = ab(expit(t[0]), expit(t[1]))
        lp = stats.betabinom.logpmf(k, n, a, b)
        if truncated:
            lp = lp - _logcdf_trunc(c, n, a, b, uniq)
        v = -lp.sum()
        return v if np.isfinite(v) else 1e12
    res = minimize(nll, [logit(mu0), logit(1e-3)], method='Nelder-Mead',
                   options={'xatol': 1e-4, 'fatol': 1e-3, 'maxiter': maxiter})
    if not np.isfinite(res.fun):
        return mu0, 0.0
    return float(expit(res.x[0])), float(expit(res.x[1]))

def run(mu, rho, n_pos, n_rep, cov_hi, seed=0):
    rng = np.random.default_rng(seed)
    a, b = ab(mu, rho)
    acc = {k: [] for k in ('nm', 'nr', 'tm', 'tr', 'kept', 'tail')}
    for _ in range(n_rep):
        n = np.rint(10 ** rng.uniform(np.log10(50), np.log10(cov_hi), n_pos)).astype(np.int64)
        p = rng.beta(a, b, size=n_pos)
        k = rng.binomial(n, p)
        keep = (k / n) < TRUNC_FRAC
        acc['kept'].append(keep.mean())
        cc = np.floor(TRUNC_FRAC * n).astype(np.int64)
        acc['tail'].append(float(np.mean(stats.betabinom.sf(cc, n, a, b))))
        kk, nn = k[keep], n[keep]
        if kk.sum() == 0:
            continue
        m1, r1 = fit(kk, nn, truncated=False)
        m2, r2 = fit(kk, nn, truncated=True)
        acc['nm'].append(m1); acc['nr'].append(r1)
        acc['tm'].append(m2); acc['tr'].append(r2)
    return {k: np.asarray(v) for k, v in acc.items()}

def ci(x):
    return x.mean(), 1.96 * x.std(ddof=1) / np.sqrt(len(x))

REALISTIC = [('mismatch', 8.01e-04, 1.10e-03),
             ('deletion', 1.32e-04, 5.35e-03),
             ('rt_stop',  6.58e-04, 6.68e-03)]
STRESS = [('stress-A', 2.00e-02, 5.00e-02),
          ('stress-B', 3.00e-02, 1.00e-01)]

def report(settings, n_pos, n_rep, cov_hi):
    hdr = (f"{'setting':<10}{'param':<5}{'true':>11}{'naive fit':>22}{'bias':>9}"
           f"{'truncated fit':>22}{'bias':>9}")
    print(hdr); print('-' * len(hdr))
    for name, mu, rho in settings:
        t0 = time.time()
        r = run(mu, rho, n_pos, n_rep, cov_hi)
        for par, true, kn, kt in (('mu', mu, 'nm', 'tm'), ('rho', rho, 'nr', 'tr')):
            nm, nc = ci(r[kn]); tm, tc = ci(r[kt])
            print(f"{name:<10}{par:<5}{true:>11.3e}{nm:>14.4e} +/-{nc:<6.1e}"
                  f"{100*(nm-true)/true:>8.2f}%{tm:>14.4e} +/-{tc:<6.1e}"
                  f"{100*(tm-true)/true:>8.2f}%")
        for par, kn, kt in (('mu', 'nm', 'tm'), ('rho', 'nr', 'tr')):
            d = 100 * (r[kt] - r[kn]) / r[kn]          # paired, same dataset
            dm, dc = ci(d)
            print(f"{'':<10}{par:<5}paired truncated-vs-naive difference: "
                  f"{dm:+.3f}% +/- {dc:.3f}% (95% CI of the mean, n={len(d)})")
        print(f"{'':<10}retained {100*r['kept'].mean():.4f}% of positions; "
              f"mean analytic mass above the cut = {r['tail'].mean():.3e}  "
              f"[{time.time()-t0:.0f}s, {n_rep} replicate datasets]\n")

if __name__ == '__main__':
    n_rep = int(sys.argv[1]) if len(sys.argv) > 1 else 15
    print(f"Truncation at observed rate < {TRUNC_FRAC}; Nelder-Mead maxiter={MAXITER} "
          f"(production settings)\n")
    print("=== REALISTIC: measured channel nulls, coverage log-uniform 50-5,000 ===")
    report(REALISTIC, n_pos=2000, n_rep=n_rep, cov_hi=5000)
    print("=== STRESS: mu/rho large enough that the cut removes real mass ===")
    report(STRESS, n_pos=2000, n_rep=n_rep, cov_hi=5000)
