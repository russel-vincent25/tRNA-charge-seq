# The background null: truncation bias and failed fits

Branch `fix/background-null`. Covers two statistical defects in
`trnaseq/modifications/modification_caller.py`.

The background model is a per-channel, per-sample beta-binomial null
(`ChannelBackground`). Site calling is an OR across three channels
(`mismatch`, `deletion`, `rt_stop`), each tested against its own null.
Nothing here changes detection thresholds, the per-enzyme channel priors, or
the OR logic.

---

## Defect 1 — the null was fitted on a truncated sample with no correction

`estimate_channel_backgrounds` keeps only positions whose observed rate is
below `0.5 * threshold` (the "null bulk"), so real modifications do not
inflate the null. That guard is right. But the retained sample was then
handed to an ordinary beta-binomial likelihood, as if it were complete.
Truncating the upper tail biases both `mu` and `rho` **downward**, which makes
downstream p-values **anti-conservative**.

### The fix

`_fit_beta_binomial(k, n, truncation=...)` now optionally fits the correctly
normalised truncated likelihood. Each observation's truncation point depends
on its own coverage:

```
c_i       = floor(truncation * n_i)
loglik_i  = betabinom.logpmf(k_i, n_i, a, b) - log betabinom.cdf(c_i, n_i, a, b)
```

Degenerate cases are guarded: `c_i < 0` and `c_i < k_i` (no usable
normalisation) are dropped; a CDF at or near 0 is clipped to `1e-300` so the
optimiser sees a large finite penalty rather than a NaN it cannot escape; a
non-finite objective returns `1e12`. The CDF is evaluated once per *distinct*
coverage and broadcast back, which is what makes the fit affordable at all.

### Does it matter? — simulation

`docs/simulations/background_null_truncation_sim.py`. Method: draw
`k_i ~ BetaBinom(n_i, a, b)` at known `(mu, rho)`, coverage log-uniform over
50–5,000, 2,000 positions per dataset; truncate exactly as the pipeline does
(retain rate `< 0.05`); fit the naive and the truncated likelihood **on the
same retained sample**; repeat over 30 independent datasets. Nelder-Mead at
production settings (`xatol 1e-4`, `fatol 1e-3`, `maxiter 400`). Uncertainty
is the 95% CI of the mean over the 30 datasets. Because both fits see the
same data, the paired difference is the sensitive statistic.

Reproduce with:

```
python docs/simulations/background_null_truncation_sim.py 30
```

**At the measured null rates** (h-tRNAseq atlas, stage 6):

| channel  | param | true     | naive bias | truncated bias | paired difference (truncated − naive) |
|----------|-------|----------|-----------:|---------------:|--------------------------------------:|
| mismatch | mu    | 8.01e-04 |     +0.10% |         +0.10% |                    +0.005% ± 0.001%   |
| mismatch | rho   | 1.10e-03 |     +0.64% |         +0.66% |                    +0.015% ± 0.001%   |
| deletion | mu    | 1.32e-04 |     +2.73% |         +3.18% |                    +0.392% ± 0.164%   |
| deletion | rho   | 5.35e-03 |     +6.14% |         +6.99% |                    +0.604% ± 0.242%   |
| rt_stop  | mu    | 6.58e-04 |     −0.21% |         +0.33% |                    +0.519% ± 0.074%   |
| rt_stop  | rho   | 6.68e-03 |     −1.20% |         −0.34% |                    +0.834% ± 0.115%   |

The cut retains 99.98–100.00% of positions. The analytic mass above the cut
is 4.8e-06 (mismatch), 1.0e-05 (deletion) and 9.9e-05 (rt_stop) — the null
sits 20–40 sigma below the truncation point, so there is essentially nothing
to correct for.

**Stress regime**, where the cut genuinely removes mass, confirming the
mechanism is real and the implementation is correct:

| setting            | param | true  | naive bias | truncated bias | mass above cut |
|--------------------|-------|-------|-----------:|---------------:|---------------:|
| mu .02 / rho .05   | mu    | 0.020 |    −50.59% |         −2.03% |          12.6% |
| mu .02 / rho .05   | rho   | 0.050 |    −60.09% |         −1.74% |                |
| mu .03 / rho .10   | mu    | 0.030 |    −70.99% |        −10.80% |          19.1% |
| mu .03 / rho .10   | rho   | 0.100 |    −76.28% |        −10.98% |                |

### Verdict — correction implemented, **off by default**

The diagnosis is correct in principle: the naive fit on truncated data is
biased down by 50–76% once the cut removes real mass, exactly the
anti-conservative direction. At *this pipeline's* null rates it is not:
the correction moves `mu` by at most +0.52% and `rho` by at most +0.83%,
and always **upward** (the conservative direction). Those shifts are an
order of magnitude smaller than the fit's own replicate-to-replicate
scatter (per-dataset SD 10–35%), and they are far too small to move a BH
q-value across any plausible alpha.

So `truncation_correction=False` is the default and production behaviour is
unchanged. The correction costs roughly 100x the fitting time (a CDF
evaluation per distinct coverage per optimiser step) for a sub-1% shift.

**Turn it on** (`estimate_channel_backgrounds(..., truncation_correction=True)`)
if the calling thresholds are lowered, if a channel's fitted `mu` rises
towards the cut, or if a new RT/chemistry pushes `rho` into the 1e-2 range.
A practical trigger: the correction starts to matter once the analytic mass
above the cut exceeds about 1e-2. Note also that at `rho = 0.1` even the
truncated fit is 11% low at `maxiter 400`; that regime needs a tighter
optimiser tolerance as well.

---

## Defect 2 — a failed fit silently degraded to a binomial

`_fit_beta_binomial` returned `(mu0, 0.0)` when Nelder-Mead did not yield a
finite objective. `rho = 0` **is** the binomial: it understates the variance
of a per-position null rate and overstates significance. Nothing recorded
that it had happened. The convergence check was also only
`np.isfinite(res.fun)` — `res.success` was never inspected.

### What changed

1. **Convergence is checked properly**: `res.success and np.isfinite(res.fun)
   and res.fun < 1e12`. A Nelder-Mead failure with a finite objective — the
   common case, and the exact hole in the old code — is now caught.

2. **`ChannelBackground` records the fit**, alongside the existing `source`
   field (which says *which positions* were used and is unchanged):

   - `converged: bool`
   - `fit_status: str` — `'ok'`, `'nonconvergence'`, `'no_events'`,
     `'default'`, `'pooled_rho'`

3. **Warnings** go through the module logger, naming sample and channel.

4. **Status is propagated into the outputs**:

   - `modification_calls` (per call, for the dominant channel):
     `background_rho`, `background_fit_status`, `background_converged`
   - `modification_summary.csv` (per sample, per channel):
     `bg_fit_status_{channel}`, `bg_converged_{channel}`
   - stage 6 logs a `WARNING` line counting degraded samples per channel

   Filter on `background_converged` before quoting a call's significance.

### The fallback: not `rho = 0`

`rho = 0` is the one value that is certainly wrong — it asserts no
overdispersion at all, in a model whose entire purpose is that per-position
null rates *are* overdispersed. Two safer fallbacks are implemented, in
order:

**1. Method of moments, within the sample** (`_moment_rho`). Under the common
correlation model `Var(k_i) = n_i mu (1 - mu) [1 + (n_i - 1) rho]`, the
Pearson chi-square statistic `S = sum (k_i - n_i mu)^2 / (n_i mu (1 - mu))`
has expectation `N + rho * sum(n_i - 1)`, giving
`rho_hat = (S - N) / sum(n_i - 1)`, clipped to `[0, 0.5]`. This is the
standard ANOVA/quasi-likelihood moment estimator for clustered binomial data
with unequal cluster sizes (Kleinman, *JASA* 1973). It is closed-form, always
available, needs no other sample, and cannot silently produce a binomial.
Being computed on the same truncated sample it is itself mildly
conservative-leaning low, but it is bounded away from 0 whenever the data are
overdispersed at all.

**2. Pooled `rho` across samples** (`pool_channel_backgrounds`). Samples in a
run share library chemistry, RT and batch, so the **median `rho` of that
channel's converged fits** is a better estimate than any single-sample
rescue. Stage 6 applies this after fitting every sample. Only `rho` is
pooled: `mean` is genuinely sample-specific (depth and batch move it) and a
failed fit still pins it down through the pooled event rate. A pooled null is
re-flagged `fit_status='pooled_rho'` and **keeps `converged=False`**, so it
stays visibly degraded downstream. If no sample's fit converged for a
channel, the method-of-moments value stands and the warning says so.

`'no_events'` gets `rho = 0` **from the fit** — an all-zero sample carries no
information about overdispersion, and the `0.5 / total_n` half-event
pseudocount is already a conservative floor on the mean. It is flagged
`converged=False`, so pooling then treats it like any other degraded null and
lends it the channel's pooled `rho` when a converged sibling sample exists.
That only ever makes it more conservative.

---

## Tests

`tests/test_background_null.py` (16 tests).

- **Defect 1**: the truncated fit recovers known `(mu, rho)` in the stress
  regime where the naive fit is >50% low (averaged over four datasets, since
  a single fit has ~15% scatter); the correction is inert at each measured
  channel's null rate; degenerate truncation points do not crash; the flag
  reaches `estimate_channel_backgrounds`.
- **Defect 2**: a `res.success == False` with a *finite* objective is flagged
  rather than silently returned; a non-finite objective likewise; the
  fallback `rho` is `> 0`; `ChannelBackground` carries the status and the
  warning names the sample and channel; `'ok'`, `'no_events'` and
  `'default'` are each distinguished; `_moment_rho` detects overdispersion
  and returns ~0 on binomial data; pooling replaces degraded `rho` with the
  converged median, does not mutate its input, leaves converged fits alone,
  and warns when there is nothing to pool from; the status reaches
  `finalize_calls` output columns.

Full suite: 307 passed, 20 skipped.
