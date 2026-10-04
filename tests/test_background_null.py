"""Regression tests for the per-channel background null model.

Two defects are covered:

1. The null is fitted on an upper-truncated sample (only positions below half
   the calling threshold are retained). Fitting that with an ordinary
   beta-binomial likelihood biases mu and rho **downward**, which makes
   downstream p-values anti-conservative. ``_fit_beta_binomial(...,
   truncation=...)`` fits the correctly normalised truncated likelihood.

   The simulation in ``docs/BACKGROUND_NULL_TRUNCATION.md`` shows the bias is
   numerically inert at this pipeline's measured null rates, so the
   correction is OFF by default. Both facts are pinned here: the correction
   works where truncation bites, and it changes nothing where it does not.

2. A failed fit used to return ``rho = 0`` -- a plain binomial -- silently.
   Non-convergence must now be detected (``res.success`` as well as
   finiteness), flagged on the :class:`ChannelBackground`, warned about, and
   backed off to a method-of-moments rho rather than 0.
"""

import logging

import numpy as np
import pandas as pd
import pytest
from scipy import stats
from scipy.optimize import OptimizeResult

from trnaseq.modifications import modification_caller as mc
from trnaseq.modifications.modification_caller import (
    CHANNELS,
    ChannelBackground,
    ModificationCaller,
    _fit_beta_binomial,
    _moment_rho,
    estimate_channel_backgrounds,
    pool_channel_backgrounds,
)

# The pipeline keeps positions with rate < 0.5 * threshold; thresholds are 0.10
TRUNC = 0.05

# Measured nulls of the three channels (h-tRNAseq atlas, stage 6)
MEASURED = {
    'mismatch': (8.01e-04, 1.10e-03),
    'deletion': (1.32e-04, 5.35e-03),
    'rt_stop': (6.58e-04, 6.68e-03),
}


def _draw(mu, rho, n_pos, seed, cov_lo=50, cov_hi=5000):
    """A beta-binomial sample truncated exactly as the pipeline truncates."""
    rng = np.random.default_rng(seed)
    a, b = mu * (1 - rho) / rho, (1 - mu) * (1 - rho) / rho
    n = np.rint(10 ** rng.uniform(np.log10(cov_lo), np.log10(cov_hi), n_pos)).astype(np.int64)
    k = rng.binomial(n, rng.beta(a, b, size=n_pos))
    keep = (k / n) < TRUNC
    return k[keep], n[keep]


# ---------------------------------------------------------------------------
# Defect 1 -- truncated likelihood
# ---------------------------------------------------------------------------

class TestTruncatedLikelihood:

    def test_truncated_fit_recovers_params_where_naive_fails(self):
        """Where the cut removes real mass, only the truncated fit is unbiased.

        mu = 0.02, rho = 0.05 leaves ~13% of the mass above the cut. The naive
        likelihood then under-estimates BOTH parameters by >50% -- the
        anti-conservative direction -- while the truncated likelihood recovers
        them. Averaged over four independent datasets because a single fit at
        this sample size has ~15% replicate scatter (see
        docs/BACKGROUND_NULL_TRUNCATION.md).
        """
        mu, rho = 0.02, 0.05
        naive, trunc = [], []
        for seed in (7, 8, 9, 10):
            k, n = _draw(mu, rho, n_pos=2000, seed=seed)
            assert len(k) < 2000, "this regime must actually truncate"
            nm, nr, ns = _fit_beta_binomial(k, n)
            tm, tr, ts = _fit_beta_binomial(k, n, truncation=TRUNC)
            assert ns == 'ok' and ts == 'ok'
            naive.append((nm, nr))
            trunc.append((tm, tr))
        naive_mu, naive_rho = np.mean(naive, axis=0)
        trunc_mu, trunc_rho = np.mean(trunc, axis=0)

        # Naive is biased DOWN -- understates the null, overstates significance
        assert naive_mu < 0.6 * mu
        assert naive_rho < 0.6 * rho
        # Truncated recovers both
        assert trunc_mu == pytest.approx(mu, rel=0.25)
        assert trunc_rho == pytest.approx(rho, rel=0.30)
        # and is unambiguously the better estimate of each
        assert abs(trunc_mu - mu) < abs(naive_mu - mu)
        assert abs(trunc_rho - rho) < abs(naive_rho - rho)

    @pytest.mark.parametrize('channel', sorted(MEASURED))
    def test_correction_is_inert_at_measured_null_rates(self, channel):
        """Justifies ``truncation_correction=False`` being the default.

        At the measured (mu, rho) the cut sits far into the tail (<1e-4 of
        the mass), so the two likelihoods must agree to well within the
        Monte-Carlo noise of the fit itself.
        """
        mu, rho = MEASURED[channel]
        k, n = _draw(mu, rho, n_pos=2000, seed=11)
        naive_mu, naive_rho, _ = _fit_beta_binomial(k, n)
        trunc_mu, trunc_rho, _ = _fit_beta_binomial(k, n, truncation=TRUNC)
        assert trunc_mu == pytest.approx(naive_mu, rel=0.02)
        assert trunc_rho == pytest.approx(naive_rho, rel=0.05)
        # and the correction never moves the null DOWN (anti-conservative)
        assert trunc_mu >= naive_mu * 0.99

    def test_degenerate_truncation_points_do_not_crash(self):
        """c_i < 0, c_i below the observed count, and tiny n_i are survivable."""
        k = np.array([0, 1, 0, 2, 0], dtype=np.int64)
        n = np.array([3, 5, 10, 8, 4000], dtype=np.int64)   # floor(.05*n) = 0,0,0,0,200
        mu, rho, status = _fit_beta_binomial(k, n, truncation=TRUNC)
        assert np.isfinite(mu) and np.isfinite(rho)
        assert 0.0 <= rho < 1.0 and 0.0 < mu < 1.0
        assert status in ('ok', 'nonconvergence', 'no_events')

    def test_truncation_flag_reaches_estimate_channel_backgrounds(self):
        rng = np.random.default_rng(3)
        pscm, ref = {}, {}
        for i in range(20):
            mat = np.zeros((76, 8))
            mat[:, 2] = mat[:, 6] = 2000
            mat[:, 7] = rng.binomial(2000, 1e-3, 76)
            pscm[f't{i}'] = mat
            ref[f't{i}'] = {'seq': 'G' * 76}
        plain = estimate_channel_backgrounds(pscm, ref)['rt_stop']
        corrected = estimate_channel_backgrounds(
            pscm, ref, truncation_correction=True)['rt_stop']
        assert plain.mean == pytest.approx(corrected.mean, rel=0.05)


# ---------------------------------------------------------------------------
# Defect 2 -- failed fits must be flagged, never silently binomial
# ---------------------------------------------------------------------------

class TestFitStatus:

    @staticmethod
    def _force_failure(monkeypatch, success=False, fun=1.0):
        def fake_minimize(f, x0, **kw):
            return OptimizeResult(x=np.asarray(x0, dtype=float), fun=fun,
                                  success=success, nit=0, status=2,
                                  message='forced')
        monkeypatch.setattr(mc, 'minimize', fake_minimize)

    def test_nonconvergence_is_flagged_not_silently_returned(self, monkeypatch):
        """``res.success == False`` must be caught even when ``res.fun`` is finite.

        This is the exact hole in the old code: it checked only
        ``np.isfinite(res.fun)``, so a Nelder-Mead failure returned
        ``(mu0, 0.0)`` -- a binomial null -- with nothing recorded.
        """
        self._force_failure(monkeypatch, success=False, fun=1.0)
        k, n = _draw(1e-3, 1e-2, n_pos=800, seed=5)
        mu, rho, status = _fit_beta_binomial(k, n)
        assert status == 'nonconvergence'
        assert rho > 0, "fallback must not be a plain binomial (rho == 0)"
        assert mu == pytest.approx(k.sum() / n.sum(), rel=1e-9)

    def test_nonfinite_objective_also_flagged(self, monkeypatch):
        self._force_failure(monkeypatch, success=True, fun=np.inf)
        k, n = _draw(1e-3, 1e-2, n_pos=800, seed=5)
        _, rho, status = _fit_beta_binomial(k, n)
        assert status == 'nonconvergence' and rho > 0

    def test_background_carries_status_and_warns(self, monkeypatch, caplog):
        self._force_failure(monkeypatch, success=False, fun=1.0)
        rng = np.random.default_rng(2)
        pscm, ref = {}, {}
        for i in range(15):
            mat = np.zeros((76, 8))
            mat[:, 2] = mat[:, 6] = 1500
            mat[:, 7] = rng.binomial(1500, 2e-3, 76)
            pscm[f't{i}'] = mat
            ref[f't{i}'] = {'seq': 'G' * 76}
        with caplog.at_level(logging.WARNING, logger=mc.__name__):
            bgs = estimate_channel_backgrounds(pscm, ref, sample_label='S512_R1')
        bg = bgs['rt_stop']
        assert bg.fit_status == 'nonconvergence'
        assert bg.converged is False
        assert bg.rho > 0
        text = caplog.text
        assert 'S512_R1' in text and 'rt_stop' in text and 'did not converge' in text

    def test_converged_fit_is_marked_ok(self):
        rng = np.random.default_rng(1)
        pscm, ref = {}, {}
        mu, rho = 1e-3, 1e-2
        a, b = mu * (1 - rho) / rho, (1 - mu) * (1 - rho) / rho
        for i in range(60):
            mat = np.zeros((76, 8))
            mat[:, 2] = mat[:, 6] = 3000
            mat[:, 7] = rng.binomial(3000, rng.beta(a, b, 76))
            pscm[f't{i}'] = mat
            ref[f't{i}'] = {'seq': 'G' * 76}
        bg = estimate_channel_backgrounds(pscm, ref)['rt_stop']
        assert bg.fit_status == 'ok' and bg.converged is True

    def test_no_events_is_its_own_status(self):
        mat = np.zeros((20, 8))
        mat[:, 2] = mat[:, 6] = 1000
        bgs = estimate_channel_backgrounds(
            {'Synthetic_a': mat}, {'Synthetic_a': {'seq': 'G' * 20}})
        assert bgs['deletion'].fit_status == 'no_events'
        assert bgs['deletion'].converged is False
        assert bgs['deletion'].source == 'synthetic'   # source is unchanged

    def test_default_background_is_flagged(self):
        bgs = estimate_channel_backgrounds({}, {})
        for c in CHANNELS:
            assert bgs[c].fit_status == 'default'
            assert bgs[c].converged is False

    def test_moment_rho_detects_overdispersion(self):
        k, n = _draw(1e-3, 2e-2, n_pos=4000, seed=9)
        mu = k.sum() / n.sum()
        assert _moment_rho(k, n, mu) > 5e-3
        # and is ~0 on a genuinely binomial sample
        rng = np.random.default_rng(4)
        nb = np.full(4000, 2000, dtype=np.int64)
        kb = rng.binomial(2000, 1e-3, 4000)
        assert _moment_rho(kb, nb, kb.sum() / nb.sum()) < 2e-3


class TestPooling:

    def test_pooled_rho_replaces_degraded_fits(self, caplog):
        samples = {
            'A': {'mismatch': ChannelBackground(mean=8e-4, rho=1.0e-3, converged=True),
                  'deletion': ChannelBackground(mean=1e-4, rho=5.0e-3, converged=True),
                  'rt_stop': ChannelBackground(mean=6e-4, rho=6.0e-3, converged=True)},
            'B': {'mismatch': ChannelBackground(mean=9e-4, rho=1.2e-3, converged=True),
                  'deletion': ChannelBackground(mean=1e-4, rho=7.0e-3, converged=True),
                  'rt_stop': ChannelBackground(mean=7e-4, rho=8.0e-3, converged=True)},
            'C': {'mismatch': ChannelBackground(mean=8e-4, rho=0.0, converged=False,
                                                fit_status='nonconvergence'),
                  'deletion': ChannelBackground(mean=1e-4, rho=6.0e-3, converged=True),
                  'rt_stop': ChannelBackground(mean=6e-4, rho=7.0e-3, converged=True)},
        }
        with caplog.at_level(logging.WARNING, logger=mc.__name__):
            out = pool_channel_backgrounds(samples)
        fixed = out['C']['mismatch']
        assert fixed.rho == pytest.approx(np.median([1.0e-3, 1.2e-3]))
        assert fixed.fit_status == 'pooled_rho'
        assert fixed.converged is False       # still visibly degraded
        assert fixed.mean == 8e-4             # mean stays sample-specific
        assert 'pooled' in caplog.text
        # converged fits are untouched, and the input is not mutated
        assert out['A']['mismatch'].rho == 1.0e-3
        assert samples['C']['mismatch'].rho == 0.0

    def test_no_converged_fit_to_borrow_from(self, caplog):
        samples = {'A': {'mismatch': ChannelBackground(
            mean=8e-4, rho=0.0, converged=False, fit_status='nonconvergence')}}
        with caplog.at_level(logging.WARNING, logger=mc.__name__):
            out = pool_channel_backgrounds(samples)
        assert out['A']['mismatch'].fit_status == 'nonconvergence'
        assert 'no converged fit' in caplog.text


class TestStatusReachesOutputs:

    def test_finalize_calls_carries_background_provenance(self):
        bgs = {
            'mismatch': ChannelBackground(mean=8e-4, rho=1e-3,
                                          converged=True, fit_status='ok'),
            'deletion': ChannelBackground(mean=1e-4, rho=5e-3, converged=False,
                                          fit_status='pooled_rho'),
            'rt_stop': ChannelBackground(mean=1e-3, rho=1e-2,
                                         converged=True, fit_status='ok'),
        }
        caller = ModificationCaller(organism='Escherichia coli',
                                    channel_backgrounds=bgs)
        rows = [{'position': p, 'coverage': 2000, 'correct_nt': 'G',
                 'mismatch_rate': 8e-4, 'gap_rate': 1e-4, 'rt_stop_pct': 0.1}
                for p in range(2, 42)]
        rows.append({'position': 100, 'coverage': 2000, 'correct_nt': 'G',
                     'mismatch_rate': 8e-4, 'gap_rate': 0.30, 'rt_stop_pct': 0.1})
        calls = caller.call_all('tRNA-x', pd.DataFrame(rows))

        assert not calls.empty
        for col in ('background_fit_status', 'background_converged',
                    'background_rho', 'background_error_rate'):
            assert col in calls.columns
        row = calls.iloc[0]
        assert row['dominant_channel'] == 'deletion'
        assert row['background_fit_status'] == 'pooled_rho'
        assert bool(row['background_converged']) is False
        assert row['background_rho'] == 5e-3
