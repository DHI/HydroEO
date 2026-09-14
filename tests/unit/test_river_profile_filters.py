"""Unit tests for river-profile filter stages (synthetic arrays, no network)."""

import numpy as np
import pytest

from HydroEO.utils.filters import river_profile_filters as rpf

pytestmark = pytest.mark.unit


def test_robust_line_ignores_outlier():
    x = np.linspace(0, 1000, 50)
    y_true = 2.0 * x + 5.0
    y = y_true.copy()
    y[10] += 500.0  # single large outlier

    a, b, r = rpf.robust_line(x, y, huber_k=1.5, iters=3)

    assert a == pytest.approx(2.0, abs=0.05)
    assert b == pytest.approx(5.0, abs=5.0)
    assert abs(r[10]) > abs(r[0]) * 10


def test_sliding_window_matches_brute_force():
    rng = np.random.default_rng(0)
    xv = np.sort(rng.uniform(0, 1000, 30))
    half_width = 75.0

    for j, j0, j1 in rpf.sliding_window(xv, half_width):
        expected = np.where(np.abs(xv - xv[j]) <= half_width)[0]
        assert j0 == expected.min()
        assert j1 == expected.max()


def test_hampel_1d_meter_masks_spike():
    x = np.arange(0, 40) * 1000.0
    y = np.full(x.shape, 10.0, dtype=np.float32)
    y[20] = 25.0  # single spike

    out, changed = rpf.hampel_1d_meter(
        x, y, win_m=10_000.0, k=3.0, action="mask", min_valid=8, replace_mode="trend"
    )

    assert changed[20]
    assert np.isnan(out[20])
    assert not changed[0]
    assert out[0] == pytest.approx(10.0)


def test_hampel_1d_meter_too_few_points_returns_unchanged():
    x = np.arange(0, 3) * 1000.0
    y = np.array([1.0, 2.0, 3.0], dtype=np.float32)

    out, changed = rpf.hampel_1d_meter(
        x, y, win_m=10_000.0, k=3.0, action="mask", min_valid=8, replace_mode="trend"
    )

    np.testing.assert_array_equal(out, y)
    assert not changed.any()


def test_rolling_quantile_trend_follows_linear_signal():
    x = np.arange(0, 100) * 1000.0
    y = (0.001 * x + 3.0).astype(np.float32)

    out = rpf.rolling_quantile_trend(x, y, win_m=10_000.0, q=0.5, min_valid=8)

    finite = np.isfinite(out)
    assert finite.sum() > 0
    np.testing.assert_allclose(out[finite], y[finite], atol=0.2)


def test_rolling_density_counts_neighbors():
    x = np.arange(0, 10) * 1000.0
    y = np.ones_like(x)

    counts = rpf.rolling_density(x, y, total_win_m=2_000.0)

    # +/-1000m window around an interior point includes itself + one neighbor each side
    assert counts[5] == 3
    # edge point only has neighbors on one side
    assert counts[0] == 2


def test_rolling_density_ignores_nan_values():
    x = np.arange(0, 10) * 1000.0
    y = np.ones_like(x)
    y[4] = np.nan

    counts = rpf.rolling_density(x, y, total_win_m=2_000.0)

    # density is a property of the neighborhood, not conditioned on the
    # point's own value: window around x=4000 has 2 finite neighbors (3000,
    # 5000) but not itself (NaN)
    assert counts[4] == 2
    # window around x=3000 (2000, 3000, 4000) also loses the NaN at 4000
    assert counts[3] == 2
    # an untouched interior point still sees all 3 neighbors as finite
    assert counts[7] == 3


def test_dilate_mask_grows_true_regions():
    mask = np.array([False, False, True, False, False])

    out = rpf.dilate_mask(mask, k=1)

    np.testing.assert_array_equal(out, [False, True, True, True, False])


def test_dilate_mask_noop_when_k_zero():
    mask = np.array([False, True, False])
    out = rpf.dilate_mask(mask, k=0)
    np.testing.assert_array_equal(out, mask)


def test_lsq_spline_fill_only_touches_gaps():
    x = np.linspace(0, 100_000, 60)
    y = (5.0 + 0.0002 * x).astype(np.float32)
    y_in = y.copy()
    gap = slice(25, 30)
    y_in[gap] = np.nan

    y_out, y_fit = rpf.lsq_spline_fill(x, y_in, k=3, target_spacing=20_000.0)

    kept = np.ones(x.shape, bool)
    kept[gap] = False
    np.testing.assert_array_equal(y_out[kept], y[kept])
    assert np.isfinite(y_out[gap]).all()
    np.testing.assert_allclose(y_out[gap], y[gap], atol=0.5)


def test_lsq_spline_fill_returns_input_when_too_few_points():
    x = np.array([0.0, 1.0])
    y = np.array([1.0, 2.0], dtype=np.float32)

    y_out, y_fit = rpf.lsq_spline_fill(x, y, k=3)

    np.testing.assert_array_equal(y_out, y)
    assert np.isnan(y_fit).all()


def test_vertical_density_soft_clamp_pulls_in_layover_return():
    x = np.linspace(0, 17_000, 60)
    y = (10.0 + 0.0001 * x).astype(np.float32)
    y_in = y.copy()
    y_in[30] += 2.0  # vertically stacked (layover) return, same along-river bin

    out = rpf.vertical_density_soft_clamp(
        x,
        y_in,
        bin_width=18_000.0,
        y_bin=0.1,
        fixed_halfw=0.2,
        mad_factor=1.8,
        min_count=4,
        min_mode_count=2,
        slope_gain=0.35,
    )

    assert abs(out[30] - y[30]) < abs(y_in[30] - y[30])
