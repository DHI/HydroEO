"""Filter stages for longitudinal river WSE profile cleaning.

Each function operates on a 1-D along-river distance array ``x`` (metres)
and a co-located value array ``y`` (e.g. water-surface elevation), and is
independent of any Project/config state so it can be unit tested directly
with synthetic arrays. Used by :mod:`HydroEO.satellites.swot.river_profile`.
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import make_lsq_spline
from scipy.ndimage import binary_dilation


def robust_line(xw, yw, huber_k: float = 1.5, iters: int = 2):
    """Weighted least squares line with Huber IRLS re-weighting."""
    w = np.ones_like(yw, float)
    X = np.vstack([xw, np.ones_like(xw)]).T
    a, b = 0.0, 0.0
    for _ in range(max(0, iters)):
        WX = X * w[:, None]
        XtWX = WX.T @ X
        XtWy = WX.T @ yw
        try:
            a, b = np.linalg.solve(XtWX, XtWy)
        except np.linalg.LinAlgError:
            a, b = np.polyfit(xw, yw, 1)
        r = yw - (a * xw + b)
        med = np.nanmedian(r)
        mad = np.nanmedian(np.abs(r - med)) + 1e-6
        scale = 1.4826 * mad
        u = np.abs(r) / (huber_k * scale)
        w = np.where(u <= 1.0, 1.0, 1.0 / u)
    r = yw - (a * xw + b)
    return a, b, r


def sliding_window(xv, half_width: float):
    """
    ``xv`` must be sorted ascending. Yields ``(j, j0, j1)`` for every index
    ``j``, where ``xv[j0:j1+1]`` are exactly the points with
    ``xv[j]-half_width <= x <= xv[j]+half_width``.

    True O(n) two-pointer scan: both edges only ever move forward.
    """
    n = xv.size
    j0 = 0
    j1 = -1
    for j in range(n):
        xj = xv[j]
        while j0 < j and xv[j0] < xj - half_width:
            j0 += 1
        if j1 < j:
            j1 = j
        while j1 + 1 < n and xv[j1 + 1] <= xj + half_width:
            j1 += 1
        yield j, j0, j1


def vertical_density_soft_clamp(
    x,
    y_in,
    bin_width: float,
    y_bin: float,
    fixed_halfw: float,
    mad_factor: float,
    min_count: int,
    min_mode_count: int,
    slope_gain: float,
    huber_k: float = 1.5,
):
    """
    Detrend within along-river bins and softly clamp outliers to the
    dominant vertical mode of each bin (handles vertically-stacked layover
    returns without hard-masking them).
    """
    x = np.asarray(x, float)
    y = np.asarray(y_in, np.float32)
    out = y.copy()
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < min_count:
        return out
    xs_all = x[m]
    ys_all = y[m]
    order = np.argsort(xs_all)
    xs_all = xs_all[order]
    ys_all = ys_all[order]
    idx_all = np.where(m)[0][order]
    x0, x1 = xs_all.min(), xs_all.max()
    n_bins = max(1, int(np.ceil((x1 - x0) / float(bin_width))))
    for b in range(n_bins):
        lo = x0 + b * bin_width
        hi = lo + bin_width
        # last bin is closed on the right so a point exactly at x1 isn't
        # dropped when (x1 - x0) is an exact multiple of bin_width
        sel = (xs_all >= lo) & (xs_all < hi if b < n_bins - 1 else xs_all <= hi)
        if sel.sum() < min_count:
            continue
        xs = xs_all[sel]
        ys = ys_all[sel]
        ids = idx_all[sel]
        xc = xs.mean()
        xx = xs - xc

        a, b0, _ = robust_line(xx, ys, huber_k=huber_k, iters=2)
        y_detr = ys - (a * xx + b0)

        y_min, y_max = np.nanmin(y_detr), np.nanmax(y_detr)
        if not np.isfinite(y_min) or not np.isfinite(y_max) or y_max <= y_min:
            continue
        edges = np.arange(y_min, y_max + y_bin, y_bin, dtype=float)
        if edges.size < 3:
            continue
        counts, edges = np.histogram(y_detr, bins=edges)
        if counts.max() < min_mode_count:
            continue
        k_mode = np.argmax(counts)
        mode_center = 0.5 * (edges[k_mode] + edges[k_mode + 1])
        in_mode = (y_detr >= edges[k_mode]) & (y_detr <= edges[k_mode + 1])
        if in_mode.sum() >= 3:
            local_med = np.nanmedian(y_detr[in_mode])
            local_mad = np.nanmedian(np.abs(y_detr[in_mode] - local_med))
        else:
            local_med = np.nanmedian(y_detr)
            local_mad = np.nanmedian(np.abs(y_detr - local_med))
        half_adapt = mad_factor * (1.4826 * local_mad)
        half_slope = slope_gain * abs(a) * (bin_width * 0.5)
        halfw = max(fixed_halfw, half_adapt) + half_slope
        y_detr_clamped = np.clip(y_detr, mode_center - halfw, mode_center + halfw)
        out[ids] = (a * xx + b0 + y_detr_clamped).astype(np.float32)
    return out


def hampel_1d_meter(
    x,
    y_in,
    win_m: float,
    k: float,
    action: str,
    min_valid: int,
    replace_mode: str,
    huber_k: float = 1.5,
    huber_iters: int = 2,
):
    """
    Distance-windowed Hampel filter: flags points whose residual from a
    robust local trend exceeds ``k`` scaled MADs, then masks or replaces
    them per ``action``/``replace_mode``.
    """
    x = np.asarray(x, float)
    y = np.asarray(y_in, np.float32)
    out = y.copy()
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < min_valid:
        return out, np.zeros_like(y, bool)
    idx_all = np.arange(y.size)[valid]
    xv = x[valid]
    yv = y[valid]
    order = np.argsort(xv)
    idx_all = idx_all[order]
    xv = xv[order]
    yv = yv[order]
    changed = np.zeros_like(y, bool)

    for j, j0, j1 in sliding_window(xv, win_m):
        if j1 - j0 + 1 < min_valid:
            continue
        yy = yv[j0 : j1 + 1]
        xx = xv[j0 : j1 + 1] - xv[j]
        a, b, r = robust_line(xx, yy, huber_k=huber_k, iters=huber_iters)
        r_med = np.nanmedian(r)
        mad = np.nanmedian(np.abs(r - r_med)) + 1e-6
        thr = k * 1.4826 * mad
        res_j = yv[j] - b
        if np.abs(res_j) > thr:
            k_full = idx_all[j]
            changed[k_full] = True
            if action == "mask":
                out[k_full] = np.nan
            else:
                out[k_full] = np.float32(b if replace_mode == "trend" else np.nanmedian(yy))
    return out.astype(np.float32), changed


def rolling_quantile_trend(
    x,
    y_in,
    win_m: float,
    q: float,
    min_valid: int,
    robust_iters: int = 2,
    huber_k: float = 1.5,
):
    """Local robust-trend + rolling quantile of the residual, per window."""
    x = np.asarray(x, float)
    y = np.asarray(y_in, np.float32)
    out = np.full_like(y, np.nan, np.float32)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < min_valid:
        return out
    idx = np.arange(x.size)[m]
    xv = x[m]
    yv = y[m]
    order = np.argsort(xv)
    xv = xv[order]
    yv = yv[order]
    idx = idx[order]

    for j, j0, j1 in sliding_window(xv, win_m):
        ww = yv[j0 : j1 + 1]
        if ww.size < min_valid:
            continue
        xx = xv[j0 : j1 + 1] - xv[j]
        a, b, r = robust_line(xx, ww, huber_k=huber_k, iters=robust_iters)
        qres = np.nanpercentile(r, q * 100.0).astype(np.float32)
        out[idx[j]] = np.float32(b + qres)
    return out


def rolling_density(x, y, total_win_m: float):
    """Count of finite ``y`` values within +/-``total_win_m``/2 of each ``x[i]``."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = x.size
    order = np.argsort(x)
    xv = x[order]
    finite_v = np.isfinite(y[order]).astype(float)
    cum = np.concatenate([[0.0], np.cumsum(finite_v)])
    half = float(total_win_m) * 0.5

    counts_sorted = np.empty(n, float)
    for j, j0, j1 in sliding_window(xv, half):
        counts_sorted[j] = cum[j1 + 1] - cum[j0]

    counts = np.empty(n, float)
    counts[order] = counts_sorted
    return counts


def dilate_mask(mask, k: int = 0):
    mask = np.asarray(mask, bool)
    if k <= 0:
        return mask
    return binary_dilation(mask, structure=np.ones(2 * k + 1, dtype=bool))


def _dedup_xy(xs, ys):
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    order = np.argsort(xs)
    xs = xs[order]
    ys = ys[order]
    ux, uy = [], []
    i = 0
    n = xs.size
    while i < n:
        j = i + 1
        while j < n and xs[j] == xs[i]:
            j += 1
        block = ys[i:j]
        block = block[np.isfinite(block)]
        if block.size:
            ux.append(xs[i])
            uy.append(float(np.median(block)))
        i = j
    return np.asarray(ux, float), np.asarray(uy, float)


def _build_knot_vector(xu, k: int, inner_knots: int | None, tgt_spacing: float):
    if xu.size < (k + 2):
        return None
    x0, x1 = float(xu[0]), float(xu[-1])
    if inner_knots is None:
        n_inner = max(0, int(np.floor((x1 - x0) / float(tgt_spacing))) - 1)
        n_inner = min(n_inner, max(0, xu.size // 8))
    else:
        n_inner = max(0, int(inner_knots))
    if n_inner > 0:
        t_inner = np.quantile(xu, np.linspace(0, 1, n_inner + 2)[1:-1])
    else:
        t_inner = np.array([], float)
    t = np.r_[np.repeat(x0, k + 1), t_inner, np.repeat(x1, k + 1)]
    return t


def lsq_spline_fill(
    x,
    y_in,
    k: int,
    counts=None,
    inner_knots: int | None = None,
    target_spacing: float = 18_000.0,
    weight_scheme: str = "by_density",
):
    """Fit an LSQ spline to the valid points and paste it only into NaN gaps."""
    x = np.asarray(x, float)
    y = np.asarray(y_in, float).copy()
    m = np.isfinite(y)
    if np.count_nonzero(m) < (k + 2):
        return y.astype(np.float32), np.full_like(y, np.nan, np.float32)

    xu, yu = _dedup_xy(x[m], y[m])
    if xu.size < (k + 2):
        return y.astype(np.float32), np.full_like(y, np.nan, np.float32)

    t = _build_knot_vector(xu, k=k, inner_knots=inner_knots, tgt_spacing=target_spacing)
    if t is None:
        return y.astype(np.float32), np.full_like(y, np.nan, np.float32)

    if counts is not None and weight_scheme == "by_density":
        idx = np.searchsorted(x, xu)
        idx = np.clip(idx, 0, counts.size - 1)
        w = np.asarray(counts[idx], float)
        w = np.maximum(w, 1.0)
    else:
        w = None

    try:
        spl = make_lsq_spline(xu, yu, t, k=k, w=w)
    except Exception:
        spl = make_lsq_spline(xu, yu, t, k=k)

    y_fit = spl(x)
    inside = (x >= xu[0]) & (x <= xu[-1])
    holes = inside & (~np.isfinite(y))
    y[holes] = y_fit[holes]

    y_out = y.astype(np.float32)
    y_fit = y_fit.astype(np.float32)
    y_fit[~inside] = np.nan
    return y_out, y_fit
