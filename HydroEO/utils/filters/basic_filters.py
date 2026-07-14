"""simple filters that can be applied to sat timeseries objects"""

import logging

import numpy as np
import pandas as pd
from scipy.stats import theilslopes
from sklearn.svm import SVR

from datetime import datetime

logger = logging.getLogger(__name__)

_EARTH_RADIUS_KM = 6371.0088


def _haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance (km) between two points given in degrees."""
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, lat2, lon2))
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return 2 * _EARTH_RADIUS_KM * np.arcsin(np.sqrt(a))


def _along_track_distance_km(lats, lons):
    """
    Cumulative along-track distance (km) for points already ordered along
    the ground track (see _order_along_track for how that order is chosen).
    """
    lats = np.asarray(lats, dtype=float)
    lons = np.asarray(lons, dtype=float)
    if len(lats) < 2:
        return np.zeros_like(lats)
    step = _haversine_km(lats[:-1], lons[:-1], lats[1:], lons[1:])
    return np.concatenate([[0.0], np.cumsum(step)])


def _order_along_track(lats, lons):
    """
    Return the index order that best approximates along-track order for a
    single pass, given only lat/lon (no reliable sub-daily timestamp to sort
    by). Sorts by whichever of latitude/longitude spans the larger range
    within this group, since a single pass is normally close to monotonic in
    that coordinate (avoids failing near-equatorial, near-east-west passes
    where latitude barely changes).
    """
    lat_range = np.ptp(lats) if len(lats) else 0
    lon_range = np.ptp(lons) if len(lons) else 0
    sort_vals = lats if lat_range >= lon_range else lons
    return np.argsort(sort_vals)


def _resolve_pass_groups(df, date_key, pass_key=None, platform_key=None, orbit_key=None):
    """
    Resolve a grouping key identifying individual physical passes/crossings,
    generic across missions (this function knows nothing about which real
    column names any given mission uses -- that mapping happens when the
    Timeseries object is constructed).

    Resolved per ROW, not per dataframe -- this matters as soon as sources
    from different missions are concatenated together (see
    Timeseries.concat), since one mission's rows may have a valid pass_key
    while another's are NaN. Each row falls back independently through:
      1. `pass_key`, if non-null for that row.
      2. (date_key, platform_key, orbit_key), if platform_key and orbit_key
         are both non-null for that row.
      3. date_key alone.

    Returns a pandas Series of group labels (strings), same index as df.
    """
    n = len(df)
    group = pd.Series(pd.NA, index=df.index, dtype=object)

    if pass_key and pass_key in df.columns:
        has_pass = df[pass_key].notna()
        group.loc[has_pass] = "pass:" + df.loc[has_pass, pass_key].astype(str)

    remaining = group.isna()
    if (
        remaining.any()
        and platform_key
        and orbit_key
        and platform_key in df.columns
        and orbit_key in df.columns
    ):
        has_po = remaining & df[platform_key].notna() & df[orbit_key].notna()
        group.loc[has_po] = (
            "po:"
            + df.loc[has_po, date_key].astype(str)
            + "_"
            + df.loc[has_po, platform_key].astype(str)
            + "_"
            + df.loc[has_po, orbit_key].astype(str)
        )

    remaining = group.isna()
    if remaining.any():
        group.loc[remaining] = "date:" + df.loc[remaining, date_key].astype(str)

    return group


def elevation_filter(timeseries, height_range):
    min_height, max_height = height_range
    timeseries.df = timeseries.df.loc[timeseries.df[timeseries.height_key] > min_height]
    timeseries.df = timeseries.df.loc[timeseries.df[timeseries.height_key] < max_height]

    timeseries.df = timeseries.df.reset_index(drop=True)
    timeseries.df = timeseries.df.sort_values(by=timeseries.date_key)


def mad_filter(timeseries, threshold=5):
    # calculate support to make sure we can make a statistical decision, otheriwse leave
    if len(timeseries.df) >= 30:
        # calculate standard deviation and remove obvious outliers
        med = np.median(timeseries.df.height)
        abs_dev = np.abs(timeseries.df.height - med)
        mad = np.median(abs_dev)
        timeseries.df = timeseries.df.loc[abs_dev < threshold * mad]

        timeseries.df = timeseries.df.reset_index(drop=True)
        timeseries.df = timeseries.df.sort_values(by=timeseries.date_key)

    return timeseries


def daily_mad_error(
    timeseries, window_km=None, reg_weight=0.1, reg_default=0.5,
    min_window_points=4, max_theilsen_points=60,
):
    """
    Assign each observation an error estimate for use as the Kalman filter's
    observation uncertainty.

    Rows with a pre-supplied, mission-native formal uncertainty (see
    timeseries.preset_error_key, e.g. SWOT's wse_u) skip computed ADM
    entirely and use that value directly -- SWOT's LakeSP WSE is already
    one integrated value per crossing, not raw along-track points, so
    there is nothing to compute local dispersion or a local trend over.
    All other rows fall through to one of two modes:

    window_km=None (default): original behavior. Groups observations by
        calendar day and uses the absolute deviation from the daily median
        height as the error, regularized by day-level observation count.

    window_km=<value>: along-track, distance-based sliding window per
        physical pass (see _resolve_pass_groups for how a pass is
        identified). Within +/- window_km of each point, fits a robust
        (Theil-Sen) local linear trend against along-track distance and uses
        the deviation from that local trend as the error -- rather than the
        local median -- so a real, physical along-track slope (e.g. a river
        gradient) isn't mistaken for noise. Falls back to a local median
        within the window if fewer than min_window_points fall inside it,
        and further to reg_default alone if the pass has only 1 point.
        Requires timeseries.lat_key/lon_key columns; if they are missing or
        empty, this silently falls back to the window_km=None behavior with
        a warning, rather than failing the whole pipeline.

    Parameters
    ----------
    timeseries : Timeseries
    window_km : float, optional
        Along-track half-width (km) of the sliding window. None (default)
        keeps the original day-grouped-median behavior.
    reg_weight : float, optional
        Regularization numerator, divided by the local point count.
    reg_default : float, optional
        Regularization used when a group/window has only 1 point.
    min_window_points : int, optional
        Minimum points required inside a window to fit a local Theil-Sen
        trend rather than falling back to a local median. Default 4.
    max_theilsen_points : int, optional
        Cap on points fed to a single Theil-Sen fit. Theil-Sen is O(k^2) in
        window size (pairwise slopes), so dense along-track data (e.g.
        ATL13) can make an uncapped fit the dominant cost. Windows larger
        than this are evenly subsampled down to this many points rather
        than truncated, to keep coverage across the whole window. Default
        60 (60^2 = 3600 pairs per fit, vs. e.g. 800^2 = 640000 uncapped).
    """
    # sort inplace the timeseries object
    timeseries.df = timeseries.df.sort_values(by=timeseries.date_key).reset_index(
        drop=True
    )

    full_df = timeseries.df
    date_key = timeseries.date_key
    height_key = timeseries.height_key
    error_key = timeseries.error_key
    lat_key = getattr(timeseries, "lat_key", None)
    lon_key = getattr(timeseries, "lon_key", None)
    preset_error_key = getattr(timeseries, "preset_error_key", None)

    # Split off rows that already carry a formal, mission-supplied
    # uncertainty -- they bypass everything below.
    has_preset = pd.Series(False, index=full_df.index)
    if preset_error_key and preset_error_key in full_df.columns:
        has_preset = full_df[preset_error_key].notna()

    error_full = pd.Series(np.nan, index=full_df.index, dtype=float)
    if has_preset.any():
        error_full.loc[has_preset] = (
            full_df.loc[has_preset, preset_error_key].astype(float).clip(lower=1e-6)
        )

    df = full_df.loc[~has_preset].copy()

    if df.empty:
        timeseries.df[error_key] = error_full.values
        return timeseries

    have_coords = (
        window_km is not None
        and lat_key
        and lon_key
        and lat_key in df.columns
        and lon_key in df.columns
        and df[lat_key].notna().any()
        and df[lon_key].notna().any()
    )

    if window_km is not None and not have_coords:
        logger.warning(
            "daily_mad_error: window_km given but lat_key/lon_key ('%s'/'%s') "
            "not found or empty -- falling back to day-grouped-median ADM.",
            lat_key,
            lon_key,
        )

    if not have_coords:
        # ----- original day-grouped-median behavior -----
        day_key = df[date_key].dt.floor("D")

        medval = df.groupby(day_key).median(numeric_only=True)[height_key]
        day_grp = df.groupby(day_key).count()[height_key]

        reg = reg_weight / day_grp
        reg[day_grp == 1] = reg_default

        med_map = medval.to_dict()
        reg_map = reg.to_dict()
        df["med"] = day_key.map(med_map)
        df["reg"] = day_key.map(reg_map)

        error = np.abs(df[height_key] - df["med"]) + df["reg"]
        error = (
            error.replace([np.inf, -np.inf], np.nan)
            .fillna(reg_default)
            .clip(lower=1e-6)
        )

        error_full.loc[df.index] = error.values
        timeseries.df[error_key] = error_full.values
        return timeseries

    # ----- windowed, slope-aware path -----
    groups = _resolve_pass_groups(
        df,
        date_key,
        getattr(timeseries, "pass_key", None),
        getattr(timeseries, "platform_key", None),
        getattr(timeseries, "orbit_key", None),
    )

    error = pd.Series(np.nan, index=df.index, dtype=float)

    for _, idx in df.groupby(groups).groups.items():
        idx = np.asarray(idx)
        sub = df.loc[idx]
        n = len(idx)

        if n == 1:
            error.loc[idx] = reg_default
            continue

        order = _order_along_track(sub[lat_key].values, sub[lon_key].values)
        idx_sorted = idx[order]
        lat = sub[lat_key].values[order]
        lon = sub[lon_key].values[order]
        heights = sub[height_key].values[order]
        dist = _along_track_distance_km(lat, lon)

        for i in range(n):
            window_mask = np.abs(dist - dist[i]) <= window_km
            n_window = int(window_mask.sum())

            if n_window >= min_window_points:
                win_h = heights[window_mask]
                win_d = dist[window_mask]
                if n_window > max_theilsen_points:
                    # Theil-Sen is O(k^2) in window size (pairwise slopes);
                    # dense along-track data (e.g. ATL13) can otherwise make
                    # this the dominant cost. Evenly subsample rather than
                    # truncate, to keep coverage across the whole window.
                    sel = np.linspace(
                        0, n_window - 1, max_theilsen_points
                    ).round().astype(int)
                    win_h = win_h[sel]
                    win_d = win_d[sel]
                slope, intercept, _, _ = theilslopes(win_h, win_d)
                local_ref = intercept + slope * dist[i]
            else:
                local_ref = np.median(heights[window_mask])

            reg = reg_weight / n_window if n_window > 1 else reg_default
            error.loc[idx_sorted[i]] = abs(heights[i] - local_ref) + reg

    error = error.replace([np.inf, -np.inf], np.nan).fillna(reg_default).clip(lower=1e-6)
    error_full.loc[df.index] = error.values
    timeseries.df[error_key] = error_full.values

    return timeseries


def daily_mean_merge(timeseries):
    # aggregates values that occur on the same day to a mean value

    dct = {
        "number": "mean",
        "object": lambda col: col.mode() if col.nunique() == 1 else np.nan,
    }

    groupby_cols = [timeseries.date_key]
    dct = {
        k: v
        for i in [
            {
                col: agg
                for col in timeseries.df.select_dtypes(tp).columns.difference(
                    groupby_cols
                )
            }
            for tp, agg in dct.items()
        ]
        for k, v in i.items()
    }
    timeseries.df = timeseries.df.groupby(groupby_cols).agg(
        **{k: (k, v) for k, v in dct.items()}
    )

    timeseries.df[timeseries.date_key] = pd.to_datetime(timeseries.df.index)
    timeseries.df = timeseries.df.reset_index(drop=True)
    timeseries.df = timeseries.df.sort_values(by=timeseries.date_key)

    return timeseries


def hampel(timeseries, k=7, t0=3):
    """
    vals: pandas series of values from which to remove outliers
    k: size of window (including the sample; 7 is equal to 3 on either side of value)
    """
    vals = timeseries.df.loc[:, timeseries.height_key]

    # Hampel Filter
    L = 1.4826
    rolling_median = vals.rolling(k).median()
    difference = np.abs(rolling_median - vals)
    median_abs_deviation = difference.rolling(k).median()
    threshold = t0 * L * median_abs_deviation
    outlier_idx = difference > threshold
    vals[outlier_idx] = np.nan

    timeseries.df[timeseries.height_key] = vals.values

    return timeseries


def rolling_median(timeseries, window=7):
    timeseries.df[timeseries.height_key] = (
        timeseries.df[timeseries.height_key].rolling(window).median()
    )

    return timeseries


def _run_svr_linear(heights, err=0.1, epsilon=0.1, max_iter=5000):
    """
    Linear Support Vector Regression outlier filter.

    Fits a free-slope linear SVR through heights along-track and flags
    points that deviate from that fitted line by more than `err` as
    outliers. The fitted slope is used only to decide which points are
    outliers here -- it is not removed from the retained heights (that is
    a separate, later step).

    Parameters
    ----------
    heights : array
        Heights to be fit, already ordered along-track.
    err : Float, optional
        Allowed deviation from the linear fit. The default is 0.1.
    epsilon : Float, optional
        "Epsilon in the epsilon-SVR model.
        It specifies the epsilon-tube within which no penalty is associated in the training loss
        function with points predicted within a distance epsilon from the actual value."
        (from https://scikit-learn.org/stable/modules/generated/sklearn.svm.SVR.html)
        The default is .1.
    max_iter : int, optional
        Hard cap on the underlying solver's iterations. sklearn's SVR
        defaults to max_iter=-1 (no cap at all) -- on a pathological or
        unexpectedly large/ill-conditioned group, that can make a single
        fit call run for a very long time with no visible progress,
        indistinguishable from a genuine hang. Default 5000; if this is
        hit, a warning is logged and the (not fully converged, but usually
        still reasonable) fit is used rather than blocking indefinitely.

    Returns
    -------
    array
        Integer positions (0-based, within `heights`) of the points that
        survive the filter.
    """
    heights = np.asarray(heights, dtype=float)

    # sequential x axis along track
    x = np.arange(0, len(heights))
    x = np.vstack(
        [x, np.ones(len(x))]
    ).T  # Extend x data to contain another row vector of 1s
    y = heights

    # make SVR kernel and fit with confidence bounds
    svr_rbf = SVR(kernel="linear", epsilon=epsilon, max_iter=max_iter)
    rbf = svr_rbf.fit(x, y)
    if getattr(rbf, "n_iter_", 0) is not None and np.any(
        np.asarray(rbf.n_iter_) >= max_iter
    ):
        logger.warning(
            "_run_svr_linear: SVR hit max_iter=%d without full convergence "
            "on a group of %d points -- result used anyway, but check "
            "whether this group is unexpectedly large (pass-grouping "
            "issue) or has pathological/duplicate values.",
            max_iter, len(heights),
        )
    uconf = rbf.predict(x) + err
    lconf = rbf.predict(x) - err

    # return the positions of the retained points (not the values), so the
    # caller can keep every column for those rows, not just height
    return np.where((y >= lconf) & (y <= uconf))[0]


def svr_linear(timeseries, err=0.1, epsilon=0.1, max_iter=5000, warn_group_size=500):
    """
    Along-track linear-SVR outlier filter, grouped by physical pass (see
    _resolve_pass_groups) rather than calendar day, and preserving every
    column of the surviving rows (not just date/height).

    Parameters
    ----------
    timeseries : Timeseries
    err : float, optional
        Allowed deviation from the local linear fit (m). Default 0.1.
    epsilon : float, optional
        SVR epsilon-tube. Default 0.1.
    max_iter : int, optional
        Passed through to the underlying SVR fit; see _run_svr_linear.
    warn_group_size : int, optional
        Log a warning up front for any single pass/group larger than this,
        before attempting to fit it -- SVR's cost scales badly with group
        size, so an unexpectedly large group (e.g. a pass-grouping fallback
        issue lumping many points together) is worth surfacing immediately
        rather than only being discovered via a slow or capped fit.
    """
    df = timeseries.df.copy()
    date_key = timeseries.date_key
    height_key = timeseries.height_key

    groups = _resolve_pass_groups(
        df,
        date_key,
        getattr(timeseries, "pass_key", None),
        getattr(timeseries, "platform_key", None),
        getattr(timeseries, "orbit_key", None),
    )

    group_items = df.groupby(groups).groups.items()
    group_sizes = {g: len(idx) for g, idx in group_items}
    oversized = {g: n for g, n in group_sizes.items() if n > warn_group_size}
    if oversized:
        logger.warning(
            "svr_linear: %d group(s) exceed warn_group_size=%d (largest: "
            "%s with %d points) -- SVR fit cost scales badly with group "
            "size; if this is unexpected, check pass_key/platform_key/"
            "orbit_key resolution for this timeseries.",
            len(oversized), warn_group_size,
            max(oversized, key=oversized.get), max(oversized.values()),
        )

    keep_idx = []
    for _, idx in group_items:
        idx = np.asarray(idx)
        if len(idx) < 2:
            # nothing to compare a single point against; keep it as-is
            keep_idx.extend(idx)
            continue
        heights = df.loc[idx, height_key].values
        local_keep = _run_svr_linear(heights, err=err, epsilon=epsilon, max_iter=max_iter)
        keep_idx.extend(idx[local_keep])

    timeseries.df = df.loc[np.sort(np.asarray(keep_idx))].reset_index(drop=True)

    return timeseries


"""Kalman filter for estimating state of reservoir from noisy timeseries"""


def _update(obs, xk, cov_xx, height="height", error="ADM", n=1):
    """
    Update function for Kalman filter

    Parameters
    ----------
    obs : array
        DESCRIPTION.
    xk : float
        Previous prediction, x.
    cov_xx : float
        Uncertainty of x.
    height : string, optional
        DESCRIPTION. The default is 'height'.
    error : string, optional
        DESCRIPTION. The default is 'ADM'.
    n : Int, optional
        Grid size - only relevant for large lakes e.g. Not yet implemented. The default is 1.

    Returns
    -------
    xk_plus : Float
        Kalman filter updated value.
    cov_xx_plus : Float
        Covariance of the updated value.

    """
    obs = obs.copy()
    obs = obs.replace([np.inf, -np.inf], np.nan)
    obs = obs.dropna(subset=[height, error])

    # No valid observations for this epoch: keep predicted state unchanged.
    if obs.empty:
        return xk, cov_xx

    obs_list = obs[height].values
    m = len(obs_list)

    if n == 1:
        # Closed-form scalar update: combining m independent observations
        # (variances slk) with a Gaussian prior (xk, cov_xx) has an exact
        # closed form (precision-weighted average) -- mathematically
        # identical to the general matrix path below for n=1 (verified via
        # fuzz testing against it), but O(m) instead of O(m^3) from
        # np.linalg.pinv on an m x m matrix. n=1 is the only case used
        # anywhere in this codebase today (n>1 grid support is not yet
        # implemented -- see docstring), so this is the common path.
        slk = np.maximum((obs[error].values) ** 2, 1e-12)
        xk_s = float(np.asarray(xk).reshape(()))
        cov_s = float(np.asarray(cov_xx).reshape(()))
        prec = 1.0 / cov_s + np.sum(1.0 / slk)
        xk_plus = (xk_s / cov_s + np.sum(obs_list / slk)) / prec
        cov_xx_plus = 1.0 / prec
        return np.array(xk).reshape(1, 1) * 0 + xk_plus, np.array(cov_xx).reshape(1, 1) * 0 + cov_xx_plus

    lk = np.array(obs_list).reshape(m, n)
    Ak = np.ones((m, n))

    # compute Kalman matrix (weights of the innovation)
    slk = np.maximum((obs[error].values) ** 2, 1e-12)
    cov_lk = np.diag(slk)
    inv = Ak * cov_xx * np.transpose(Ak) + cov_lk
    if not np.all(np.isfinite(inv)):
        return xk, cov_xx
    inv = inv + np.eye(m) * 1e-8
    Kk = cov_xx * np.dot(np.transpose(Ak), np.linalg.pinv(inv))

    # update x
    xk_plus = xk + np.dot(Kk, lk - Ak * xk)

    # update sigma
    cov_xx_plus = cov_xx * (1 - np.dot(Kk, Ak))

    return xk_plus, cov_xx_plus


def _pred(xk, cov_xx, n=1, system_noise=0, system_noise_unc=0.0005):
    """
    Prediction function for Kalman filter

    Parameters
    ----------
    xk : float
        Updated prediction, x.
    cov_xx : float
        Uncertainty of x.
    n : Int, optional
        Grid size - only relevant for large lakes e.g. Not yet implemented. The default is 1.
    system_noise : float, optional
        System noise. The default is 0.
    system_noise_unc : float, optional
        Uncertainty (variance, m^2) of the system noise. The default is
        0.0005 (5 cm^2), matching DAHITI (Schwatke et al., 2015).

    Returns
    -------
    xk_plus : Float
        Kalman filter predicted value.
    cov_xx_plus : Float
        Covariance of the predicted value.

    """
    # Setup dynamic model with no deterministic contribution
    thetak = np.identity(n)
    hatk = np.identity(n)
    Qk = np.array(system_noise_unc).reshape(n, n)
    qk = np.array(system_noise).reshape(n, 1)
    noise = hatk.dot(Qk).dot(hatk.T)

    # pred x
    xk_next = thetak * xk + hatk * qk
    # predict covariance matrix
    cov_xx_next = thetak.dot(cov_xx).dot(thetak.T) + noise

    return xk_next, cov_xx_next


def fit_spatial_correction_model(
    df, lat_key, lon_key, height_key, date_key, ref_lat, ref_lon,
    min_points_per_day=20, min_days=3,
):
    """
    Fit a 1D spatial deviation model along a reservoir's most stable
    spatial axis, from a dense, wide-spanning source (e.g. ICESat-2).
    Unlike apply_distance_penalty (which only inflates uncertainty), this
    is a genuine height CORRECTION: it estimates how much a crossing's
    value tends to differ from the reservoir centroid's value as a
    function of position, so that correction can be applied to ANY
    mission's points -- letting a sparse source (e.g. Sentinel-3/6, with
    only 1-2 points per crossing and no ability to fit its own local
    trend) benefit from a dense source's spatial mapping of the reservoir.

    Method: for each day with at least min_points_per_day points, remove
    that day's own median (strips out real level change over time,
    leaving just the spatial deviation pattern for that day), then fit a
    separate 1D linear slope of residual vs. position along each of two
    orthogonal local axes (east-west, north-south). The axis whose slope
    is most CONSISTENT across independent days (highest weighted-mean /
    weighted-std ratio -- like a t-statistic) is taken as the real,
    physically meaningful axis; the other is treated as noise and
    discarded. This is deliberately more conservative than trusting a
    single pooled regression across all days pooled together, which can
    show a misleadingly non-zero coefficient on the noisy axis (confirmed
    empirically: a pooled fit suggested a coefficient on the noisy axis
    of similar magnitude to the real one, while per-day fits showed it
    flipping sign inconsistently -- pooling does not reliably distinguish
    a real, stable effect from incidental within-day sampling structure).

    IMPORTANT for update stability: this function has no notion of when
    it's called or what's been fit before -- it fits fresh from whatever
    df is passed in every time. If you re-fit every run as new data
    arrives, past corrections WILL shift retroactively. Persist the
    returned model (e.g. to JSON) and reuse it across runs rather than
    refitting automatically; only refit on an explicit recalibration
    trigger. See flows.py for the persistence wrapper.

    Parameters
    ----------
    df : DataFrame
        Points from the dense source only (e.g. filter to platform=="icesat2"
        before calling), with valid lat_key/lon_key/height_key/date_key.
    ref_lat, ref_lon : float
        Reference location (e.g. the reservoir polygon's own centroid --
        see flows._reservoir_centroid) that corrected heights are
        expressed relative to.
    min_points_per_day : int, optional
        Minimum points on a given day to attempt a local fit for that day.
        Default 20.
    min_days : int, optional
        Minimum number of qualifying days required to fit a model at all.
        Default 3 -- a single day's pattern could be a real but transient
        anomaly (e.g. a wind event), not a persistent geometric feature;
        requiring several independent days is what lets us distinguish
        "persistent, worth correcting for" from "one-off."

    Returns
    -------
    dict or None
        {"axis": "x" or "y", "slope_m_per_km": float, "ref_lat": float,
         "ref_lon": float, "n_days_used": int, "diagnostics": {...}}
        None if there isn't enough dense data to fit anything -- callers
        should treat this as "no correction available," not an error.
    """
    df = df.copy()
    day_key = df[date_key].dt.floor("D")
    day_median = df.groupby(day_key)[height_key].transform("median")
    residual = df[height_key] - day_median

    x_km = (df[lon_key] - ref_lon) * 111.0 * np.cos(np.radians(ref_lat))
    y_km = (df[lat_key] - ref_lat) * 111.0

    records = []
    for day, idx in df.groupby(day_key).groups.items():
        idx = np.asarray(idx)
        if len(idx) < min_points_per_day:
            continue
        r = residual.loc[idx].values
        for axis_name, coord in (("x", x_km.loc[idx].values), ("y", y_km.loc[idx].values)):
            if np.std(coord) < 1e-6:
                continue
            slope, intercept = np.polyfit(coord, r, 1)
            pred = slope * coord + intercept
            ss_res = np.sum((r - pred) ** 2)
            ss_tot = np.sum((r - r.mean()) ** 2)
            r2 = float(1 - ss_res / ss_tot) if ss_tot > 0 else float("nan")
            records.append(dict(day=day, axis=axis_name, n=len(idx), slope=slope, r2=r2))

    if not records:
        logger.warning(
            "fit_spatial_correction_model: no day had >= %d points -- "
            "not enough dense data to fit a model.", min_points_per_day,
        )
        return None

    fits = pd.DataFrame(records)
    n_days = fits["day"].nunique()
    if n_days < min_days:
        logger.warning(
            "fit_spatial_correction_model: only %d qualifying day(s) "
            "available (need >= %d) -- not fitting a model.", n_days, min_days,
        )
        return None

    best_axis, best_slope, best_score = None, None, -np.inf
    diagnostics = {}
    for axis_name in ("x", "y"):
        sub = fits[fits.axis == axis_name]
        if sub.empty:
            continue
        w = sub.n.values
        mean_slope = float(np.average(sub.slope, weights=w))
        weighted_std = float(np.sqrt(np.average((sub.slope - mean_slope) ** 2, weights=w)))
        score = abs(mean_slope) / (weighted_std + 1e-9)
        diagnostics[axis_name] = dict(
            mean_slope_m_per_km=mean_slope, weighted_std_across_days=weighted_std,
            stability_score=float(score), n_days=int(sub["day"].nunique()),
        )
        if score > best_score:
            best_axis, best_slope, best_score = axis_name, mean_slope, score

    logger.info(
        "fit_spatial_correction_model: chose axis '%s' (score=%.2f) over "
        "the other (see diagnostics for both): %s", best_axis, best_score, diagnostics,
    )

    return {
        "axis": best_axis,
        "slope_m_per_km": best_slope,
        "ref_lat": float(ref_lat),
        "ref_lon": float(ref_lon),
        "n_days_used": int(n_days),
        "diagnostics": diagnostics,
    }


def apply_spatial_correction(timeseries, model):
    """
    Apply a fitted spatial correction model (see fit_spatial_correction_model)
    to every point (any mission) with valid lat_key/lon_key, referencing
    every crossing to the model's reference location along its dominant
    axis. Unlike apply_distance_penalty, this adjusts the height itself,
    not just its error/uncertainty.

    A no-op (with a log note, not an error) if model is None (nothing was
    fit -- e.g. not enough dense-source data yet) or if lat_key/lon_key
    aren't available on this timeseries (e.g. SWOT).
    """
    if model is None:
        logger.info("apply_spatial_correction: no model given -- skipping (no-op).")
        return timeseries

    lat_key = getattr(timeseries, "lat_key", None)
    lon_key = getattr(timeseries, "lon_key", None)
    have_coords = (
        lat_key
        and lon_key
        and lat_key in timeseries.df.columns
        and lon_key in timeseries.df.columns
        and timeseries.df[lat_key].notna().any()
    )
    if not have_coords:
        logger.info(
            "apply_spatial_correction: lat_key/lon_key not available -- "
            "skipping (no-op). Expected e.g. for SWOT."
        )
        return timeseries

    ref_lat, ref_lon = model["ref_lat"], model["ref_lon"]
    axis, slope = model["axis"], model["slope_m_per_km"]

    mask = timeseries.df[lat_key].notna() & timeseries.df[lon_key].notna()
    correction = pd.Series(0.0, index=timeseries.df.index)

    if axis == "x":
        coord = (timeseries.df.loc[mask, lon_key] - ref_lon) * 111.0 * np.cos(np.radians(ref_lat))
    else:
        coord = (timeseries.df.loc[mask, lat_key] - ref_lat) * 111.0

    correction.loc[mask] = slope * coord
    timeseries.df[timeseries.height_key] = timeseries.df[timeseries.height_key] - correction

    return timeseries


def apply_distance_penalty(timeseries, ref_lat, ref_lon, scale_per_km=0.05):
    """
    Inflate each observation's error (timeseries.error_key) by an amount
    proportional to its distance from a reference location -- a smooth,
    additive down-weighting for the Kalman filter.

    Motivation: local ADM only measures how much a crossing's points scatter
    around each other -- it says nothing about whether that crossing's
    average value is actually representative of the main reservoir body. A
    well-sampled crossing far upstream (subject to real slope effects, e.g.
    a riverine arm during drawdown) can have very LOW internal scatter (all
    its points agreeing closely with each other) while still being
    systematically biased relative to the main body -- ADM alone would
    treat that as a highly trustworthy observation, the opposite of what's
    wanted. This function adds a second, independent uncertainty term based
    purely on location, so Kalman's existing precision-weighting naturally
    discounts far-from-reference crossings without needing to guess at
    "representativeness" from point count or any other proxy.

    Call this AFTER daily_mad_error (requires timeseries.error_key to
    already be populated) and BEFORE kalman().

    Parameters
    ----------
    timeseries : Timeseries
    ref_lat, ref_lon : float
        Reference location, e.g. the reservoir polygon's own centroid (not
        a centroid of wherever satellite data happened to sample -- the
        true reservoir geometry is a more principled, sampling-independent
        reference), or a dam/outlet location if available.
    scale_per_km : float, optional
        Additional error (m) added per km of distance from the reference.
        Default 0.05 m/km -- e.g. a crossing 20 km upstream gets +1.0 m of
        inflated error. This is a real physical scale that should be tuned
        against reservoirs where the true magnitude of upstream slope bias
        is known, not trusted blindly at this default.

    Requires timeseries.lat_key/lon_key. If unavailable (e.g. SWOT, which
    has no per-observation coordinates -- see preset_error_key), this is a
    no-op with a log note, not an error: there's nothing to compute a
    distance from, and SWOT's own wse_u already supersedes ADM entirely for
    those rows regardless.
    """
    lat_key = getattr(timeseries, "lat_key", None)
    lon_key = getattr(timeseries, "lon_key", None)
    error_key = timeseries.error_key

    have_coords = (
        lat_key
        and lon_key
        and lat_key in timeseries.df.columns
        and lon_key in timeseries.df.columns
        and timeseries.df[lat_key].notna().any()
    )
    if not have_coords:
        logger.info(
            "apply_distance_penalty: lat_key/lon_key not available -- "
            "skipping (no-op). Expected e.g. for SWOT, which has no "
            "per-observation coordinates."
        )
        return timeseries

    if error_key not in timeseries.df.columns:
        logger.warning(
            "apply_distance_penalty: error_key '%s' not found -- call this "
            "after daily_mad_error, not before. Skipping (no-op).",
            error_key,
        )
        return timeseries

    mask = timeseries.df[lat_key].notna() & timeseries.df[lon_key].notna()
    dist_km = pd.Series(0.0, index=timeseries.df.index)
    if mask.any():
        dist_km.loc[mask] = _haversine_km(
            timeseries.df.loc[mask, lat_key].values,
            timeseries.df.loc[mask, lon_key].values,
            ref_lat,
            ref_lon,
        )

    timeseries.df[error_key] = timeseries.df[error_key] + scale_per_km * dist_km
    return timeseries


def kalman(timeseries, n=1, system_noise_unc=0.0005):
    """
    Run Kalman filter

    Parameters
    ----------
    timeseries : Timeseries
        Outlier-filtered timeseries to be used as input for the Kalman filter.
        Reads height_key, error_key, and date_key from the timeseries object
        itself (must contain those columns already, e.g. after daily_mad_error).
    n : Int, optional
        Grid size - only relevant for large lakes e.g. Not yet implemented. The default is 1.
    system_noise_unc : float, optional
        Variance (m^2) of the system noise injected at each prediction step,
        i.e. how much the state is allowed to drift between updates.
        The default is 0.0005 (5 cm^2), matching DAHITI (Schwatke et al., 2015).
        Increasing the value allows the filter to track
        raw observations much more closely (less smoothing) than intended.

    Returns
    -------
    xk_plus : array
        Kalman filter updated values of WSE.
    cov_xx_plus : array
        Covariance of the updated value.

    """

    df = timeseries.df.copy()
    date_key = timeseries.date_key
    height_key = timeseries.height_key
    error_key = timeseries.error_key

    # Group by calendar day, not exact timestamp. 
    day_key = df[date_key].dt.floor("D")
    dates = sorted(day_key.unique())
    xks = np.ones((n, len(dates))) * np.nan
    cov_xxs = np.ones((n, n, len(dates))) * np.nan

    # Group once up front
    grouped = {k: v for k, v in df.groupby(day_key)}

    # Initialize prediction and uncertainty matrices
    obs = grouped[dates[0]]
    lk = obs[height_key].values
    slk = (obs[error_key].values) ** 2

    # First prediction = value with smallest ADM and identity matrix of size n, n
    xks[:, 0] = lk[np.argmin(slk)]
    cov_xxs[:, :, 0] = np.identity(n)

    # Observation model
    for i, d in enumerate(dates):
        # Update
        obs = grouped[d]
        xk_plus, cov_xx_plus = _update(
            obs, xks[:, i], cov_xxs[:, :, i], height=height_key, error=error_key, n=n
        )

        # Predict
        xk1, cov_xx1 = _pred(
            xk_plus, cov_xx_plus, n=n, system_noise=0,
            system_noise_unc=system_noise_unc
        )
        xks[:, i] = xk_plus
        cov_xxs[:, :, i] = cov_xx_plus
        try:
            xks[:, i + 1] = xk1
            cov_xxs[:, :, i + 1] = cov_xx1
        except IndexError:
            pass

    # return a new dataframe with the filtered data
    df_kalman = pd.DataFrame(
        {date_key: dates, height_key: xks[0], "cov": cov_xxs[0][0]}
    )
    return df_kalman


def _run_svr_rbf(dates, heights, err=1.0, rbf_c=10000, gamma=0.0000438, epsilon=0.1,
                  max_iter=-1, max_fit_points=None):
    """
    Radial Base Function outlier filtering post-Kalman filter.
    This is run at virtual station level

    Parameters
    ----------
    dates : array
        Sorted dates of altimetry observations.
    heights : array
        Predicted (and updated) water surface elevation.
    err : Float, optional
        Observation uncertainty. The default is 1 m. As used in DAHITI for rivers
        (0.1 m can be used for lakes)
    C : TYPE, optional
        Regularization parameter in sklearn.svm.RBF. The default here is 10000 (as in DAHITI)
    gamma : TYPE, optional
        Kernel coefficient for ‘rbf’ in sklearn.svm.RBF
        Alternative values include float, 'scale' or 'auto'
        The default here is 0.0000438 (as in DAHITI)
    epsilon : TYPE, optional
        "Epsilon in the epsilon-SVR model.
        It specifies the epsilon-tube within which no penalty is associated in the training loss
        function with points predicted within a distance epsilon from the actual value."
        (from https://scikit-learn.org/stable/modules/generated/sklearn.svm.SVR.html)
        The default is .1.
    max_iter : int, optional
        Passed to the underlying SVR fit. UNLIKE _run_svr_linear, this
        defaults to -1 (unbounded) -- confirmed empirically that capping
        this specific configuration (RBF kernel + rbf_c=10000, very little
        regularization) produces non-monotonic, unreliable output quality:
        e.g. max_iter=5000 kept 10/5760 points, 50000 kept 827/5760, and
        100000 kept only 30/5760 -- only fully unbounded (23s on that real
        test case) gave the correct 2154/5760. Do not cap this without
        checking output row counts, not just wall-clock time. Use
        max_fit_points instead for a speedup that doesn't have this risk.
    max_fit_points : int, optional
        If set and there are more than this many points, the SVR is fit on
        an evenly-subsampled subset of this size (not truncated -- spread
        across the whole series), then applied via .predict() to ALL
        original points for the outlier decision. This reduces the cost
        driver (problem size going into the O(n^2)-ish solver) directly,
        rather than truncating iterations on the full problem -- much
        safer than capping max_iter for this kernel/C combination, since
        the underlying trend being fit is smooth and a representative
        subsample captures it about as well as the full set. None (default)
        disables this and fits on every point, as before.

    Returns
    -------
    rbf_filter : array
        Filter for Kalman filter predictions of WSE
    """
    # Day since start
    x = np.array(
        [0]
        + list(
            np.cumsum(
                pd.to_datetime(dates)[1:] - pd.to_datetime(dates)[:-1]
            ).days.astype("float32")
        )
    )
    # Kalman heights
    y = heights
    x = np.vstack([x, np.ones(len(x))]).T

    if max_fit_points and len(x) > max_fit_points:
        fit_sel = np.linspace(0, len(x) - 1, max_fit_points).round().astype(int)
        fit_x, fit_y = x[fit_sel], y[fit_sel]
    else:
        fit_x, fit_y = x, y

    svr_rbf = SVR(kernel="rbf", C=rbf_c, gamma=gamma, epsilon=epsilon, max_iter=max_iter)
    rbf = svr_rbf.fit(fit_x, fit_y)
    if getattr(rbf, "n_iter_", 0) is not None and np.any(
        np.asarray(rbf.n_iter_) >= max_iter > 0
    ):
        logger.warning(
            "_run_svr_rbf: SVR hit max_iter=%d without full convergence on "
            "%d fit points -- result used anyway.", max_iter, len(fit_x),
        )
    corr = rbf.predict(x)

    uconf = corr + err
    lconf = corr - err

    rbf_filter = np.where((y >= lconf) & (y <= uconf))

    return rbf_filter


def _year_fraction(dt):
    start = datetime(dt.year, 1, 1).toordinal()
    year_length = datetime(dt.year + 1, 1, 1).toordinal() - start
    return dt.year + float(dt.toordinal() - start) / year_length


def svr_radial(timeseries, err=1.0, rbf_c=10000, gamma=0.0000438, epsilon=0.1,
                max_iter=-1, max_fit_points=None):
    """
    Radial-basis post-Kalman outlier filter, run at virtual station level.

    Parameters
    ----------
    timeseries : Timeseries
        Kalman-filtered timeseries to be filtered.
    err : float, optional
        Observation uncertainty / half-width of the confidence band (m).
        DAHITI uses about 1 m for rivers and 0.1 m for lakes -- since
        flow.py runs separate river vs. reservoir flows, pass the
        appropriate value explicitly from there rather than relying on
        this default. The default here (1.0) is the more conservative
        (river-like) choice.
    rbf_c : float, optional
        Regularization parameter C of the RBF SVR. The default is 10000,
        matching DAHITI (lower values make the fit more regularized 
        / less able to follow the data than intended).
    gamma : float, optional
        RBF kernel coefficient. The default is 0.0000438, as in DAHITI.
    epsilon : float, optional
        Epsilon-tube of the SVR (see sklearn.svm.SVR). The default is 0.1.
    max_iter : int, optional
        Passed through to the underlying SVR fit; see _run_svr_rbf for why
        this defaults to -1 (unbounded) rather than being capped.
    max_fit_points : int, optional
        Passed through to _run_svr_rbf -- the safe speedup for this stage
        (fit on an evenly-subsampled subset, apply to all points), rather
        than capping max_iter.

    Returns
    -------
    vs : DataFrame
        Filtered (date, height) virtual station timeseries.
    """
    df = timeseries.df.copy()
    date_key = timeseries.date_key
    height_key = timeseries.height_key

    rbf_filter = _run_svr_rbf(
        df[date_key].values,
        df[height_key].values,
        err=err,
        rbf_c=rbf_c,
        gamma=gamma,
        epsilon=epsilon,
        max_iter=max_iter,
        max_fit_points=max_fit_points,
    )

    nb_obs = df.groupby(date_key).count().reset_index()

    nb_obs["nb_obs"] = nb_obs[date_key]

    # Create new dataframe:
    vs = pd.DataFrame(
        {
            date_key: df[date_key].values[rbf_filter],
            height_key: df[height_key].values[rbf_filter],
        }
    )

    return vs