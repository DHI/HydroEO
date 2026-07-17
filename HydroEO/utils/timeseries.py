from dataclasses import dataclass
import logging
import warnings
import os
import pandas as pd

import HydroEO.utils.filters.basic_filters as fltrs
from HydroEO.utils import general

logger = logging.getLogger(__name__)


def _parse_bin_days(time_bin):
    """Parse a simple '<N>D' style string into an integer number of days."""
    s = str(time_bin).strip().upper()
    if s.endswith("D"):
        return int(s[:-1]) if s[:-1] else 1
    raise ValueError("time_bin must be a string like '1D' or '5D'")


@dataclass
class Timeseries:
    df: pd.DataFrame
    date_key: str = "date"
    height_key: str = "height"
    error_key: str = "ADM"
    lat_key: str = "lat"
    lon_key: str = "lon"
    pass_key: str = "pass"
    platform_key: str = "platform"
    orbit_key: str = "orbit"
    preset_error_key: str = None

    def __post_init__(self):
        if self.height_key not in self.df.columns:
            logger.warning(
                "Height column '%s' missing from dataframe", self.height_key
            )
            # error_key is intentionally not checked here: it's populated
            # later in the pipeline (by daily_mad_error), not expected to
            # exist at construction time. Functions that actually need it
            # (daily_mad_error, kalman) will raise a clear error at the
            # point they read it if it's genuinely still missing then.

        if self.date_key not in self.df.columns:
            if isinstance(self.df.index, pd.DatetimeIndex):
                self.df[self.date_key] = self.df.index
            else:
                logger.warning(
                    "date_key '%s' is not in dataframe and index is not a "
                    "DatetimeIndex, so it could not be auto-populated",
                    self.date_key,
                )

    def clean(self, filters: list, filter_params: dict = None):
        filter_params = filter_params or {}

        ##### Ensure that provided filters are supported
        supported_filters = [
            "elevation",
            "MAD",
            "daily_mean",
            "hampel",
            "rolling_median",
        ]

        for filter in filters:
            if filter not in supported_filters:
                warnings.warn(
                    f"{filter} is not in list of supported filters and will not be applied"
                )
        filters = [filter for filter in filters if filter in supported_filters]

        ##### Apply filters
        if "elevation" in filters:
            # edits timeseries object in place
            fltrs.elevation_filter(
                self,
                height_range=(
                    filter_params.get("elevation_min_m", 0.0),
                    filter_params.get("elevation_max_m", 8000.0),
                ),
            )

        if "daily_mean" in filters:
            fltrs.daily_mean_merge(self)

        if "MAD" in filters:
            fltrs.mad_filter(self, threshold=filter_params.get("mad_threshold", 5.0))

        if "hampel" in filters:
            fltrs.hampel(self)

        if "rolling_median" in filters:
            fltrs.rolling_median(self)

    def bias_correct(self, platform_key=None, orbit_key=None, group_by="platform_orbit",
                      time_bin="1D", min_overlap=3, priority=None,
                      centroid_distance_warn_km=5.0):
        """
        Cross-calibrate multiple persistent sources within this timeseries
        onto a common datum, then remove each source's estimated constant
        offset from its raw observations (multi-mission/multi-track bias
        correction; ported from the standalone harmonize_and_merge work,
        adapted to operate on a single Timeseries.df in place rather than a
        dict of external per-source series).

        A "source" is one persistent virtual station. By default
        (group_by="platform_orbit") that's identified by (platform_key,
        orbit_key) TOGETHER -- e.g. (platform="icesat2", orbit="30") for
        one ICESat-2 beam, which genuinely is stable across every revisit.
        This is deliberately NOT the same grouping as pass_key: pass_key
        identifies one specific crossing in time (used for along-track
        slope/ADM work in basic_filters.py), while (platform_key, orbit_key)
        is meant to identify the same physical ground track revisited many
        times -- the thing that can carry a stable geoid/retracker bias
        worth estimating and removing before Kalman filtering treats every
        source as directly comparable.

        group_by="platform" ignores orbit_key entirely and groups by
        platform_key alone. Use this when there is no field that's
        genuinely constant across revisits of the same ground track --
        e.g. for Sentinel-3/6, where the available "orbit"-mapped field
        turned out to increment roughly once per repeat cycle (confirmed
        empirically: ~23.6 days/value for a 27-day S3 cycle, ~9.6 days/
        value for a ~9.9-day S6 cycle) rather than staying constant, so
        grouping by it fragments each platform into many near-single-
        observation "sources" that can't be calibrated against anything,
        and the real, visible inter-platform bias never gets corrected.
        platform-level grouping is coarser (loses any genuine within-
        platform inter-track bias) but is what's actually achievable
        without a reliable persistent-track field, and is what fixes an
        inter-mission bias that's otherwise passing straight through.

        If fewer than 2 distinct sources are present (e.g. a river target
        fed only by SWOT Hydrocron today), this is a no-op.

        Sources that never overlap in time with anything else cannot be
        calibrated; their observations are dropped from the corrected
        output, with a warning naming them so you can decide whether to
        anchor one manually (e.g. against an in-situ gauge) instead.

        Parameters
        ----------
        platform_key, orbit_key : str, optional
            Defaults to self.platform_key / self.orbit_key. orbit_key is
            ignored entirely when group_by="platform".
        group_by : {"platform_orbit", "platform"}, optional
            See above. Default "platform_orbit" (unchanged default
            behavior); pass "platform" for Sentinel-3/6 given the finding
            above.
        time_bin : str, optional
            "<N>D"-style window used to line up observations in time across
            sources before comparing them. Default "1D".
        min_overlap : int, optional
            Minimum overlapping time bins required to estimate a source's
            bias against the running combined reference. Default 3.
        priority : list[str], optional
            Preferred anchor source id(s), formatted "{platform}_{orbit}"
            (or just "{platform}" when group_by="platform"), if you trust
            one source's calibration more than the others. Defaults to
            whichever source has the most observations.
        centroid_distance_warn_km : float, optional
            If lat_key/lon_key are available, each source's centroid
            (mean lat/lon of its own observations) is computed and compared
            to the anchor's centroid -- purely as transparency metadata,
            NOT a spatial correction (this method still only compares
            sources temporally). A warning is logged if a source's centroid
            is more than this many km from the anchor's, since the
            estimated bias for that source may be partly real spatial
            signal rather than pure calibration offset. Default 5.0.

        After calling, self.bias_correct_diagnostics holds a dict with
        "anchor", "biases", "centroids", "distance_km_from_anchor", and
        "unanchored", for inspection/export by the caller.
        """
        platform_key = platform_key or self.platform_key
        orbit_key = orbit_key or self.orbit_key

        df = self.df.copy()

        if group_by not in ("platform_orbit", "platform"):
            raise ValueError("group_by must be 'platform_orbit' or 'platform'")

        if platform_key not in df.columns or df[platform_key].isna().all():
            logger.warning(
                "bias_correct: platform_key ('%s') not usable on this "
                "timeseries; skipping (no-op).", platform_key,
            )
            return self

        if group_by == "platform_orbit":
            if orbit_key not in df.columns or df[orbit_key].isna().all():
                logger.warning(
                    "bias_correct: orbit_key ('%s') not usable on this "
                    "timeseries; skipping (no-op). Pass group_by='platform' "
                    "if you want platform-level grouping instead.",
                    orbit_key,
                )
                return self
            source_id = df[platform_key].astype(str) + "_" + df[orbit_key].astype(str)
        else:
            source_id = df[platform_key].astype(str)

        if source_id.nunique() < 2:
            logger.info(
                "bias_correct: only one source present (%s); nothing to "
                "harmonize.", source_id.iloc[0] if len(source_id) else "n/a"
            )
            return self

        df = df.assign(_source=source_id)
        dates = pd.to_datetime(df[self.date_key])

        # Centroid tracking: NOT a spatial correction -- bias_correct still
        # only compares sources temporally, exactly as before. This purely
        # records where each source's observations are actually located, so
        # a caller can tell whether an estimated "bias" might be partly real
        # spatial signal (e.g. wind setup, inflow gradient) rather than pure
        # calibration offset, for a large/elongated reservoir where
        # different missions cross at genuinely different locations.
        lat_key = getattr(self, "lat_key", None)
        lon_key = getattr(self, "lon_key", None)
        have_coords = bool(
            lat_key and lon_key and lat_key in df.columns and lon_key in df.columns
        )
        centroids = {}
        if have_coords:
            for src, g in df.groupby("_source"):
                valid = g[[lat_key, lon_key]].dropna()
                if len(valid):
                    centroids[src] = (
                        float(valid[lat_key].mean()), float(valid[lon_key].mean())
                    )

        n_days = _parse_bin_days(time_bin)
        epoch = pd.Timestamp("1970-01-01")
        bin_id = ((dates - epoch).dt.days // n_days) * n_days
        df = df.assign(_bin=epoch + pd.to_timedelta(bin_id, unit="D"))

        frames = {}
        for src, g in df.groupby("_source"):
            binned = g.groupby("_bin")[self.height_key].median().to_frame("height")
            binned["n_points"] = g.groupby("_bin")[self.height_key].count()
            frames[src] = binned

        # Anchor = source with the most non-empty time bins, a close proxy
        # for "most individual observation dates" given date-based binning.
        anchor = next((s for s in (priority or []) if s in frames), None)
        if anchor is None:
            anchor = max(frames, key=lambda s: len(frames[s]))

        harmonized = {anchor: frames[anchor]}
        biases = {anchor: 0.0}
        remaining = {s: f for s, f in frames.items() if s != anchor}

        def combined_reference():
            allh = pd.concat(
                {s: f["height"] for s, f in harmonized.items()}, axis=1, sort=False
            )
            w = pd.concat(
                {s: f["n_points"] for s, f in harmonized.items()}, axis=1, sort=False
            )
            w = w.where(allh.notna())
            return (allh * w).sum(axis=1) / w.sum(axis=1)

        changed = True
        while remaining and changed:
            changed = False
            ref = combined_reference()
            for src in list(remaining):
                b = remaining[src]
                common = ref.index.intersection(b.index)
                common = common[ref.loc[common].notna() & b.loc[common, "height"].notna()]
                if len(common) >= min_overlap:
                    bias = float((b.loc[common, "height"] - ref.loc[common]).median())
                    biases[src] = bias
                    harmonized[src] = b
                    del remaining[src]
                    changed = True

        if remaining:
            logger.warning(
                "bias_correct: sources %s never overlapped in time with "
                "anything else and could not be calibrated; their "
                "observations are dropped from the corrected output.",
                list(remaining.keys()),
            )

        keep_mask = df["_source"].isin(biases.keys())
        bias_values = df["_source"].map(biases)

        corrected = self.df.loc[keep_mask.values].copy()
        corrected[self.height_key] = (
            corrected[self.height_key].values - bias_values[keep_mask].values
        )

        distance_km_from_anchor = {}
        if have_coords and anchor in centroids:
            anchor_lat, anchor_lon = centroids[anchor]
            for src in biases:
                if src in centroids:
                    src_lat, src_lon = centroids[src]
                    distance_km_from_anchor[src] = float(
                        fltrs._haversine_km(anchor_lat, anchor_lon, src_lat, src_lon)
                    )
            far = {s: d for s, d in distance_km_from_anchor.items()
                   if s != anchor and d > centroid_distance_warn_km}
            if far:
                logger.warning(
                    "bias_correct: source(s) %s are >%.0f km from the anchor "
                    "source's centroid (anchor='%s'). The estimated bias for "
                    "these may be partly real spatial signal (e.g. wind "
                    "setup, inflow gradient) rather than pure calibration "
                    "offset -- treat with caution for large/elongated "
                    "targets. Distances (km): %s",
                    list(far.keys()), centroid_distance_warn_km, anchor, far,
                )

        logger.info("bias_correct: estimated biases (m): %s", biases)
        if centroids:
            logger.info("bias_correct: source centroids (lat, lon): %s", centroids)

        self.bias_correct_diagnostics = {
            "anchor": anchor,
            "biases": biases,
            "centroids": centroids,
            "distance_km_from_anchor": distance_km_from_anchor,
            "unanchored": list(remaining.keys()),
        }

        self.df = corrected.reset_index(drop=True)
        return self

    def merge(
        self,
        save_progress=False,
        dir=".\\merged_progress",
        window_km=None,
        max_theilsen_points=60,
        spatial_correction_model=None,
        bias_time_bin="1D",
        bias_min_overlap=3,
        bias_priority=None,
        bias_group_by="platform_orbit",
        bias_centroid_warn_km=5.0,
        ref_lat=None,
        ref_lon=None,
        distance_penalty_scale_per_km=None,
        svr_linear_err=0.1,
        svr_linear_epsilon=0.1,
        svr_linear_max_iter=5000,
        svr_radial_max_iter=-1,
        svr_radial_err=1.0,
        svr_radial_rbf_c=10000,
        svr_radial_gamma=0.0000438,
        svr_radial_epsilon=0.1,
    ):
        """
        Run the full cleaning/merging pipeline: along-track outlier
        rejection, multi-source bias correction, ADM error estimation,
        Kalman filtering, and a final RBF outlier pass.

        Parameters worth tuning per target type (e.g. river vs. reservoir,
        typically passed in from flow-level code):
            window_km : along-track half-width (km) for the windowed ADM.
                None keeps the original day-grouped-median ADM. Rivers with
                a real slope benefit from a real window (e.g. ~0.5-1.5 km);
                small/round lakes may not need one at all.
            svr_radial_err : final RBF confidence band (m). DAHITI-style
                defaults are ~1.0 for rivers, ~0.1 for lakes/reservoirs.
            bias_time_bin : widen this (e.g. "10D"-"15D") if your sources
                have sparse, non-coincident revisit patterns (e.g. ICESat-2
                vs. a 27-day Sentinel repeat) -- at "1D" they may never be
                seen as overlapping and get dropped by bias_correct.
            bias_group_by : "platform_orbit" groups sources by (platform,
                orbit) -- only meaningful if orbit_key is a genuinely
                stable, repeating ground-track identifier. If it isn't
                (e.g. confirmed for Sentinel-3/6: the available field
                increments roughly once per repeat cycle rather than
                staying constant), use "platform" instead, which groups by
                platform alone and is what actually resolves an inter-
                mission bias rather than fragmenting each platform into
                many single-crossing "sources" that can't be calibrated.
            bias_centroid_warn_km : logs a warning (and records distances
                in self.bias_correct_diagnostics) when a source's
                observations are centered more than this far from the
                anchor source's -- NOT a spatial correction, just a flag
                that the estimated bias for that source may be partly real
                spatial signal rather than pure calibration offset, worth a
                closer look for large/elongated targets. Default 5.0 km.
            ref_lat, ref_lon, distance_penalty_scale_per_km : if all three
                are given, applies apply_distance_penalty after
                daily_mad_error -- inflates each observation's error based
                on distance from (ref_lat, ref_lon), e.g. the reservoir
                polygon's own centroid, so Kalman naturally down-weights
                crossings far from the main reservoir body (e.g. far
                upstream, subject to real slope bias that local ADM alone
                cannot detect since it only measures local scatter, not
                representativeness). None (default) skips this step
                entirely -- opt-in, not a silent behavior change.
            spatial_correction_model : if given (see
                basic_filters.fit_spatial_correction_model), applies
                apply_spatial_correction right after svr_linear -- a
                genuine height correction (not just error inflation) using
                a spatial deviation model fit from a dense source (e.g.
                ICESat-2), letting sparse missions benefit from it too.
                Applied early, before bias_correct, so bias_correct
                estimates pure platform calibration offset rather than
                having spatial sampling differences between missions
                contaminate the bias estimate. None (default) skips this
                step entirely. See flows.py for the fit-once/persist-to-
                disk wrapper -- refitting this model every run would make
                past corrections shift retroactively as new data arrives.
        See basic_filters.svr_linear/daily_mad_error/svr_radial for the
        rest of these parameters.

        IMPORTANT: svr_linear_max_iter (for svr_linear) and svr_radial_max_iter
        (for svr_radial) are deliberately SEPARATE parameters, not shared.
        Confirmed empirically that capping svr_radial's RBF+high-C solver
        produces non-monotonic, unreliable output quality (e.g. max_iter=
        5000 kept 10/5760 points on real data; only unbounded gave the
        correct 2154/5760) -- svr_radial_max_iter must stay -1 (unbounded)
        by default. Do not consolidate these into one shared parameter.
        """
        # make a folder for saving steps of the timeseries cleaning process
        general.ifnotmakedirs(dir)

        # spatial_correction.csv/distance_penalty.csv are only written
        # below if that feature is actually active for THIS run (both
        # are opt-in, off by default) -- but ifnotmakedirs only creates
        # the directory if it's missing, it never clears existing
        # content. Without this, a run where the feature is OFF would
        # simply skip writing a new file, leaving a stale one from
        # whenever it was last ON sitting here looking like current
        # output -- confirmed as a real, reported point of confusion.
        for _stale_candidate in ("spatial_correction.csv", "distance_penalty.csv"):
            _stale_path = os.path.join(dir, _stale_candidate)
            if os.path.exists(_stale_path):
                os.remove(_stale_path)

        # run the SVR linear outlier filtering
        self = fltrs.svr_linear(self, err=svr_linear_err, epsilon=svr_linear_epsilon, max_iter=svr_linear_max_iter)
        if save_progress:
            self.export_csv(os.path.join(dir, "svr_linear.csv"))

        # optional: correct heights using a spatial deviation model (e.g.
        # fit from dense ICESat-2 coverage) BEFORE bias_correct, so the
        # bias estimate reflects pure platform calibration offset rather
        # than being contaminated by different missions sampling
        # spatially different parts of the reservoir
        if spatial_correction_model is not None:
            self = fltrs.apply_spatial_correction(self, spatial_correction_model)
            if save_progress:
                self.export_csv(os.path.join(dir, "spatial_correction.csv"))

        # cross-calibrate multiple sources (missions/tracks) onto a common
        # datum, if more than one is present; a no-op otherwise (e.g. a
        # river target fed only by SWOT Hydrocron today)
        self = self.bias_correct(
            time_bin=bias_time_bin, min_overlap=bias_min_overlap,
            priority=bias_priority, group_by=bias_group_by,
            centroid_distance_warn_km=bias_centroid_warn_km,
        )
        if save_progress:
            self.export_csv(os.path.join(dir, "bias_correct.csv"))

        # run the ADM running filter
        self = fltrs.daily_mad_error(
            self, window_km=window_km, max_theilsen_points=max_theilsen_points
        )
        if save_progress:
            self.export_csv(os.path.join(dir, "daily_mad_error.csv"))

        # optional: inflate error by distance from a reference location
        # (e.g. reservoir centroid), so Kalman down-weights crossings that
        # may be hydraulically unrepresentative (e.g. far upstream) even if
        # ADM alone would treat them as highly precise
        if ref_lat is not None and ref_lon is not None and distance_penalty_scale_per_km is not None:
            self = fltrs.apply_distance_penalty(
                self, ref_lat=ref_lat, ref_lon=ref_lon,
                scale_per_km=distance_penalty_scale_per_km,
            )
            if save_progress:
                self.export_csv(os.path.join(dir, "distance_penalty.csv"))

        # run Kalman filter
        df_kalman = fltrs.kalman(self)
        self = Timeseries(
            df_kalman,
            date_key=self.date_key,
            height_key=self.height_key,
            error_key=self.error_key,
            lat_key=self.lat_key,
            lon_key=self.lon_key,
            pass_key=self.pass_key,
            platform_key=self.platform_key,
            orbit_key=self.orbit_key,
            preset_error_key=self.preset_error_key,
        )
        if save_progress:
            self.export_csv(os.path.join(dir, "kalman.csv"))

        # a radial base svr to get the final timeseries
        df_rbf = fltrs.svr_radial(
            self,
            err=svr_radial_err,
            rbf_c=svr_radial_rbf_c,
            gamma=svr_radial_gamma,
            epsilon=svr_radial_epsilon,
            max_iter=svr_radial_max_iter,
        )
        self = Timeseries(
            df_rbf,
            date_key=self.date_key,
            height_key=self.height_key,
            error_key=self.error_key,
            lat_key=self.lat_key,
            lon_key=self.lon_key,
            pass_key=self.pass_key,
            platform_key=self.platform_key,
            orbit_key=self.orbit_key,
            preset_error_key=self.preset_error_key,
        )
        if save_progress:
            self.export_csv(os.path.join(dir, "svr_radial.csv"))

        return self

    def export_csv(self, path):
        self.df.to_csv(path)


def concat(
    timeseries: list = [],
    main_date_key="date",
    main_height_key="height",
    main_lat_key="lat",
    main_lon_key="lon",
    main_pass_key="pass",
    main_platform_key="platform",
    main_orbit_key="orbit",
    main_preset_error_key="preset_error",
):
    """
    Concatenate several Timeseries objects (typically one per mission/product)
    into a single Timeseries on a common set of column names.

    Different missions may use different underlying column names for the same
    conceptual field (e.g. Sentinel-6's "file_name" vs. ICESat-2's lack of an
    equivalent pass identifier) -- each source's key attributes (ts.pass_key,
    ts.platform_key, etc.) are read and renamed onto the common main_*_key
    names here, rather than assuming a fixed literal column name. A column is
    only carried through for sources that actually have it; sources missing
    a given key simply get NaN for that column after concatenation.

    This includes preset_error_key: a mission that already supplies its own
    formal per-observation uncertainty (e.g. SWOT's wse_u) needs that column
    preserved through concatenation, or it silently disappears the moment
    it's combined with missions that don't have one -- exactly the same
    class of bug as lat/lon/pass being dropped without this handling.
    """
    optional_keys = {
        main_lat_key: "lat_key",
        main_lon_key: "lon_key",
        main_pass_key: "pass_key",
        main_platform_key: "platform_key",
        main_orbit_key: "orbit_key",
        main_preset_error_key: "preset_error_key",
    }

    df_list = []

    # create a single timeseries object from the multiple timeseries
    for ts in timeseries:
        rename_map = {ts.date_key: main_date_key, ts.height_key: main_height_key}
        keep = [ts.date_key, ts.height_key]

        for main_name, attr in optional_keys.items():
            src_col = getattr(ts, attr, None)
            if src_col and src_col in ts.df.columns:
                keep.append(src_col)
                rename_map[src_col] = main_name

        df = ts.df[keep].rename(columns=rename_map)
        df_list.append(df)

    # concatenate dfs
    df = pd.concat(df_list, ignore_index=True)

    # turn this combined df into a timeseries object and clean
    ts = Timeseries(
        df,
        date_key=main_date_key,
        height_key=main_height_key,
        lat_key=main_lat_key,
        lon_key=main_lon_key,
        pass_key=main_pass_key,
        platform_key=main_platform_key,
        orbit_key=main_orbit_key,
        preset_error_key=(
            main_preset_error_key if main_preset_error_key in df.columns else None
        ),
    )

    # return the merged timeseries
    return ts