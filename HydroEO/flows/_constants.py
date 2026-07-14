"""Shared constants for the merge/clean pipeline (both reservoirs and rivers).

Split out of the original flows.py so both _clean_engine.py and
_merge_engine.py (and their reservoir/river wrapper modules) can import
these without needing to import each other.
"""

PRODUCT_TIMESERIES_KEYS = {
    "sentinel3": dict(
        lat_key="lat", lon_key="lon", pass_key="file_name",
        platform_key="platform", orbit_key="relative_orbit",
    ),
    # NOTE: sentinel6 still uses "pass" as orbit_key -- NOT verified to be
    # unstable the way it was for sentinel3 (confirmed empirically: on
    # real data, "pass" was unique-per-crossing for every S3A/S3B visit,
    # i.e. not stable at all, while "relative_orbit" genuinely repeated
    # across multiple visits -- e.g. S3B crossed via 2 distinct stable
    # configurations, with real biases of -0.14m and +0.22m that a
    # platform-only grouping was averaging into one misleading +0.04m).
    # Sentinel-6 may have the same "pass" instability and may also have
    # its own "relative_orbit"-equivalent column, but this hasn't been
    # checked against real S6 data -- don't assume the same fix applies
    # without verifying first.
    "sentinel6": dict(
        lat_key="lat", lon_key="lon", pass_key="file_name",
        platform_key="platform", orbit_key="pass",
    ),
    "icesat2": dict(
        lat_key="lat", lon_key="lon", pass_key="pass",
        platform_key="platform", orbit_key="beam",
    ),
    "swot": dict(
        platform_key="platform", orbit_key="orbit", preset_error_key="wse_u",
    ),
}

# Default .merge() tuning for reservoirs, mirroring the shape of
# processing_options (a project-level dict of pipeline parameters) but
# applied once per reservoir rather than per-product, since .merge() runs
# on the already-combined multi-product timeseries. Override via
# prj.merging_options in project config; falls back to these reservoir-
# appropriate defaults if that attribute isn't set.
DEFAULT_RESERVOIR_MERGING_OPTIONS = {
    "window_km": 1.5,
    "svr_linear_err": 0.1,
    "svr_linear_epsilon": 0.1,
    # Both updated from DAHITI's lake-tuned defaults (err=0.1, gamma=
    # 0.0000438) based on real reservoir data validated this session --
    # the lake-tuned gamma implied a ~151-day smoothing lengthscale, far
    # too coarse for a reservoir with real multi-week transitions (see
    # the svr_radial oversmoothing discussion). err=1.0 (river-like,
    # rather than the stricter lake value) and gamma x50 (~21-day
    # lengthscale instead of ~151 days) let the trend actually track
    # real fast changes instead of rejecting them as if they were noise.
    "svr_radial_err": 1.0,
    "svr_radial_rbf_c": 10000,
    "svr_radial_gamma": 0.0000438 * 50,
    "svr_radial_epsilon": 0.1,
    # Confirmed on real data across two reservoirs: revisit sparsity varies a
    # lot (e.g. one reservoir had icesat2/S3A/S3B visiting only 7/14/13
    # distinct days all year). At "10D"/3, sparse sources can fail to ever
    # find 3 overlapping bins and get dropped as unanchored ENTIRELY (not
    # just trimmed) -- confirmed: this silently dropped 2 of 3 missions
    # (icesat2, S3B) for one real reservoir. "20D"/1 recovered all of it.
    # Widening is monotonically safe against data loss (a wider window can
    # only find equal-or-more overlapping bins, never fewer) -- the
    # tradeoff is a very wide bin could blur real water-level change within
    # the window into the bias estimate; 20D is a modest widening, not an
    # extreme one.
    "bias_time_bin": "20D",
    "bias_min_overlap": 1,
    # Confirmed empirically on real data: "platform_orbit" (using
    # orbit_key -- now sentinel3's verified-stable "relative_orbit"
    # column, see PRODUCT_TIMESERIES_KEYS) reveals genuine within-platform
    # bias heterogeneity that "platform" alone was masking. One real
    # reservoir's S3B crosses via two distinct, independently stable
    # configurations (5 days on one, 8 on the other) with biases of
    # -0.14m and +0.22m respectively -- "platform" grouping averaged
    # these into one misleading +0.04m. Same pattern for ICESat-2's
    # beams (orbit_key="beam"): per-beam biases ranged 0.06-0.18m under
    # "platform_orbit", collapsed to one number under "platform". Total
    # kept-row count was IDENTICAL either way on the reservoir tested
    # (3182/4901) -- this is a precision gain, not a data-loss risk, at
    # least for sentinel3/icesat2. NOTE: sentinel6 still uses "pass" as
    # orbit_key (unverified whether it's stable or has a
    # relative_orbit-equivalent -- see PRODUCT_TIMESERIES_KEYS) -- if
    # it's actually unstable like sentinel3's old "pass" mapping was,
    # "platform_orbit" could fragment sentinel6 into single-crossing
    # sources. Recheck against real sentinel6 data before trusting this
    # default for a project relying heavily on sentinel6.
    "bias_group_by": "platform_orbit",
    # Not a spatial correction -- just flags (and records in
    # ts.bias_correct_diagnostics) when a source's observations are
    # centered far from the anchor's, since for a large/elongated
    # reservoir some of the estimated bias could be real spatial signal.
    # Worth a closer look per-reservoir if this fires, not an error.
    "bias_centroid_warn_km": 5.0,
    # Off by default -- inflates Kalman input error by distance from the
    # reservoir polygon's own centroid, addressing crossings that may be
    # hydraulically unrepresentative (e.g. far upstream, subject to real
    # slope bias) even when ADM alone reports them as highly precise. Set
    # to a real value (m of extra error per km of distance) to enable --
    # the right scale depends on the true magnitude of upstream slope bias
    # for your reservoirs, which needs empirical tuning, not a guessed
    # default.
    "distance_penalty_scale_per_km": None,
    # Off by default -- a genuine height correction (not just error
    # inflation) using a spatial deviation model fit once from a dense
    # source (default ICESat-2) and persisted to disk per reservoir (see
    # _get_or_fit_spatial_correction_model) so past corrections don't
    # shift retroactively as new data arrives. Turn on once you've
    # confirmed (as we did empirically) that the target reservoir shows a
    # real, day-to-day-consistent spatial deviation pattern -- fitting
    # requires several qualifying dense-source days (see
    # fit_spatial_correction_model's min_days), and silently does nothing
    # if there isn't enough dense-source data yet.
    "use_spatial_correction": False,
    "spatial_correction_dense_source": "icesat2",
}


DEFAULT_RIVER_MERGING_OPTIONS = {
    # Mostly a starting point copied from the reservoir defaults and NOT
    # independently validated against real river data the way the
    # reservoir defaults were validated this session -- river dynamics
    # differ genuinely (e.g. a real, expected along-reach gradient), so
    # do not assume the rest of these are correct without checking.
    # svr_radial_err/gamma below ARE an explicit exception (set directly,
    # not copied): gamma x100 (~15-day lengthscale, vs DAHITI's ~151-day
    # lake value) and err=1.0, matching the same oversmoothing reasoning
    # as the reservoir defaults, just with a shorter lengthscale given
    # rivers can change faster still.
    "window_km": 1.5,
    "svr_linear_err": 0.1,
    "svr_linear_epsilon": 0.1,
    "svr_radial_err": 1.0,
    "svr_radial_rbf_c": 10000,
    "svr_radial_gamma": 0.0000438 * 100,
    "svr_radial_epsilon": 0.1,
    "bias_time_bin": "20D",
    "bias_min_overlap": 1,
    # Same reasoning/evidence as the reservoir default (see
    # DEFAULT_RESERVOIR_MERGING_OPTIONS) for switching from "platform" to
    # "platform_orbit" -- but this is carried over, not independently
    # verified against real river data. A single river target (node/reach)
    # is a much smaller footprint than a reservoir, so it's genuinely
    # unclear whether the same within-platform configuration split
    # (e.g. S3B's two distinct crossing geometries) would even occur at
    # this scale -- check real per-target bias diagnostics once river
    # data exists before trusting this.
    "bias_group_by": "platform_orbit",
    "bias_centroid_warn_km": 5.0,
    # Off by default, same reasoning as reservoirs. NOTE: an earlier
    # version of this comment claimed a river target's footprint is
    # "much smaller than a reservoir" -- that's wrong for reaches
    # specifically (confirmed ~10km typical length, comparable to or
    # larger than many reservoirs), so distance_penalty/spatial
    # correction may matter just as much for reaches as for reservoirs.
    # It remains true that these tools address spread WITHIN one
    # target's own crossing footprint, never the natural gradient
    # BETWEEN different targets, which should never be "corrected away".
    "distance_penalty_scale_per_km": None,
    "use_spatial_correction": False,
    "spatial_correction_dense_source": "icesat2",
    # Off by default. ONLY meaningful when
    # prj.rivers.target_id_col == "reach_id" -- reference-corrects
    # non-SWOT crossings (ICESat-2/Sentinel-3/6) to what they'd read at
    # the reach's geometric midpoint, using SWOT's own directly-measured
    # "slope" field (see _fit_reach_slope_correction/
    # _apply_reach_slope_correction). Requires "slope" to be present in
    # mission_options["swot"]["hydrocron_fields"]["reaches"]. The
    # midpoint-referenced assumption for SWOT's own reach WSE is an
    # evidence-based inference from the RiverSP processing chain, not a
    # fact directly confirmed in SWOT's documentation -- and the sign of
    # the correction has not been empirically verified against real
    # data in this session. Validate both before trusting this in
    # production.
    "use_reach_slope_correction": False,
}


