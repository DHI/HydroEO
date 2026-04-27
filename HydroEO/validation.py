"""Configuration validation for HydroEO projects."""

import os
import datetime
import logging

from HydroEO.satellites.icesat2 import (
    ATL13_DEFAULT_FIELDS,
    SR_ATL13_VALID_ANCILLARY_FIELDS,
)

logger = logging.getLogger(__name__)

SUPPORTED_TRACK_KEYS = ["gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"]
SUPPORTED_CLEAN_FILTERS = ["elevation", "MAD", "daily_mean", "hampel", "rolling_median"]

DEFAULT_SWOT_HYDROCRON_FIELDS = {
    "nodes": [
        "node_id",
        "node_q",
        "reach_id",
        "time_str",
        "wse",
        "wse_u",
        "p_wse",
        "geoid_hght",
        "sword_version",
        "solid_tide",
        "load_tidef",
        "pole_tide",
        "width",
        "width_u",
        "p_width",
        "xovr_cal_q",
        "rdr_sig0",
        "xovr_cal_c",
        "dark_frac",
    ],
    "reaches": [
        "reach_id",
        "reach_q",
        "time_str",
        "wse",
        "wse_u",
        "slope",
        "slope_u",
        "slope2",
        "slope2_u",
        "width",
        "width_u",
        "geoid_hght",
        "solid_tide",
        "load_tidef",
        "pole_tide",
        "p_wse",
        "p_width",
    ],
}

DEFAULT_SWOT_QUALITY_FILTERS = {
    "nodes": {"max_q": 2},
    "reaches": {"max_q": 2},
}


def is_valid_date_tuple(value):
    """Validate that value is a [year, month, day] list with valid date values."""
    if not isinstance(value, list) or len(value) != 3:
        return False
    if not all(isinstance(v, int) for v in value):
        return False
    try:
        datetime.date(*value)
    except Exception:
        return False
    return True


def validate_config(
    config,
    supported_clean_filters=None,
    default_swot_hydrocron_fields=None,
    default_swot_quality_filters=None,
):
    """Validate loaded config and report all discovered issues at once.

    Parameters
    ----------
    config : dict
        The configuration dictionary to validate.
    supported_clean_filters : list, optional
        List of supported cleaning filters. Defaults to SUPPORTED_CLEAN_FILTERS.
    default_swot_hydrocron_fields : dict, optional
        Default SWOT Hydrocron fields. Defaults to DEFAULT_SWOT_HYDROCRON_FIELDS.
    default_swot_quality_filters : dict, optional
        Default SWOT quality filters. Defaults to DEFAULT_SWOT_QUALITY_FILTERS.

    Returns
    -------
    bool
        True if config is valid.

    Raises
    ------
    ValueError
        If config is invalid, raises with detailed list of issues.
    """
    if supported_clean_filters is None:
        supported_clean_filters = SUPPORTED_CLEAN_FILTERS
    if default_swot_hydrocron_fields is None:
        default_swot_hydrocron_fields = DEFAULT_SWOT_HYDROCRON_FIELDS
    if default_swot_quality_filters is None:
        default_swot_quality_filters = DEFAULT_SWOT_QUALITY_FILTERS

    issues = []
    cfg = config or {}

    if "project" not in cfg or not isinstance(cfg["project"], dict):
        issues.append("Missing required section 'project'.")
    elif not cfg["project"].get("main_dir"):
        issues.append("Missing required key 'project.main_dir'.")

    has_reservoirs = "reservoirs" in cfg
    has_rivers = "rivers" in cfg

    if has_reservoirs and has_rivers:
        issues.append(
            "Sections 'reservoirs' and 'rivers' are mutually exclusive. Configure only one."
        )
    if not has_reservoirs and not has_rivers:
        issues.append(
            "Missing required section: provide either 'reservoirs' or 'rivers'."
        )

    if has_reservoirs:
        if not isinstance(cfg["reservoirs"], dict):
            issues.append("Section 'reservoirs' must be a mapping of key/value pairs.")
        else:
            if not cfg["reservoirs"].get("path"):
                issues.append("Missing required key 'reservoirs.path'.")
            elif not os.path.exists(cfg["reservoirs"]["path"]):
                issues.append(
                    f"Path in 'reservoirs.path' does not exist: {cfg['reservoirs']['path']}"
                )

            if not cfg["reservoirs"].get("id_key"):
                issues.append("Missing required key 'reservoirs.id_key'.")

    if has_rivers:
        if not isinstance(cfg["rivers"], dict):
            issues.append("Section 'rivers' must be a mapping of key/value pairs.")
        else:
            rivers_cfg = cfg["rivers"]
            has_aoi_path = bool(rivers_cfg.get("aoi_path"))
            has_node_numbers = "node_numbers" in rivers_cfg
            has_reach_numbers = "reach_numbers" in rivers_cfg
            provided_inputs = sum([has_aoi_path, has_node_numbers, has_reach_numbers])

            if provided_inputs == 0:
                issues.append(
                    "Provide exactly one rivers input source: 'rivers.aoi_path' or 'rivers.node_numbers' or 'rivers.reach_numbers'."
                )
            elif provided_inputs > 1:
                issues.append(
                    "'rivers.aoi_path', 'rivers.node_numbers', and 'rivers.reach_numbers' are mutually exclusive. Provide only one."
                )

            if has_aoi_path:
                path = rivers_cfg["aoi_path"]
                if not os.path.exists(path):
                    issues.append(f"Path in 'rivers.aoi_path' does not exist: {path}")
                elif not path.lower().endswith((".shp", ".gpkg")):
                    issues.append(
                        "'rivers.aoi_path' must reference a '.shp' or '.gpkg' file."
                    )

                if not rivers_cfg.get("id_key"):
                    issues.append(
                        "Missing required key 'rivers.id_key' when 'rivers.aoi_path' is provided."
                    )

                continent_key = rivers_cfg.get("continent_key")
                if continent_key not in ["af", "as", "eu", "na", "oc", "sa"]:
                    issues.append(
                        "'rivers.continent_key' is required with 'rivers.aoi_path' and must be one of ['af', 'as', 'eu', 'na', 'oc', 'sa']."
                    )

                feature_type = rivers_cfg.get("feature_type")
                if feature_type not in ["nodes", "reaches"]:
                    issues.append(
                        "'rivers.feature_type' is required with 'rivers.aoi_path' and must be one of ['nodes', 'reaches']."
                    )

                buffer_meters = rivers_cfg.get("buffer_meters")
                if buffer_meters is not None and (
                    not isinstance(buffer_meters, (int, float)) or buffer_meters < 0
                ):
                    issues.append(
                        "'rivers.buffer_meters' must be None, 0, or a positive number."
                    )

            if has_node_numbers:
                node_numbers = rivers_cfg.get("node_numbers")
                if (
                    not isinstance(node_numbers, list)
                    or len(node_numbers) == 0
                    or any(not isinstance(v, int) for v in node_numbers)
                ):
                    issues.append(
                        "'rivers.node_numbers' must be a non-empty list of integers."
                    )

                if not rivers_cfg.get("id"):
                    issues.append(
                        "Missing required key 'rivers.id' when 'rivers.node_numbers' is provided."
                    )

            if has_reach_numbers:
                reach_numbers = rivers_cfg.get("reach_numbers")
                if (
                    not isinstance(reach_numbers, list)
                    or len(reach_numbers) == 0
                    or any(not isinstance(v, int) for v in reach_numbers)
                ):
                    issues.append(
                        "'rivers.reach_numbers' must be a non-empty list of integers."
                    )

                if not rivers_cfg.get("id"):
                    issues.append(
                        "Missing required key 'rivers.id' when 'rivers.reach_numbers' is provided."
                    )

    for mission in ["swot", "icesat2", "sentinel3", "sentinel6"]:
        if mission not in cfg:
            continue

        mission_cfg = cfg[mission]
        if not isinstance(mission_cfg, dict):
            issues.append(f"Section '{mission}' must be a mapping of key/value pairs.")
            continue

        for key in ["download", "process"]:
            if key not in mission_cfg:
                issues.append(f"Missing required key '{mission}.{key}'.")
            elif not isinstance(mission_cfg[key], bool):
                issues.append(f"'{mission}.{key}' must be a boolean value.")

        for key in ["startdate", "enddate"]:
            if key not in mission_cfg:
                issues.append(f"Missing required key '{mission}.{key}'.")
            elif not is_valid_date_tuple(mission_cfg[key]):
                issues.append(
                    f"'{mission}.{key}' must be [year, month, day] with valid integer values."
                )

        if is_valid_date_tuple(mission_cfg.get("startdate")) and is_valid_date_tuple(
            mission_cfg.get("enddate")
        ):
            if datetime.date(*mission_cfg["startdate"]) > datetime.date(
                *mission_cfg["enddate"]
            ):
                issues.append(
                    f"'{mission}.startdate' cannot be after '{mission}.enddate'."
                )

        filters = mission_cfg.get("processing_filters", ["elevation", "MAD"])
        if not isinstance(filters, list) or any(
            not isinstance(v, str) for v in filters
        ):
            issues.append(f"'{mission}.processing_filters' must be a list of strings.")
        else:
            invalid_filters = [v for v in filters if v not in supported_clean_filters]
            if invalid_filters:
                issues.append(
                    f"Invalid filters in '{mission}.processing_filters': {invalid_filters}. "
                    f"Valid filters are: {supported_clean_filters}."
                )

        min_h = mission_cfg.get("elevation_min_m", 0.0)
        max_h = mission_cfg.get("elevation_max_m", 8000.0)
        if not isinstance(min_h, (int, float)) or not isinstance(max_h, (int, float)):
            issues.append(
                f"'{mission}.elevation_min_m' and '{mission}.elevation_max_m' must be numeric."
            )
        elif min_h >= max_h:
            issues.append(
                f"'{mission}.elevation_min_m' must be smaller than '{mission}.elevation_max_m'."
            )

        mad_threshold = mission_cfg.get("mad_threshold", 5.0)
        if not isinstance(mad_threshold, (int, float)) or mad_threshold <= 0:
            issues.append(f"'{mission}.mad_threshold' must be a positive number.")

        if mission == "swot":
            max_distance = mission_cfg.get("pld_match_max_distance_m", 100.0)
            if not isinstance(max_distance, (int, float)) or max_distance < 0:
                issues.append(
                    "'swot.pld_match_max_distance_m' must be a non-negative number."
                )

            excluded_obs = mission_cfg.get("exclude_obs_id_values", ["no_data"])
            if not isinstance(excluded_obs, list) or any(
                not isinstance(v, str) for v in excluded_obs
            ):
                issues.append("'swot.exclude_obs_id_values' must be a list of strings.")

            hydrocron_fields = mission_cfg.get(
                "hydrocron_fields", default_swot_hydrocron_fields
            )
            if not isinstance(hydrocron_fields, dict):
                issues.append("'swot.hydrocron_fields' must be a mapping.")
            else:
                extra_keys = sorted(
                    key
                    for key in hydrocron_fields.keys()
                    if key not in ["nodes", "reaches"]
                )
                if extra_keys:
                    issues.append(
                        "'swot.hydrocron_fields' only accepts 'nodes' and 'reaches' keys. "
                        f"Unexpected keys: {extra_keys}."
                    )

                for feature_type in ["nodes", "reaches"]:
                    if feature_type not in hydrocron_fields:
                        issues.append(
                            f"Missing required key 'swot.hydrocron_fields.{feature_type}'."
                        )
                        continue

                    fields = hydrocron_fields[feature_type]
                    if not isinstance(fields, list) or any(
                        not isinstance(value, str) for value in fields
                    ):
                        issues.append(
                            f"'swot.hydrocron_fields.{feature_type}' must be a list of strings."
                        )

            quality_filters = mission_cfg.get(
                "quality_filters", default_swot_quality_filters
            )
            if not isinstance(quality_filters, dict):
                issues.append("'swot.quality_filters' must be a mapping.")
            else:
                extra_keys = sorted(
                    key
                    for key in quality_filters.keys()
                    if key not in ["nodes", "reaches"]
                )
                if extra_keys:
                    issues.append(
                        "'swot.quality_filters' only accepts 'nodes' and 'reaches' keys. "
                        f"Unexpected keys: {extra_keys}."
                    )

                for feature_type in ["nodes", "reaches"]:
                    if feature_type not in quality_filters:
                        issues.append(
                            f"Missing required key 'swot.quality_filters.{feature_type}'."
                        )
                        continue

                    feature_filters = quality_filters[feature_type]
                    if not isinstance(feature_filters, dict):
                        issues.append(
                            f"'swot.quality_filters.{feature_type}' must be a mapping."
                        )
                        continue

                    max_q = feature_filters.get("max_q")
                    if not isinstance(max_q, int):
                        issues.append(
                            f"'swot.quality_filters.{feature_type}.max_q' must be an integer."
                        )

        if mission == "icesat2":
            track_keys = mission_cfg.get("track_keys", SUPPORTED_TRACK_KEYS)
            if not isinstance(track_keys, list) or any(
                not isinstance(v, str) for v in track_keys
            ):
                issues.append("'icesat2.track_keys' must be a list of strings.")
            else:
                invalid_track_keys = [
                    k for k in track_keys if k not in SUPPORTED_TRACK_KEYS
                ]
                if invalid_track_keys:
                    issues.append(
                        f"Invalid entries in 'icesat2.track_keys': {invalid_track_keys}. "
                        f"Valid options are: {SUPPORTED_TRACK_KEYS}."
                    )

            fields = mission_cfg.get("atl13_fields", ATL13_DEFAULT_FIELDS)
            if not isinstance(fields, list) or any(
                not isinstance(v, str) for v in fields
            ):
                issues.append("'icesat2.atl13_fields' must be a list of strings.")
            else:
                invalid_fields = [
                    v for v in fields if v not in SR_ATL13_VALID_ANCILLARY_FIELDS
                ]
                if invalid_fields:
                    valid_fields = sorted(SR_ATL13_VALID_ANCILLARY_FIELDS)
                    issues.append(
                        "Invalid ATL13 fields in 'icesat2.atl13_fields': "
                        f"{invalid_fields}. Valid options are: {valid_fields}."
                    )

            atl13_cfg = mission_cfg.get("atl13", {})
            if not isinstance(atl13_cfg, dict):
                issues.append("'icesat2.atl13' must be a mapping of key/value pairs.")
            else:
                if "pass_invalid" in atl13_cfg and not isinstance(
                    atl13_cfg["pass_invalid"], bool
                ):
                    issues.append("'icesat2.atl13.pass_invalid' must be a boolean.")
                for _beam_key in ("beams", "spots"):
                    if _beam_key in atl13_cfg:
                        _v = atl13_cfg[_beam_key]
                        if not isinstance(_v, list) or any(
                            not isinstance(x, str) for x in _v
                        ):
                            issues.append(
                                f"'icesat2.atl13.{_beam_key}' must be a list of strings."
                            )

        if mission in ["sentinel3", "sentinel6"]:
            sigma0_max = mission_cfg.get("sigma0_max", 1e5)
            if not isinstance(sigma0_max, (int, float)) or sigma0_max <= 0:
                issues.append(f"'{mission}.sigma0_max' must be a positive number.")

            subset_file_id = mission_cfg.get(
                "subset_file_id", "enhanced_measurement.nc"
            )
            if not isinstance(subset_file_id, str) or subset_file_id.strip() == "":
                issues.append(f"'{mission}.subset_file_id' must be a non-empty string.")

            download_threads = mission_cfg.get("download_threads", 1)
            if not isinstance(download_threads, int) or download_threads <= 0:
                issues.append(
                    f"'{mission}.download_threads' must be a positive integer."
                )

    if issues:
        raise ValueError("Invalid configuration:\n - " + "\n - ".join(issues))

    return True
