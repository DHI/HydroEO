"""Reservoirs: extraction + clean + merge + dfs0-export orchestration.

create_reservoirs_timeseries and everything it (transitively) calls --
_extract_reservoirs_timeseries, its three per-mission workers,
_clean_reservoirs_timeseries, _merge_reservoirs_timeseries, and
_export_cleaned_to_dfs0 -- are tested together via
patch.object(flows, "_name") in tests/unit/test_flows.py and
tests/unit/test_timeseries.py, so they all live in this one module.
`mikeio` is imported here (rather than at the package level only) because
_export_cleaned_to_dfs0 is the only place it's actually called, and tests
patch it as `flows.mikeio` -- see flows/__init__.py's re-export.
"""

import logging
import os

import mikeio
import pandas as pd
from tqdm import tqdm

from HydroEO.satellites import swot, icesat2, sentinel
from HydroEO.utils import general
from ._clean_engine import _clean_timeseries
from ._merge_engine import _merge_timeseries

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from HydroEO.project import Project

logger = logging.getLogger(__name__)


def create_reservoirs_timeseries(prj: "Project") -> None:
    """Extract, clean, and merge timeseries for reservoirs.

    Parameters
    ----------
    prj : Project
        Project instance with reservoirs configuration
    """
    if not hasattr(prj, "reservoirs"):
        return

    # Extract raw observations from downloaded files (skips reservoirs whose
    # gpkg already exists, unless prj.reservoirs.overwrite_extraction=True --
    # confirmed this re-read/re-extraction was the dominant real-world cost,
    # far more than anything in clean()/merge())
    _extract_reservoirs_timeseries(
        prj, overwrite=getattr(prj.reservoirs, "overwrite_extraction", False)
    )

    # Clean observations with filters
    _clean_reservoirs_timeseries(prj)

    # Export to dfs0 if enabled
    if getattr(prj.reservoirs, "export_to_dfs0", False):
        _export_cleaned_to_dfs0(prj)

    # Merge multi-mission timeseries
    _merge_reservoirs_timeseries(prj)


def _extract_reservoirs_timeseries(prj: "Project", overwrite: bool = False) -> None:
    """Extract timeseries observations from raw downloaded files.

    Parameters
    ----------
    overwrite : bool, optional
        If False (default), any reservoir/mission whose output .gpkg
        already exists is skipped entirely rather than re-read and
        re-extracted. Confirmed on real data that this re-extraction --
        not clean()/merge() -- was the dominant cost in real end-to-end
        runs (orders of magnitude larger than the merge pipeline itself).
        Set True to force re-extraction (e.g. new raw downloads arrived).
    """
    if "icesat2" in prj.to_process:
        _extract_icesat2_observations(prj, overwrite=overwrite)

    if "sentinel3" in prj.to_process:
        _extract_sentinel_observations(prj, "sentinel3", "S3", overwrite=overwrite)

    if "sentinel6" in prj.to_process:
        _extract_sentinel_observations(prj, "sentinel6", "S6", overwrite=overwrite)

    if "swot" in prj.to_process:
        _extract_swot_observations(prj, overwrite=overwrite)


def _extract_icesat2_observations(prj: "Project", overwrite: bool = False) -> None:
    """Extract ICESat-2 ATL13 observations for each reservoir.

    ICESat-2's raw data is a single atl13.parquet per reservoir that is
    rewritten in full (the entire configured date range) on every
    download -- unlike SWOT/Sentinel, there's no per-granule file list
    to diff against a "processed" log. Re-extraction is instead gated
    on whether the parquet is newer than the last extraction (e.g.
    after project.update() downloaded through a later date), rather
    than simply "does icesat2.gpkg already exist" -- otherwise data
    added by update() would never reach the merged timeseries without
    overwrite=True forcing a full reprocess of every reservoir.
    """
    available_ids = [
        id
        for id in prj.reservoirs.download_gdf[prj.reservoirs.id_key]
        if os.path.exists(
            os.path.join(prj.dirs["icesat2_processed"], f"{id}", "atl13.parquet")
        )
    ]
    if not available_ids:
        logger.warning("No ICESat-2 downloads found; skipping timeseries extraction.")
        return

    ids_to_process = []
    if overwrite:
        ids_to_process = available_ids
    else:
        skip_count = 0
        for id in available_ids:
            parquet_path = os.path.join(
                prj.dirs["icesat2_processed"], f"{id}", "atl13.parquet"
            )
            dst_path = os.path.join(
                prj.dirs["output"], f"{id}", "raw_observations", "icesat2.gpkg"
            )
            if not os.path.exists(dst_path):
                ids_to_process.append(id)
            elif os.path.getmtime(parquet_path) > os.path.getmtime(dst_path):
                ids_to_process.append(id)
            else:
                skip_count += 1
        if skip_count:
            logger.info(
                "ICESat-2 extraction: skipping %d reservoir(s) with an existing "
                "icesat2.gpkg no older than its source parquet (pass "
                "overwrite=True to force re-extraction).",
                skip_count,
            )
        if not ids_to_process:
            return

    empty_ids = []
    for id in tqdm(ids_to_process, desc="Extracting ICESat-2 ATL13 product"):
        sub_gdf = prj.reservoirs.download_gdf.loc[
            prj.reservoirs.download_gdf[prj.reservoirs.id_key] == id
        ]
        download_dir = os.path.join(prj.dirs["icesat2_processed"], f"{id}")
        dst_dir = os.path.join(prj.dirs["output"], f"{id}", "raw_observations")
        general.ifnotmakedirs(dst_dir)
        dst_path = os.path.join(dst_dir, "icesat2.gpkg")

        try:
            icesat2.extract_observations(
                src_dir=download_dir,
                dst_path=dst_path,
                features=sub_gdf,
                atl13_fields=prj.mission_options.get("icesat2", {}).get("atl13_fields"),
                track_keys=prj.mission_options.get("icesat2", {}).get("track_keys"),
            )
        except Exception as exc:
            logger.warning("Failed to extract ICESat-2 for %s: %s", id, exc)

        if not os.path.exists(dst_path):
            empty_ids.append(id)

    if empty_ids:
        logger.warning(
            "ICESat-2 timeseries empty for: %s (no observations passed the spatial filter or the download returned no data)",
            ", ".join(str(i) for i in empty_ids),
        )


def _extract_sentinel_observations(
    prj: "Project", mission_key: str, product: str, overwrite: bool = False
) -> None:
    """Extract Sentinel-3 or Sentinel-6 observations for each reservoir.

    Uses a per-reservoir "already extracted" file log (see
    HydroEO.satellites.sentinel.preprocess.extract_observations's
    processed_log_path) rather than a per-reservoir skip based on
    whether {mission_key}.gpkg already exists -- so new subset files
    downloaded after the first extraction (e.g. via project.update())
    are picked up and appended on the next run instead of being
    silently ignored until overwrite=True forces a full reprocess of
    every reservoir.
    """
    available_ids = [
        id
        for id in prj.reservoirs.download_gdf[prj.reservoirs.id_key]
        if os.path.exists(os.path.join(prj.dirs[mission_key], f"{id}"))
    ]
    if not available_ids:
        logger.warning(
            "No %s downloads found; skipping timeseries extraction.", mission_key
        )
        return

    empty_ids = []
    for id in tqdm(available_ids, desc=f"Extracting Sentinel-{product} product"):
        sub_gdf = prj.reservoirs.download_gdf.loc[
            prj.reservoirs.download_gdf[prj.reservoirs.id_key] == id
        ]
        download_dir = os.path.join(prj.dirs[mission_key], f"{id}")
        dst_dir = os.path.join(prj.dirs["output"], f"{id}", "raw_observations")
        general.ifnotmakedirs(dst_dir)
        dst_path = os.path.join(dst_dir, f"{mission_key}.gpkg")
        processed_log_path = os.path.join(download_dir, "extracted.log")

        try:
            sentinel.extract_observations(
                src_dir=download_dir,
                dst_path=dst_path,
                features=sub_gdf,
                sigma0_max=prj.mission_options.get(mission_key, {}).get(
                    "sigma0_max", 1e5
                ),
                processed_log_path=processed_log_path,
                overwrite=overwrite,
            )
        except Exception as exc:
            logger.warning("Failed to extract %s for %s: %s", mission_key, id, exc)

        if not os.path.exists(dst_path):
            empty_ids.append(id)

    if empty_ids:
        logger.warning(
            "Sentinel-%s timeseries empty for: %s (no observations passed the spatial filter or the download returned no data)",
            product,
            ", ".join(str(i) for i in empty_ids),
        )


def _extract_swot_observations(prj: "Project", overwrite: bool = False) -> None:
    """Extract SWOT Lake SP observations for all reservoirs.

    Uses a project-wide "already extracted" file log (see
    HydroEO.satellites.swot.preprocess.extract_observations's
    processed_log_path) rather than a per-reservoir skip based on
    whether swot.gpkg already exists. SWOT granule shapefiles
    accumulate forever in one shared download directory across ALL
    reservoirs (see satellites.swot.preprocess.subset_by_id), so
    whether a granule is "new" can only be decided at the file level,
    not per reservoir -- and doing so lets new downloads (e.g. via
    project.update()) get picked up and appended on the next run
    instead of requiring overwrite=True to reprocess every reservoir
    from the full granule history.
    """
    download_dir = prj.dirs["swot"]
    if not os.path.exists(download_dir):
        logger.warning("No SWOT downloads found; skipping timeseries extraction.")
        return

    features = prj.reservoirs.download_gdf
    id_key = prj.reservoirs.id_key
    processed_log_path = os.path.join(download_dir, "extracted.log")

    empty_ids = swot.extract_observations(
        src_dir=download_dir,
        dst_dir=prj.dirs["output"],
        dst_file_name="swot.gpkg",
        features=features,
        id_key=id_key,
        exclude_obs_id_values=prj.mission_options.get("swot", {}).get(
            "exclude_obs_id_values", ["no_data"]
        ),
        processed_log_path=processed_log_path,
        overwrite=overwrite,
    )
    if empty_ids:
        logger.warning(
            "SWOT timeseries empty for: %s (no observations matched the prior lake ID or all were excluded)",
            ", ".join(str(i) for i in empty_ids),
        )


# Per-product mapping from generic Timeseries key attributes to the actual
# column names each mission's extractor writes. Sentinel-3 shares
# Sentinel-6's extractor/schema (same sentinel.extract_observations
# function, see _extract_sentinel_observations).
#
# ****************************************************************************
# TERMINOLOGY TRAP -- read before touching orbit_key/pass_key for Sentinel:
# The raw Sentinel-3/6 data has TWO similarly-named but opposite-meaning
# columns:
#   - "orbit": the absolute revolution counter. Unique on every single
#     crossing, never repeats. USELESS as orbit_key (bias_correct needs a
#     persistent identifier to accumulate overlap against -- grouping by
#     something that's different every time means every "source" has
#     exactly 1 observation and nothing can ever be calibrated: this is
#     exactly the bug that caused every S3A/S3B track to be dropped as
#     unanchored in practice).
#   - "pass": the satellite-engineering term for the STABLE, REPEATING
#     ground track number (same value every ~27-day repeat cycle for
#     S3A/S3B). This is what orbit_key actually needs.
# Confusingly, our own framework's `pass_key` means the OPPOSITE thing (one
# specific, one-time crossing -- e.g. file_name) from what "pass" means in
# the satellite data itself (the repeating track). Do not be tempted to
# point pass_key at the raw "pass" column -- file_name is correct there.
# ****************************************************************************
#
# ICESat-2's orbit_key is "beam" (the persistent ground track/virtual
# station) -- cycle_number only matters as an ingredient of the compound
# "pass" column built at extraction time (see
# HydroEO.satellites.icesat2.preprocess.extract_observations). SWOT's
# LakeSP product is already one integrated WSE per crossing with its own
# formal uncertainty (wse_u), so it needs neither lat/lon nor pass_key --
# see preset_error_key, and daily_mad_error's handling of it.
def _clean_reservoirs_timeseries(prj: "Project") -> None:
    """Apply quality filters to extracted reservoir timeseries."""
    _clean_timeseries(prj, "reservoirs")


def _merge_reservoirs_timeseries(prj: "Project") -> None:
    """Merge multi-mission timeseries into combined datasets, for reservoirs."""
    _merge_timeseries(prj, "reservoirs")


def _export_cleaned_to_dfs0(prj: "Project") -> None:
    """Export cleaned timeseries observations to dfs0 format.

    Parameters
    ----------
    prj : Project
        Project instance with reservoirs configuration

    Notes
    -----
    Exports height data from cleaned CSV observations to dfs0 format
    for each product in each reservoir's cleaned_observations folder.
    """
    ids_with_cleaned = [
        id
        for id in prj.reservoirs.download_gdf[prj.reservoirs.id_key]
        if os.path.exists(
            os.path.join(prj.dirs["output"], f"{id}", "cleaned_observations")
        )
    ]

    if not ids_with_cleaned:
        logger.warning(
            "No cleaned observations found for any reservoir; skipping dfs0 export."
        )
        return

    for id in ids_with_cleaned:
        cleaned_dir = os.path.join(prj.dirs["output"], f"{id}", "cleaned_observations")

        for product in tqdm(
            prj.to_process, desc="Exporting cleaned observations to dfs0"
        ):
            csv_path = os.path.join(cleaned_dir, f"{product}.csv")

            if not os.path.exists(csv_path):
                continue

            try:
                # Load and prepare data
                df = pd.read_csv(csv_path)

                # Parse and set datetime index
                df["date"] = pd.to_datetime(
                    df.date, format="mixed", utc=True
                ).dt.tz_convert(None)
                df = df.set_index("date")
                df = df.sort_index()

                # Remove rows with NaN heights
                df = df.dropna(subset=["height"])

                if len(df) == 0:
                    logger.warning(
                        "No valid height data for %s in %s after cleaning",
                        product,
                        id,
                    )
                    continue

                # Create Dataset and export to dfs0
                items = {"height": mikeio.ItemInfo(mikeio.EUMType.Water_Level)}
                ds = mikeio.from_pandas(df[["height"]], items=items)

                # Write dfs0 file
                dfs0_path = os.path.join(cleaned_dir, f"{product}.dfs0")
                ds.to_dfs(dfs0_path)

                logger.debug("Exported dfs0: %s", dfs0_path)

            except Exception as exc:
                logger.error(
                    "Failed to export %s to dfs0 for %s: %s",
                    product,
                    id,
                    exc,
                )