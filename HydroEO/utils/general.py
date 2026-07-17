"""General Utilities to aid in other modules"""

import os
import re
import shutil
import zipfile
from typing import Union
from tqdm import tqdm


def normalize_path(path_str: str) -> str:
    """Make a path from a config file robust to mixed / and \ separators,
    regardless of what OS it was written on or is being read on."""
    path_str = path_str.strip().strip('"').strip("'")  # strip stray quotes too
    unified = re.sub(r'[\\/]+', '/', path_str)          # collapse any run of slashes to a single /
    return os.path.normpath(unified)                    # let os.path convert to native separators



def ifnotmakedirs(dir: str):
    if not os.path.exists(dir):
        os.makedirs(dir)


def unzip_dir_files(dir: str, dest_dir: str, show_progress=False):
    # ensure the destination exists
    ifnotmakedirs(dest_dir)

    # walk throught the files in the directory and check for zips
    for file in tqdm(
        os.listdir(dir),
        desc=f"Unzipping files in {os.path.basename(dir)}",
        unit="file",
        disable=not show_progress,
    ):
        if file.endswith(".zip"):
            # unzip and save to a temporary folder in the main dir
            with zipfile.ZipFile(os.path.join(dir, file), "r") as zip_ref:
                zip_ref.extractall(dest_dir)


def unzip_dir_files_with_ext(dir: str, dest_dir: str, ext: str, show_progress=False):
    # push all ext inputs to a list of inputs
    if isinstance(ext, str):
        ext = [ext]

    # ensure the destination exists
    ifnotmakedirs(dest_dir)

    # walk throught the files in the directory and check for zips
    for file in tqdm(
        os.listdir(dir),
        desc=f"Unzipping files in {os.path.basename(dir)}",
        unit="file",
        disable=not show_progress,
    ):
        if file.endswith(".zip"):
            # unzip and save to a temporary folder in the main dir
            with zipfile.ZipFile(os.path.join(dir, file), "r") as zip_ref:
                # list files in archive
                members = zip_ref.namelist()
                for member in members:
                    # see if the file is not in the list of extentions to keep, and delete if not
                    member_ext = "." + member.split(".")[-1]
                    if member_ext in ext:
                        zip_ref.extract(member, dest_dir)


def remove_non_exts(dir: str, ext: Union[str, list]):
    # push all inputs to a list of inputs
    if isinstance(ext, str):
        ext = [ext]

    # cycle through all files in directory
    for name in os.listdir(dir):
        # check if file extention is not in list
        name_ext = "." + name.split(".")[-1]
        if name_ext not in ext:
            # determine if its a file or directory and delete accordingly
            item = os.path.join(dir, name)

            if os.path.isdir(item):
                shutil.rmtree(item)

            elif os.path.isfile(item):
                os.remove(item)


def read_id_log(log_path: str) -> set:
    """Read a newline-delimited log of identifiers into a set.

    Matches the `downloaded.log` convention already used by
    HydroEO.satellites.swot._download.download and
    HydroEO.satellites.sentinel.download to track which granules have
    already been fetched. Used the same way here to track which raw
    files have already been read during extraction (see
    HydroEO.satellites.swot.preprocess.extract_observations and
    HydroEO.satellites.sentinel.preprocess.extract_observations).

    Returns an empty set if the log doesn't exist yet.
    """
    if not os.path.exists(log_path):
        return set()
    with open(log_path, "r") as f:
        return {line.rstrip() for line in f if line.strip()}


def append_id_log(log_path: str, ids) -> None:
    """Append identifiers to a newline-delimited log, creating it if needed."""
    ids = list(ids)
    if not ids:
        return
    log_dir = os.path.dirname(log_path)
    if log_dir:
        ifnotmakedirs(log_dir)
    with open(log_path, "a") as f:
        for i in ids:
            f.write(f"{i}\n")


def write_id_log(log_path: str, ids) -> None:
    """Replace a newline-delimited log wholesale with the given identifiers.

    Used after a forced full re-extraction (overwrite=True), so the log
    reflects exactly the files that were just (re)read, rather than
    keeping stale entries from before the rebuild or leaving the log
    partially out of sync with what's actually on disk.
    """
    log_dir = os.path.dirname(log_path)
    if log_dir:
        ifnotmakedirs(log_dir)
    with open(log_path, "w") as f:
        for i in ids:
            f.write(f"{i}\n")


def append_and_dedupe_gpkg(dst_path: str, new_gdf, subset=None):
    """Merge newly-extracted observations into an existing GeoPackage.

    If `dst_path` doesn't exist yet, `new_gdf` is written as-is. If it
    does exist, the existing contents are read, concatenated with
    `new_gdf`, de-duplicated, and written back -- so repeated
    incremental extraction runs (see satellites.swot.preprocess and
    satellites.sentinel.preprocess extract_observations) accumulate
    history instead of losing whatever was previously extracted.

    Parameters
    ----------
    subset : list[str], optional
        Non-geometry columns to de-duplicate on (e.g. ["obs_id"] for a
        stable per-observation identifier). Defaults to every
        non-geometry column when no such identifier is available --
        rows are only considered duplicates if they match everywhere
        except geometry, which is a safe (if slightly conservative)
        fallback since genuinely new observations essentially never
        collide on every other field by chance.

    Returns
    -------
    GeoDataFrame
        The combined, de-duplicated result that was written to disk.
    """
    import geopandas as gpd
    import pandas as pd

    if os.path.exists(dst_path):
        existing = gpd.read_file(dst_path)
        combined = pd.concat([existing, new_gdf], ignore_index=True)
        dedup_cols = subset or [c for c in combined.columns if c != "geometry"]
        combined = combined.drop_duplicates(subset=dedup_cols, keep="last").reset_index(
            drop=True
        )
        combined = gpd.GeoDataFrame(combined, geometry="geometry", crs=new_gdf.crs)
    else:
        combined = new_gdf

    combined.to_file(dst_path, driver="GPKG")
    return combined


def center_longitude(lon_org):
    # This assumes the longitude is provided in degrees east of the greenwhich meridian

    lon_converted = lon_org[lon_org > 180]
    lon_converted = -(
        360 - lon_converted
    )  # only edit the values greater than 180 degrees

    lon_org[lon_org > 180] = (
        lon_converted  # insert convereted values to the original array
    )

    return lon_org
