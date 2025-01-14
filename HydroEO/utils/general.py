"""General Utilities to aid in other modules"""

import os
import shutil
import zipfile
from typing import Union
from tqdm import tqdm


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
