"""General Utilities to aid in other modules"""

import os
import shutil
import zipfile

def ifnotmakedirs(dir : str):
    if not os.path.exists(dir):
        os.makedirs(dir)

def unzip_dir_files(dir : str, dest_dir:str):

    # ensure the destination exists
    ifnotmakedirs(dest_dir)

    # walk throught the files in the directory and check for zips
    for file in os.listdir(dir):
        if file.endswith('.zip'):

            # unzip and save to a temporary folder in the main dir
            with zipfile.ZipFile(os.path.join(dir, file), 'r') as zip_ref:
                zip_ref.extractall(dest_dir)

def unzip_dir_files_with_ext(dir : str, dest_dir:str, ext:str):

    # ensure the destination exists
    ifnotmakedirs(dest_dir)

    # walk throught the files in the directory and check for zips
    for file in os.listdir(dir):
        if file.endswith('.zip'):

            # unzip and save to a temporary folder in the main dir
            with zipfile.ZipFile(os.path.join(dir, file), 'r') as zip_ref:

                # list files in archive
                members = zip_ref.namelist()
                for member in members:
                    if member.endswith(ext):
                        zip_ref.extract(member, dest_dir)


def remove_non_exts(dir:str, ext:str):

    for name in os.listdir(dir):
        if not name.endswith(ext):

            item = os.path.join(dir, name)

            if os.path.isdir(item):
                shutil.rmtree(item)

            elif os.path.isfile(item):
                os.remove(item)