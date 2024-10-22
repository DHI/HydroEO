"""General Utilities to aid in other modules"""

import os
import zipfile

def ifnotmakedirs(dir : str):
    if not os.path.exists(dir):
        os.makedirs(dir)

def unzip_dir_files(dir : str, dest_dir:str):

    # ensure the destination exists
    ifnotmakedirs(dest_dir)

    # walk throught he files in the directory and check for zips
    for file in os.listdir(dir):
        if file.split('.')[-1] == 'zip':

            # unzip and save to a temporary folder in the main dir
            with zipfile.ZipFile(os.path.join(dir, file), 'r') as zip_ref:
                zip_ref.extractall(dest_dir)