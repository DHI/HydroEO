"""General Utilities to aid in other modules"""

import os
import typing

def ifnotmakedirs(dir : str):
    if not os.path.exists(dir):
        os.makedirs(dir)
