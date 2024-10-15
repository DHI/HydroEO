"""The project module"""

import os
import typing
import yaml

import numpy as np


class Project:

    def __init__(self, name: str, config: str = 0.0) -> None:
        """Create a new simulator

        Parameters
        ----------
        name: str
            name of altimetry project
        config: path to configuration file
            lower bound

        Examples
        --------
        >>> proj = Project(name="", low=2.0, high=10.0)
        """
        self._name = name
        self._config = config

    def initialize(self) -> None:

        pass