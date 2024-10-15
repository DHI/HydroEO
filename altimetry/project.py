"""The project module"""

# import os
# import typing
# import yaml


class Project:

    def __init__(self, name: str, config: str = './data/myconfig.yaml') -> None:
        """Create a new simulator

        Parameters
        ----------
        name: str
            name of altimetry project
        config: path to configuration file
            lower bound

        Examples
        --------
        >>> proj = Project(name="Lower Mekong", config='./data/lower_mekong.yaml')
        """
        self._name = name
        self._config = config

    def initialize(self) -> None:

        pass