from dataclasses import dataclass
import warnings
import pandas as pd

from altimetry.utils.filters.basic_filters import sanity
from altimetry.utils.filters.kalman import kalman


@dataclass
class Timeseries:
    df: pd.DataFrame
    date_key: str = "date"
    height_key: str = "height"
    error_key: str = "error"

    def __post__init__(self):
        if (
            self.error_key not in self.df.columns
            or self.height_key not in self.df.columns
        ):
            print("Height or error columns missing")
            return

        if self.date_key not in self.df.columns:
            if isinstance(self.df.index, pd.DatetimeIndex):
                self.df[self.date_key] = self.df.index
            else:
                print(
                    "date key is not in dataframe but index has been set as date column"
                )
                return

    def filter(self, filters: list):
        # Ensure that provided filters are supported
        supported_filters = ["sanity", "kalman"]

        for filter in filters:
            if filter not in supported_filters:
                warnings.warn(
                    f"{filter} is not in list of supported filters and will not be applied"
                )
        filters = [filter for filter in filters if filter in supported_filters]

        if "sanity" in filters:
            df = fltrs.basic_filters.sanity()
