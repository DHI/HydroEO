from dataclasses import dataclass
import warnings
import pandas as pd

from altimetry.utils.filters.basic_filters import elevation_filter, daily_mean_filter
# from altimetry.utils.filters.kalman import kalman


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

    def clean(self, filters: list):
        ##### Ensure that provided filters are supported
        supported_filters = ["elevation", "daily_mean"]

        for filter in filters:
            if filter not in supported_filters:
                warnings.warn(
                    f"{filter} is not in list of supported filters and will not be applied"
                )
        filters = [filter for filter in filters if filter in supported_filters]

        ##### Apply filters
        if "elevation" in filters:
            # edits timeseries object in place
            elevation_filter(self, height_range=(0, 8000))

        if "daily_mean" in filters:
            daily_mean_filter(self)

    def export_csv(self, path):
        self.df.to_csv(path)


# def merge(timeseries :list = []):

#     df_list = []

#     # create a single timeseries object from the multiple timeseries
#     for ts in timeseries:
