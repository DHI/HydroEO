from dataclasses import dataclass
import warnings
import pandas as pd

import altimetry.utils.filters.basic_filters as fltrs


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
        supported_filters = [
            "elevation",
            "MAD",
            "daily_mean",
            "hampel",
            "rolling_median",
        ]

        for filter in filters:
            if filter not in supported_filters:
                warnings.warn(
                    f"{filter} is not in list of supported filters and will not be applied"
                )
        filters = [filter for filter in filters if filter in supported_filters]

        ##### Apply filters
        if "elevation" in filters:
            # edits timeseries object in place
            fltrs.elevation_filter(self, height_range=(0, 8000))

        if "daily_mean" in filters:
            fltrs.daily_mean_merge(self)

        if "MAD" in filters:
            fltrs.mad_filter(self, threshold=4)

        if "hampel" in filters:
            fltrs.hampel(self)

        if "rolling_median" in filters:
            fltrs.rolling_median(self)

    def bias_correct(self, orbit_key="orbit", product_key="product"):
        # bias correct between orbits and products
        # TODO: implement bias correction

        # for now pretend like we bias corrected and put out in the correct output
        pass

    def merge(self):
        # run the SVR linear outlier filtering
        print("SVR")
        fltrs.svr_linear(self)

        # run the big outlier filtering on the whole set of values
        print("MAD filter")
        fltrs.mad_filter(self, threshold=5)

        # run the ADM running filter
        print("Daily MAD filter")
        fltrs.daily_mad_error(self)

        # here we should run the kalman filter
        print("Kalman")
        df_kalman = fltrs.kalman(self)
        self = Timeseries(df_kalman, date_key=self.date_key, height_key=self.height_key)

        # a radial base svr to get the final timeseries
        df_rbf = fltrs.svr_radial(self)
        self = Timeseries(df_rbf, date_key=self.date_key, height_key=self.height_key)

    def export_csv(self, path):
        self.df.to_csv(path)


def concat(timeseries: list = [], main_date_key="date", main_height_key="height"):
    df_list = []

    # create a single timeseries object from the multiple timeseries
    for ts in timeseries:
        df = ts.df[[ts.date_key, ts.height_key]]
        df = df.rename(
            columns={ts.date_key: main_date_key, ts.height_key: main_height_key}
        )

        df_list.append(df)

    # concatenate dfs
    df = pd.concat(df_list)

    # turn this combined df into a timeseries object and clean
    ts = Timeseries(df, date_key=main_date_key, height_key=main_height_key)

    # return the merged timeseries
    return ts
