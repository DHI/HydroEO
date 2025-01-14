from dataclasses import dataclass
import warnings
import os
import pandas as pd

import HydroEO.utils.filters.basic_filters as fltrs
from HydroEO.utils import general


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
            fltrs.mad_filter(self, threshold=5)

        if "hampel" in filters:
            fltrs.hampel(self)

        if "rolling_median" in filters:
            fltrs.rolling_median(self)

    def bias_correct(self, orbit_key="orbit", product_key="platform"):
        # bias correct between orbits and products
        # TODO: implement bias correction

        # we treat each orbit as a single entity regardless of platform
        # entity = (
        #     self.df[product_key].astype(str)
        #     + "_"
        #     + self.df[orbit_key].astype(int).astype(str)
        # )

        # # first find the orbit with the most overlap with other observations
        # for e in entity.unique():
        #     print(e)

        # for platform in bak.vs_id.unique():
        #     subset = bak.loc[bak.vs_id == platform]
        #     orb_mean = subset.groupby(['orbit']).mean().height.reset_index()

        #     for orb in subset.orbit.unique():
        #         subset_orb = bak.loc[bak.orbit == orb]
        #         orb_mean = subset_orb.groupby(['orbit']).mean().height.reset_index()
        #         # print(orb, [[n, [pd.to_datetime(d).month for d in subset_orb.date.unique()].count(n) / len(subset_orb) * 100] for n in range(1,13)])
        #         print(platform, orb, len(subset_orb), [[n, [pd.to_datetime(d).year for d in subset_orb.date.unique()].count(n) / len(subset_orb) * 100] for n in range(2017, 2023)])

        #     nb_obs = subset.groupby(['orbit']).count().height.reset_index()
        #     ind = nb_obs.idxmax(axis=0)['height']
        #     ref_orbit = nb_obs.loc[ind, 'orbit']

        # for now pretend like we bias corrected and put out in the correct output
        pass

    def merge(self, save_progress=False, dir=".\\merged_progress"):
        # make a folder for saving steps of the timeseries cleaning process
        general.ifnotmakedirs(dir)

        # run the SVR linear outlier filtering
        self = fltrs.svr_linear(self)
        if save_progress:
            self.export_csv(os.path.join(dir, "svr_linear.csv"))

        # run the big outlier filtering on the whole set of values
        # fltrs.mad_filter(self, threshold=6)
        # if save_progress:
        #     self.export_csv(os.path.join(dir, "mad_filter.csv"))

        # run the ADM running filter
        self = fltrs.daily_mad_error(self)
        if save_progress:
            self.export_csv(os.path.join(dir, "daily_mad_error.csv"))

        # here we should run the kalman filter TODO: edit the kalman and svr_radial functions to edit the timeseries in place or edit the above to return new edited timeseries
        df_kalman = fltrs.kalman(self)
        self = Timeseries(df_kalman, date_key=self.date_key, height_key=self.height_key)
        if save_progress:
            self.export_csv(os.path.join(dir, "kalman.csv"))

        # a radial base svr to get the final timeseries
        df_rbf = fltrs.svr_radial(self)
        self = Timeseries(df_rbf, date_key=self.date_key, height_key=self.height_key)
        if save_progress:
            self.export_csv(os.path.join(dir, "svr_radial.csv"))

        return self

    def export_csv(self, path):
        self.df.to_csv(path)


def concat(timeseries: list = [], main_date_key="date", main_height_key="height"):
    df_list = []

    # create a single timeseries object from the multiple timeseries
    for ts in timeseries:
        keep = [ts.date_key, ts.height_key]
        if "orbit" in ts.df.columns:
            keep = keep + ["orbit"]
        if "platform" in ts.df.columns:
            keep = keep + ["platform"]

        df = ts.df[keep]
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
