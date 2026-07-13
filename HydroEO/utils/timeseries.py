from dataclasses import dataclass
import logging
import warnings
import os
import pandas as pd

import HydroEO.utils.filters.basic_filters as fltrs
from HydroEO.utils import general

logger = logging.getLogger(__name__)


def _parse_bin_days(time_bin):
    """Parse a simple '<N>D' style string into an integer number of days."""
    s = str(time_bin).strip().upper()
    if s.endswith("D"):
        return int(s[:-1]) if s[:-1] else 1
    raise ValueError("time_bin must be a string like '1D' or '5D'")


@dataclass
class Timeseries:
    df: pd.DataFrame
    date_key: str = "date"
    height_key: str = "height"
    error_key: str = "ADM"
    lat_key: str = "lat"
    lon_key: str = "lon"
    pass_key: str = "pass"
    platform_key: str = "platform"
    orbit_key: str = "orbit"
    preset_error_key: str = None

    def __post_init__(self):
        if self.height_key not in self.df.columns:
            logger.warning(
                "Height column '%s' missing from dataframe", self.height_key
            )
            # error_key is intentionally not checked here: it's populated
            # later in the pipeline (by daily_mad_error), not expected to
            # exist at construction time. Functions that actually need it
            # (daily_mad_error, kalman) will raise a clear error at the
            # point they read it if it's genuinely still missing then.

        if self.date_key not in self.df.columns:
            if isinstance(self.df.index, pd.DatetimeIndex):
                self.df[self.date_key] = self.df.index
            else:
                logger.warning(
                    "date_key '%s' is not in dataframe and index is not a "
                    "DatetimeIndex, so it could not be auto-populated",
                    self.date_key,
                )

    def clean(self, filters: list, filter_params: dict = None):
        filter_params = filter_params or {}

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
            fltrs.elevation_filter(
                self,
                height_range=(
                    filter_params.get("elevation_min_m", 0.0),
                    filter_params.get("elevation_max_m", 8000.0),
                ),
            )

        if "daily_mean" in filters:
            fltrs.daily_mean_merge(self)

        if "MAD" in filters:
            fltrs.mad_filter(self, threshold=filter_params.get("mad_threshold", 5.0))

        if "hampel" in filters:
            fltrs.hampel(self)

        if "rolling_median" in filters:
            fltrs.rolling_median(self)

    def bias_correct(self, orbit_key="orbit", product_key="platform"):
        raise NotImplementedError("bias_correct is not yet implemented")

    def merge(self, save_progress=False, dir=".\\merged_progress"):
        # make a folder for saving steps of the timeseries cleaning process
        general.ifnotmakedirs(dir)

        # run the SVR linear outlier filtering
        self = fltrs.svr_linear(self)
        if save_progress:
            self.export_csv(os.path.join(dir, "svr_linear.csv"))

        # run the ADM running filter
        self = fltrs.daily_mad_error(self)
        if save_progress:
            self.export_csv(os.path.join(dir, "daily_mad_error.csv"))

        # run Kalman filter
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
