"""simple filters that can be applied to sat timeseries objects"""

import numpy as np
import pandas as pd


def elevation_filter(timeseries, height_range):
    timeseries.df = timeseries.df.loc[timeseries.df[timeseries.height_key] > 0]
    timeseries.df = timeseries.df.loc[timeseries.df[timeseries.height_key] < 8000]

    timeseries.df = timeseries.df.reset_index(drop=True)
    timeseries.df.sort_values(by=timeseries.date_key)


def mad_filter(timeseries, threshold=5):
    # calculate support to make sure we can make a statistical decision, otheriwse leave
    if len(timeseries.df) >= 30:
        # calculate standard deviation and remove obvious outliers
        med = np.median(timeseries.df.height)
        abs_dev = np.abs(timeseries.df.height - med)
        mad = np.median(abs_dev)
        timeseries.df = timeseries.df.loc[abs_dev < threshold * mad]

        timeseries.df = timeseries.df.reset_index(drop=True)
        timeseries.df = timeseries.df.sort_values(by=timeseries.date_key)


def daily_mad_error(timeseries, reg_weight=0.1, reg_default=0.5, error="ADM"):
    # Group by day and get median
    medval = (
        timeseries.df.groupby(timeseries.df[timeseries.date_key])
        .median(numeric_only=True)
        .reset_index()
        .set_index(timeseries.df[timeseries.date_key])[timeseries.height_key]
    )
    day_grp = timeseries.df.groupby(timeseries.date_key).count()[timeseries.height_key]

    # Add regularization factor to avoid giving advantage to median:
    reg = reg_weight / day_grp
    reg[day_grp == 1] = reg_default
    timeseries.df[error] = (
        np.abs(timeseries.df[timeseries.height_key] - medval) + reg
    ).values
    # timeseries.df = timeseries.df.reset_index()


def daily_mean_merge(timeseries):
    # aggregates values that occur on the same day to a mean value

    dct = {
        "number": "mean",
        "object": lambda col: col.mode() if col.nunique() == 1 else np.nan,
    }

    groupby_cols = [timeseries.date_key]
    dct = {
        k: v
        for i in [
            {
                col: agg
                for col in timeseries.df.select_dtypes(tp).columns.difference(
                    groupby_cols
                )
            }
            for tp, agg in dct.items()
        ]
        for k, v in i.items()
    }
    timeseries.df = timeseries.df.groupby(groupby_cols).agg(
        **{k: (k, v) for k, v in dct.items()}
    )

    timeseries.df[timeseries.date_key] = pd.to_datetime(timeseries.df.index)
    timeseries.df = timeseries.df.reset_index(drop=True)
    timeseries.df.sort_values(by=timeseries.date_key)


def hampel(timeseries, k=7, t0=3):
    """
    vals: pandas series of values from which to remove outliers
    k: size of window (including the sample; 7 is equal to 3 on either side of value)
    """
    vals = timeseries.df.loc[:, timeseries.height_key]

    # Hampel Filter
    L = 1.4826
    rolling_median = vals.rolling(k).median()
    difference = np.abs(rolling_median - vals)
    median_abs_deviation = difference.rolling(k).median()
    threshold = t0 * L * median_abs_deviation
    outlier_idx = difference > threshold
    vals[outlier_idx] = np.nan

    timeseries.df[timeseries.height_key] = vals.values


def rolling_median(timeseries, window=7):
    timeseries.df[timeseries.height_key] = (
        timeseries.df[timeseries.height_key].rolling(window).median()
    )
