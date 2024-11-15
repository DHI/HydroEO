"""simple filters that can be applied to sat timeseries objects"""

import numpy as np


def elevation_filter(timeseries, height_range):
    timeseries.df = timeseries.df.loc[timeseries.df[timeseries.height_key] > 0]
    timeseries.df = timeseries.df.loc[timeseries.df[timeseries.height_key] < 8000]


def daily_mean_filter(timeseries):
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
    agg = timeseries.df.groupby(groupby_cols).agg(**{k: (k, v) for k, v in dct.items()})

    timeseries.df = agg
