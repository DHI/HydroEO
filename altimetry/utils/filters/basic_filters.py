"""simple filters that can be applied to sat timeseries objects"""

import numpy as np
import pandas as pd
from sklearn.svm import SVR

# from tqdm import tqdm

from datetime import datetime


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

    return timeseries


def daily_mad_error(timeseries, reg_weight=0.1, reg_default=0.5, error_key="ADM"):
    # sort inplace the timeseries object
    timeseries.df = timeseries.df.sort_values(by=timeseries.date_key).reset_index(
        drop=True
    )

    # exctract for calculations
    df = timeseries.df.copy()
    date_key = timeseries.date_key
    height_key = timeseries.height_key

    # assign a group to each unique day
    unique_days = df[date_key].dt.date.unique()

    # Group by day and get median
    medval = df.groupby(df[date_key]).median(numeric_only=True)[height_key]
    day_grp = df.groupby(date_key).count()[height_key]

    # Add regularization factor to avoid giving advantage to median:
    reg = reg_weight / day_grp
    reg[day_grp == 1] = reg_default

    # map the median val and regularization to the dates they belong to
    med_map = dict(zip(unique_days, medval))
    reg_map = dict(zip(unique_days, reg))
    df["med"] = df[date_key].map(med_map)
    df["reg"] = df[date_key].map(reg_map)

    # calculate error
    error = (np.abs(df[height_key] - df["med"]) + df["reg"]).values

    # add error to timeseries
    timeseries.df[error_key] = error

    return timeseries


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

    return timeseries


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

    return timeseries


def rolling_median(timeseries, window=7):
    timeseries.df[timeseries.height_key] = (
        timeseries.df[timeseries.height_key].rolling(window).median()
    )

    return timeseries


def _run_svr_linear(heights, err=0.1, epsilon=0.1):
    """
    Linear Support Vector Regression
    Fit 0-slope line through heights along-track to remove outliers.
    This method allows non-0 slope

    Parameters
    ----------
    heights : array
        Heights to be fit.
    err : Float, optional
        Allowed deviation from linear regression. The default is .01.
    epsilon : Float, optional
        "Epsilon in the epsilon-SVR model.
        It specifies the epsilon-tube within which no penalty is associated in the training loss
        function with points predicted within a distance epsilon from the actual value."
        (from https://scikit-learn.org/stable/modules/generated/sklearn.svm.SVR.html)
        The default is .1.

    Returns
    -------
    array
        Outlier filtered heights.

    """
    # sequential x axis along track
    x = np.arange(0, len(heights))
    x = np.vstack(
        [x, np.ones(len(x))]
    ).T  # Extend x data to contain another row vector of 1s
    y = heights.values  # TODO: investigate First remove slope if any

    # make SVR kernel and fit with confidence bounds
    svr_rbf = SVR(kernel="linear", epsilon=epsilon)
    rbf = svr_rbf.fit(x, y)
    uconf = rbf.predict(x) + err
    lconf = rbf.predict(x) - err

    # only keep the filtered valeus
    filtered = np.where((y >= lconf) & (y <= uconf))[0]

    return np.array(heights)[filtered]


def svr_linear(
    timeseries,
):  # TODO: this should maybe be processed on the individual product timeseries?
    df = timeseries.df
    date_key = timeseries.date_key
    height_key = timeseries.height_key

    # get the remaining heights after the filter
    lin_filt = df.groupby(date_key)[height_key].apply(_run_svr_linear).reset_index()

    # reassign the values to their date (essentially, "ungrouby")
    r = pd.DataFrame(
        {
            col: np.repeat(lin_filt[col].values, lin_filt[height_key].str.len())
            for col in lin_filt.columns.drop(height_key)
        }
    )
    r = r.assign(**{height_key: np.concatenate(lin_filt[height_key].values)})[
        lin_filt.columns
    ]

    # Reset the time series to the filtered values
    timeseries.df = r

    return timeseries


"""Kalman filter for estimating state of reservoir from noisy timeseries"""


def _update(obs, xk, cov_xx, height="height", error="ADM", n=1):
    """
    Update function for Kalman filter

    Parameters
    ----------
    obs : array
        DESCRIPTION.
    xk : float
        Previous prediction, x.
    cov_xx : float
        Uncertainty of x.
    height : string, optional
        DESCRIPTION. The default is 'height'.
    error : string, optional
        DESCRIPTION. The default is 'ADM'.
    n : Int, optional
        Grid size - only relevant for large lakes e.g. Not yet implemented. The default is 1.

    Returns
    -------
    xk_plus : Float
        Kalman filter updated value.
    cov_xx_plus : Float
        Covariance of the updated value.

    """
    obs_list = obs[height].values
    m = len(obs_list)
    lk = np.array(obs_list).reshape(m, n)
    Ak = np.ones((m, n))

    # compute Kalman matrix (weights of the innovation)
    slk = (obs[error].values) ** 2
    cov_lk = np.diag(slk)
    inv = Ak * cov_xx * np.transpose(Ak) + cov_lk
    Kk = cov_xx * np.dot(np.transpose(Ak), np.linalg.pinv(inv))

    # update x
    xk_plus = xk + np.dot(Kk, lk - Ak * xk)

    # update sigma
    cov_xx_plus = cov_xx * (1 - np.dot(Kk, Ak))

    return xk_plus, cov_xx_plus


def _pred(xk, cov_xx, n=1, system_noise=0, system_noise_unc=0.05):
    """
    Prediction function for Kalman filter

    Parameters
    ----------
    xk : float
        Updated prediction, x.
    cov_xx : float
        Uncertainty of x.
    n : Int, optional
        Grid size - only relevant for large lakes e.g. Not yet implemented. The default is 1.
    system_noise : float, optional
        System noise. The default is 0.
    system_noise_unc : float, optional
        Uncertainty of the system noise. The default is 0.05.

    Returns
    -------
    xk_plus : Float
        Kalman filter predicted value.
    cov_xx_plus : Float
        Covariance of the predicted value.

    """
    # Setup dynamic model with no deterministic contribution
    thetak = np.identity(n)
    hatk = np.identity(n)
    Qk = np.array(system_noise_unc).reshape(n, n)
    qk = np.array(system_noise).reshape(n, 1)
    noise = hatk.dot(Qk).dot(hatk.T)

    # pred x
    xk_next = thetak * xk + hatk * qk
    # predict covariance matrix
    cov_xx_next = thetak.dot(cov_xx).dot(thetak.T) + noise

    return xk_next, cov_xx_next


def kalman(timeseries, error_key="ADM", n=1):
    """
    Run Kalman filter

    Parameters
    ----------
    time_series : DataFrame
        Outlier filtered dataframe to be used as input for Kalman filter.
        Must contain height and error columns
    height : string, optional
        Name of height column. The default is 'height_OCOG'.
    error : string, optional
        Name of error column. The default is 'ADM'.
    n : Int, optional
        Grid size - only relevant for large lakes e.g. Not yet implemented. The default is 1.

    Returns
    -------
    xk_plus : array
        Kalman filter updated values of WSE.
    cov_xx_plus : array
        Covariance of the updated value.

    """

    df = timeseries.df.copy()
    date_key = timeseries.date_key
    height_key = timeseries.height_key

    dates = sorted(df[date_key].unique())
    xks = np.ones((n, len(df[date_key].unique()))) * np.nan
    cov_xxs = np.ones((n, n, len(df[date_key].unique()))) * np.nan

    # Initialize prediction and uncertainty matrices
    obs = df.loc[df[date_key] == dates[0]]
    lk = obs[height_key].values
    slk = (obs[error_key].values) ** 2

    # First prediction = value with smallest ADM and identity matrix of size n, n
    xks[:, 0] = lk[np.argmin(slk)]
    cov_xxs[:, :, 0] = np.identity(n)

    # Observation model
    for i, d in enumerate(dates):
        # Update
        obs = df.loc[df[date_key] == d]
        xk_plus, cov_xx_plus = _update(
            obs, xks[:, i], cov_xxs[:, :, i], height=height_key, error=error_key, n=n
        )

        # Predict
        xk1, cov_xx1 = _pred(
            xk_plus, cov_xx_plus, n=n, system_noise=0, system_noise_unc=0.05
        )
        xks[:, i] = xk_plus
        cov_xxs[:, :, i] = cov_xx_plus
        try:
            xks[:, i + 1] = xk1
            cov_xxs[:, :, i + 1] = cov_xx1
        except IndexError:
            pass

    # return a new dataframe with the filtered data
    df_kalman = pd.DataFrame(
        {date_key: dates, height_key: xks[0], "cov": cov_xxs[0][0]}
    )
    return df_kalman


def _run_svr_rbf(dates, heights, err=1, rbf_c=10000, gamma=0.0000438, epsilon=0.1):
    """
    Radial Base Function outlier filtering post-Kalman filter.
    This is run at virtual station level

    Parameters
    ----------
    dates : array
        Sorted dates of altimetry observations.
    heights : array
        Predicted (and updated) water surface elevation.
    err : Float, optional
        Observation uncertainty. The default is 1 m. As used in DAHITI for rivers
        (0.1 m can be used for lakes)
    C : TYPE, optional
        Regularization parameter in sklearn.svm.RBF. The default here is 10000 (as in DAHITI)
    gamma : TYPE, optional
        Kernel coefficient for ‘rbf’ in sklearn.svm.RBF
        Alternative values include float, 'scale' or 'auto'
        The default here is 0.0000438 (as in DAHITI)
    epsilon : TYPE, optional
        "Epsilon in the epsilon-SVR model.
        It specifies the epsilon-tube within which no penalty is associated in the training loss
        function with points predicted within a distance epsilon from the actual value."
        (from https://scikit-learn.org/stable/modules/generated/sklearn.svm.SVR.html)
        The default is .1.

    Returns
    -------
    rbf_filter : array
        Filter for Kalman filter predictions of WSE
    """
    # Day since start
    x = np.array(
        [0]
        + list(
            np.cumsum(
                pd.to_datetime(dates)[1:] - pd.to_datetime(dates)[:-1]
            ).days.astype("float32")
        )
    )
    # Kalman heights
    y = heights
    x = np.vstack([x, np.ones(len(x))]).T
    # y = np.vstack([y, np.ones(len(y))]).T

    # Scale values:
    # std_x = StandardScaler()
    # std_y = StandardScaler()
    # x2 = std_x.fit_transform(x)
    # y2 = std_y.fit_transform(y)[:,0]

    # Extend x data to contain another row vector of 1s

    svr_rbf = SVR(kernel="rbf", C=rbf_c, gamma=gamma, epsilon=epsilon)
    rbf = svr_rbf.fit(x, y)
    corr = rbf.predict(x)

    uconf = corr + err
    lconf = corr - err

    rbf_filter = np.where((y >= lconf) & (y <= uconf))

    return rbf_filter


def _year_fraction(dt):
    start = datetime(dt.year, 1, 1).toordinal()
    year_length = datetime(dt.year + 1, 1, 1).toordinal() - start
    return dt.year + float(dt.toordinal() - start) / year_length


def svr_radial(timeseries):
    df = timeseries.df.copy()
    date_key = timeseries.date_key
    height_key = timeseries.height_key
    # error_key = timeseries.error_key

    rbf_filter = _run_svr_rbf(
        df[date_key].values,
        df[height_key].values,
        err=1,  # TODO: consider 0.1 for lakes after bias correction?
        rbf_c=1000,
        gamma=0.0000438,
        epsilon=0.1,
    )  # these are the default values as in DAHITI with error changed from 1 to 0.1m for lakes

    nb_obs = df.groupby(date_key).count().reset_index()
    # pass_ = df.groupby(date_key).first().reset_index()
    # coords = df.loc[df.groupby(date_key)[error_key].idxmin()].reset_index()

    nb_obs["nb_obs"] = nb_obs[date_key]
    # nb_obs[["relative orbit", "orbit", "platform"]] = pass_[
    #     ["relative_orbit", "orbit", "platform"]
    # ]
    # nb_obs[["lat", "lon"]] = coords[["lat", "lon"]]

    # Create new dataframe:
    vs = pd.DataFrame(
        {
            date_key: df[date_key].values[rbf_filter],
            height_key: df[height_key].values[rbf_filter],
        }
    )
    # vs = pd.merge(
    #     nb_obs[["day", "nb_obs", "relative orbit", "orbit", "platform", "lat", "lon"]],
    #     vs,
    #     left_on="day",
    #     right_on="date",
    # )

    # vs["decimal year"] = [
    #     str(round(_year_fraction(pd.to_datetime(uniq_d)), 8)) for uniq_d in vs.date
    # ]
    # vs["date"] = [
    #     datetime.strftime(pd.to_datetime(uniq_d), "%Y/%m/%d") for uniq_d in vs.date
    # ]

    return vs
