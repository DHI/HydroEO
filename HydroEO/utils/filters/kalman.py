"""Kalman filter for estimating state of reservoir from noisy timeseries"""

import numpy as np


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


# def kalman(time_series, date_key="day", height_key="height", error_key="ADM", n=1):
#     """
#     Run Kalman filter

#     Parameters
#     ----------
#     time_series : DataFrame
#         Outlier filtered dataframe to be used as input for Kalman filter.
#         Must contain height and error columns
#     height : string, optional
#         Name of height column. The default is 'height_OCOG'.
#     error : string, optional
#         Name of error column. The default is 'ADM'.
#     n : Int, optional
#         Grid size - only relevant for large lakes e.g. Not yet implemented. The default is 1.

#     Returns
#     -------
#     xk_plus : array
#         Kalman filter updated values of WSE.
#     cov_xx_plus : array
#         Covariance of the updated value.

#     """

#     dates = sorted(time_series[date_key]unique())
#     xks = np.ones((n, len(time_series[date_key].unique()))) * np.nan
#     cov_xxs = np.ones((n, n, len(time_series[date_key].unique()))) * np.nan

#     # Initialize prediction and uncertainty matrices
#     obs = time_series.loc[time_series[date_key] == dates[0]]
#     lk = obs[height_key].values
#     slk = (obs[error_key].values) ** 2

#     # First prediction = value with smallest ADM and identity matrix of size n, n
#     xks[:, 0] = lk[np.argmin(slk)]
#     cov_xxs[:, :, 0] = np.identity(n)

#     # Observation model
#     for i, d in enumerate(dates):
#         # Update
#         obs = time_series.loc[time_series[date_key] == d]
#         xk_plus, cov_xx_plus = _update(
#             obs, xks[:, i], cov_xxs[:, :, i], height=height_key, error=error_key, n=n
#         )

#         # Predict
#         xk1, cov_xx1 = _pred(
#             xk_plus, cov_xx_plus, n=n, system_noise=0, system_noise_unc=0.05
#         )
#         xks[:, i] = xk_plus
#         cov_xxs[:, :, i] = cov_xx_plus
#         try:
#             xks[:, i + 1] = xk1
#             cov_xxs[:, :, i + 1] = cov_xx1
#         except IndexError:
#             pass

#     return xks, cov_xxs
