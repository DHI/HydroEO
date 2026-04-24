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
