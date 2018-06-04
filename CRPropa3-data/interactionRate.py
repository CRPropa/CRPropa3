import numpy as np
from scipy.integrate import cumtrapz, romb


eV = 1.60217657e-19  # [J]
Mpc = 3.08567758e22  # [m]


def calc_rate_eps(eps, xs, gamma, field, z=0, cdf=False):
    """
    Calculate the interaction rate for given tabulated cross sections against an isotropic photon background.
    The tabulated cross sections need to be of length n = 2^i + 1 and the tabulation points log-linearly spaced.

    eps   : tabulated photon energies [J] in nucleus rest frame
    xs    : tabulated cross sections [m^2]
    gamma : (array of) nucleus Lorentz factors
    field : photon background, see photonField.py
    z     : redshift
    cdf   : calculate cumulative differential rate

    Returns :
        interaction rate 1/lambda(gamma) [1/Mpc] or
        cumulative differential rate d(1/lambda)/d(s_kin) [1/Mpc/J^2]
    """
    F = cumtrapz(x=eps, y=eps * xs, initial=0)
    n = field.getDensity(np.outer(1. / (2 * gamma), eps), z)
    if cdf:
        y = n * F / eps**2
        return cumtrapz(x=eps, y=y, initial=0) / np.expand_dims(gamma, -1) * Mpc
    else:
        y = n * F / eps
        dx = mean_log_spacing(eps)
        return romb(y, dx=dx) / gamma * Mpc


def calc_rate_s(s_kin, xs, E, field, z=0, cdf=False):
    """
    Calculate the interaction rate for given tabulated cross sections against an isotropic photon background.
    The tabulated cross sections need to be of length n = 2^i + 1 and the tabulation points log-linearly spaced.

    s_kin : tabulated (s - m**2) for cross sections [J^2]
    xs    : tabulated cross sections [m^2]
    E     : (array of) cosmic ray energies [J]
    field : photon background, see photonField.py
    z     : redshift
    cdf   : calculate cumulative differential rate

    Returns :
        interaction rate 1/lambda(gamma) [1/Mpc] or
        cumulative differential rate d(1/lambda)/d(s_kin) [1/Mpc/J^2]
    """
    F = cumtrapz(x=s_kin, y=s_kin * xs, initial=0)
    n = field.getDensity(np.outer(1. / (4 * E), s_kin), z)
    if cdf:
        y = n * F / s_kin**2
        return cumtrapz(x=s_kin, y=y, initial=0) / 2 / np.expand_dims(E, -1) * Mpc
    else:
        y = n * F / s_kin
        ds = mean_log_spacing(s_kin)
        return romb(y, dx=ds) / 2 / E * Mpc


def mean_log_spacing(x):
    """ <Delta log(x)> """
    return np.mean(np.diff(np.log(x)))


def romb_truncate(x, n):
    """ Truncate array to largest size n = 2^i + 1 """
    i = int(np.floor(np.log2(n))) + 1
    return x[0:2**i + 1]


def romb_pad_zero(x, n):
    """ Pad array with zeros """
    npad = n - len(x)
    return np.r_[x, np.zeros(npad)]


def romb_pad_logspaced(x, n):
    """ Pad array with log-linear increasing values """
    npad = n - len(x)
    dlx = np.mean(np.diff(np.log(x)))
    xpad = x[-1] * np.exp(dlx * np.arange(1, npad + 1))
    return np.r_[x, xpad]
