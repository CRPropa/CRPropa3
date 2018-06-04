"""
Calculate the energy loss rate through electron pair production
References:
(B70) Blumenthal 1970, Phys.Rev. D
(C92) Chodorowski et al. 1992, ApJ 400:181-185
"""
import numpy as np
from scipy import integrate
import photonField


eV = 1.60217657e-19  # [J]
Mpc = 3.08567758e22  # [m]
c0 = 299792458.  # [m/s]
r0 = 2.817940e-15  # classical electron radius [m]
alpha = 7.297352e-3  # fine-structure constant
me = 9.10938291e-31  # electron mass [kg]
me_c2 = me * c0**2  # electron mass in [J/c^2]
mp = 1.67262178e-27  # proton mass [kg]


def lossRate(gamma, field, z=0):
    """
    Loss rate from electron pair production with the given photon background, cf. C92, equation 3.11
    gamma   : list of nucleus Lorentz factors
    field   : photon background
    z       : redshift
    Returns : 1/gamma dgamma/dx [1/Mpc]
    """

    def phi(k):
        """
        Parametrization of the integral 3.12 (C92)
        """
        _c = np.array([0.8048, 0.1459, 1.137e-3, -3.879e-6])
        _d = np.array([-86.07, 50.96, -14.45, 8 / 3.])
        _f = np.array([2.910, 78.35, 1837])
        # phi(k) for k < 25, eq. 3.14
        if k < 25:
            return np.pi / 12 * (k - 2)**4 / (1 + sum(_c * (k - 2)**np.arange(1, 5)))
        # phi(k) for k > 25, eq. 3.18 and 3.16
        return k * sum(_d * np.log(k)**np.arange(4)) / (1 - sum(_f * k**-np.arange(1, 4)))

    def integrand(logk, gamma, field):
        """
        Integrand of equation 3.11 (C92), logarithmic version
        logk  : ln(k) = ln(2 gamma eps / (me c^2)), photon energy
        gamma : nucleus Lorentz factor
        field : photon background
        """
        k = np.exp(logk)
        eps = k * me_c2 / 2 / gamma  # photon energy [J] in lab frame
        n = field.getDensity(eps, z)  # spectral number density [1/m^3/J]
        n *= me_c2  # from substitution d eps / d k
        return n * phi(k) / k  # includes *k from substitution k -> ln(k)

    rate = np.zeros_like(gamma)
    err = np.zeros_like(gamma)

    # minimum and maximum energy of the fields photons in units of me*c^2
    epsmin = field.getEmin() / me_c2
    epsmax = field.getEmax() / me_c2

    for i, g in enumerate(gamma):
        lkmin = np.log(max(2, 2 * g * epsmin))
        lkmax = np.log(2 * g * epsmax)
        lksep = np.linspace(lkmin, lkmax, 11)[1:-1]
        rate[i], err[i] = integrate.quad(
            integrand, lkmin, lkmax, points=lksep, args=(g, field))

    # prefactor of equation 3.11 (C92) and conversion [1/s] --> [1/Mpc]
    a = alpha * r0**2 * me / mp * Mpc
    return a * rate / gamma, a * err / gamma


# -------------------------------------------------
# Generate tables for energy loss rate
# -------------------------------------------------
gamma = np.logspace(6, 14, 161)  # tabulated Lorentz factors

fields = [
    photonField.CMB(),
    photonField.EBL_Kneiske04(),
    photonField.EBL_Stecker05(),
    photonField.EBL_Franceschini08(),
    photonField.EBL_Finke10(),
    photonField.EBL_Dominguez11(),
    photonField.EBL_Gilmore12(),
    photonField.EBL_Stecker16('upper'),
    photonField.EBL_Stecker16('lower')]

for field in fields:
    print(field.name)
    rate = lossRate(gamma, field)[0]
    s = (rate > 1e-12)  # truncate if loss rate is < 10^-12 / Mpc

    fname = 'data/ElectronPairProduction/lossrate_%s.txt' % field.name
    data = np.c_[np.log10(gamma[s]), rate[s]]
    fmt = '%.2f\t%.6e'
    header = 'Loss rate for electron-pair production with the %s\nlog10(gamma)\t1/gamma dgamma/dx [1/Mpc]' % field.info
    np.savetxt(fname, data, fmt=fmt, header=header)


# -------------------------------------------------
# Reformat CRPropa2 tables of differential spectrum of secondary electrons
# This should be reimplemented for extension to the other backgrounds,
# cross-checking and documentation.
# -------------------------------------------------
d1 = np.genfromtxt('tables/EPP/pair_spectrum_cmb.table', unpack=True)
d2 = np.genfromtxt('tables/EPP/pair_spectrum_cmbir.table', unpack=True)

# amplitudes dN/dEe(Ep)
A1 = d1[2].reshape((70, 170))  # CMB
A2 = d2[2].reshape((70, 170))  # CMB + IRB (which?)
A3 = A2 - A1  # IRB only

# # normalize to 1
# A1 = (A1.T / sum(A1, axis=1)).T
# A3 = (A3.T / sum(A3, axis=1)).T

# save
np.savetxt('data/ElectronPairProduction/spectrum_CMB.txt', A1, fmt='%.5e')
np.savetxt('data/ElectronPairProduction/spectrum_IRB.txt', A3, fmt='%.5e')
