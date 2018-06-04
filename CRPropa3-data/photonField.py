import numpy as np
from os import path

cdir = path.split(__file__)[0]
datadir = path.join(cdir, 'tables/')

eV = 1.60217657e-19  # [J]
erg = 1e-7  # [J]
c0 = 299792458  # [m/s]
h = 6.62606957e-34  # [m^2 kg / s]
kB = 1.3806488e-23  # [m^2 kg / s^2 / K]
T_CMB = 2.72548  # CMB temperature [K]


# --------------------------------------------------------
# interfaces
# --------------------------------------------------------
class CMB:
    """
    Cosmic microwave background radiation
    """
    name = 'CMB'
    info = 'CMB'
    redshift = None
    def getDensity(self, eps, z=0):
        """
        Comoving spectral number density dn/deps [1/m^3/J] at given photon energy eps [J] and redshift z.
        Multiply with (1+z)^3 for the physical number density.
        """
        return 8*np.pi / c0**3 / h**3 * eps**2 / (np.exp(eps/(kB*T_CMB)) - 1)

    def getEmin(self, z=0):
        """Minimum effective photon energy in [J]"""
        return 1e-10 * eV

    def getEmax(self, z=0):
        """Maximum effective photon energy in [J]"""
        return 0.1 * eV

class EBL:
    """
    Base class for extragalactic background light (EBL) models
    """
    def __init__(self):
        self.data = {}  # dictionary {redshift : (eps, dn/deps)}

    def getDensity(self, eps, z=0):
        """
        Comoving spectral number density dn/deps [1/m^3/J] at given photon energy eps [J] and redshift z.
        Multiply with (1+z)^3 for the physical spectral number density.
        The tabulated data is interpolated linearly in log-log; no extrapolation is performed.
        """
        tab_eps, tab_n = self.data[z]
        # return np.interp(eps, tab_eps, tab_n, 0, 0)  # linear interpolation
        x  = np.log10(eps)
        tx = np.log10(tab_eps)
        ty = np.log10(tab_n)
        return 10**np.interp(x, tx, ty, -np.inf, -np.inf)  # log-log interpolation

    def getEmin(self, z=0):
        """Minimum tabulated photon energy in [J]"""
        return self.data[z][0][0]

    def getEmax(self, z=0):
        """Maximum tabulated photon energy in [J]"""
        return self.data[z][0][-1]

# --------------------------------------------------------
# EBL (optical and infrared) models
# --------------------------------------------------------
class EBL_Kneiske04(EBL):
    name = 'IRB_Kneiske04'
    info = 'cosmic infrared and optical background radiation model of Kneiske et al. 2004'
    files = datadir + 'EBL_Kneiske_2004/all_z'
    redshift = np.linspace(0, 5, 51)

    def __init__(self):
        EBL.__init__(self)
        # d[0] : eps [eV]
        # d[1-51] : n(eps), [1/m^3/eV]
        d = np.genfromtxt(self.files, unpack=True)
        eps = d[0] * eV
        n = d[1:] / eV
        for i,z in enumerate(self.redshift):
            self.data[z] = eps, n[i]

class EBL_Kneiske10(EBL):
    name = 'IRB_Kneiske10'
    info = 'cosmic infrared and optical background radiation lower limit model of Kneiske et al. 2010'
    files = datadir + 'EBL_Kneiske_2010/%.1f'
    redshift = (0, .1, .3, .8, 2)

    def __init__(self):
        EBL.__init__(self)
        for z in self.redshift:
            # x : wavelength in [mu m]
            # y : lambda I_lambda [nW/m^2/sr]
            x, y = np.genfromtxt(self.files%z, unpack=True)
            wl = (10**x * 1e-6)  # wavelength in [m]
            eps = h * c0 / wl
            n = (10**y * 1e-9) * (4 * np.pi / c0) / eps**2
            self.data[z] = eps[::-1], n[::-1]

class EBL_Dole06(EBL):
    name = 'IRB_Dole06'
    info = 'cosmic infrared and optical background radiation model of Dole et al. 2006'
    files = datadir + 'EBL_Dole_2006/0.0'
    redshift = (0)

    def __init__(self):
        EBL.__init__(self)
        # d[0] : lambda [mu m]
        # d[1] : n(eps), [W/m^2/sr]
        d = np.genfromtxt(self.files, unpack=True)
        eps = h * c0 / (d[0] * 1e-6)  # photon energy [J]
        n = d[1] * (4 * np.pi / c0) / eps**2
        self.data[0] = eps[::-1], n[::-1]

class EBL_Franceschini08(EBL):
    name = 'IRB_Franceschini08'
    info = 'cosmic infrared and optical background radiation model of Franceschini et al. 2008'
    files = datadir + 'EBL_Franceschini_2008/%1.1f'
    redshift = (0, .2, .4, .6, .8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0)

    def __init__(self):
        EBL.__init__(self)
        for z in self.redshift:
            # x : log10(eps / eV)
            # y : eps * dn/deps [1/cm^3]
            x, y = np.genfromtxt(self.files%z, unpack=True)
            eps = 10**x * eV
            n = 10**y / eps * 1e6
            n /= (1 + z)**3  # make comoving
            self.data[z] = eps, n

class EBL_Stecker05(EBL):
    name = 'IRB_Stecker05'
    info = 'cosmic infrared and optical background radiation model of Stecker at al. 2005'
    files = datadir + 'EBL_Stecker_2005/data2.txt'
    redshift = np.linspace(0, 5, 26)

    def __init__(self):
        EBL.__init__(self)
        # d[0]    : log10(eps/eV)
        # d[1-26] : log10(eps*n(eps)/cm^3), not comoving
        d = np.genfromtxt(self.files, unpack=True)
        eps = 10**d[0] * eV
        n = 10**d[1:] / eps * 1e6
        for i,z in enumerate(self.redshift):
            # convert n(eps) to comoving density
            self.data[z] = eps, n[i] / (1+z)**3

class EBL_Finke10(EBL):
    name = 'IRB_Finke10'
    info = 'cosmic infrared and optical background radiation model of Finke et al. 2010 (Model C)'
    files = datadir + 'EBL_Finke_2010/z%.2f.dat'
    redshift = np.arange(0, 5, 0.01)

    def __init__(self):
        EBL.__init__(self)
        for z in self.redshift:
            # d[0] : eps / eV
            # d[1] : comoving energy density in erg / cm^3
            d = np.genfromtxt(self.files % z, unpack=True)
            eps = d[0] * eV
            n   = d[1] * erg * 1e6 / eps**2 # [J/m^3]
            self.data[z] = eps, n

class EBL_Gilmore12(EBL):
    name = 'IRB_Gilmore12'
    info = 'cosmic infrared and optical background radiation model of Gilmore et al. 2012 (Evolving dust model, arXiv:1104.0671)'
    files = datadir + 'EBL_Gilmore_2012/eblflux_fiducial.dat'
    redshift = np.array([0, 0.015, 0.025, 0.044, 0.05, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0])

    def __init__(self):
        EBL.__init__(self)
        # d[0] : rest frame wavelength in [angstrom]
        # d[1-21] : proper flux [erg/s/cm^2/ang/sr]
        d = np.genfromtxt(self.files, unpack=True)
        eps = h * c0 / (d[0] * 1e-10)  # [J]
        n = d[1:] * erg / 1e-4 * d[0] / eps**2 * (4 * np.pi / c0) 
        for i,z in enumerate(self.redshift):
            n[i] /= (1 + z)**3  # make comoving
            self.data[z] = eps[::-1], n[i][::-1]

class EBL_Dominguez11(EBL):
    def __init__(self, which='best'):
        EBL.__init__(self)
    
        if which == 'best':
            fname = 'EBL_Dominguez_2011/ebl_dominguez11.out'
            self.name = 'IRB_Dominguez11'
        elif which == 'upper':
            fname = 'EBL_Dominguez_2011/ebl_upper_uncertainties_dominguez11.out'
            self.name = 'IRB_Dominguez11_upper'
        elif which == 'lower':
            fname = 'EBL_Dominguez_2011/ebl_lower_uncertainties_dominguez11.out'
            self.name = 'IRB_Dominguez11_lower'
        else:
            raise ValueError('EBL_Dominguez11 only provides "best", "upper" and "lower" models')

        self.info = 'cosmic infrared and optical background radiation (%s) model of Dominguez et al. 2011 (arXiv:1007.1459)' % which
        self.redshift = np.array([0, 0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.9])

        # d[0] : rest frame wavelength in [mu m]
        # d[1-18] : proper flux [nW/m^2/sr]
        d = np.genfromtxt(datadir + fname, unpack=True)
        eps = h * c0 / (d[0] * 1e-6)  # [J]
        n = d[1:] * 1e-9 / eps**2 * (4 * np.pi / c0)
        for i, z in enumerate(self.redshift):
            self.data[z] = eps[::-1], n[i][::-1]  # sort by ascending energy

class EBL_Stecker16(EBL):
    def __init__(self, which='upper'):
        EBL.__init__(self)

        if which == 'upper':
            fname = 'EBL_Stecker_2016/comoving_enerdens_up.csv'
        elif which == 'lower':
            fname = 'EBL_Stecker_2016/comoving_enerdens_lo.csv'
        else:
            raise ValueError('EBL_Stecker16 only provides "upper" and "lower" models')

        self.name = 'IRB_Stecker16_%s' % which
        self.info = 'cosmic infrared and optical background radiation (%s) model of Stecker et al. 2016 (arXiv:1605.01382)' % which
        self.redshift = np.linspace(0, 5, 501)
        
        # rows:    log10(photon energy / eV) ranging from -2.84 to 1.14 in steps of 0.01
        # columns: redshifts from 0 to 5 in steps of 0.01
        # units:   erg / Hz / cm^3 --> 
        d = np.genfromtxt(datadir + fname, delimiter=',').T
        eps = 10**np.arange(-2.84, 1.14001, 0.01) * eV
        n = d * erg / h * 1E6 / eps
        for i, z in enumerate(self.redshift):
            self.data[z] = eps, n[i]

# --------------------------------------------------------
# CRB (radio) models
# --------------------------------------------------------
class URB_Protheroe96:
    """
    Universal Radio Background from Protheroe, Bierman 1996.
    Taken from EleCa implementation.
    """
    name = "URB_Protheroe96"
    info = "URB_Protheroe96"

    def getDensity(self, eps, z=0):
        """
        Comoving spectral number density dn/deps [1/m^3/J] at given photon energy eps [J]
        """
        p0 = -2.23791e+01
        p1 = -2.59696e-01
        p2 = 3.51067e-01
        p3 = -6.80104e-02
        p4 = 5.82003e-01
        p5 = -2.00075e+00
        p6 = -1.35259e+00
        p7 = -7.12112e-01  # xbreak

        eps = np.r_[eps]
        x = np.log10(eps / h / 1e9)
        I = p0 + p1 * x + p2 * x**2 + p3 * x**3 / (np.exp(p4 * x) - 1)
        I[x > p7] += p6 + p5 * x[x > p7] - p2 * x[x > p7]**2
        I = 4 * np.pi / (h * c0) * (10**I / eps)

        I[eps < self.getEmin()] = 0
        I[eps > self.getEmax()] = 0
        return I

    def getEmin(self, z=0):
        """Minimum effective photon energy in [J]"""
        return 4.1e-12 * eV

    def getEmax(self, z=0):
        """Maximum effective photon energy in [J]"""
        return 2E-6 * eV # 0.825e-6 * eV


class CRB_ARCADE2:
    def getDensity(self, eps):
        """
        Spectral number density dn/deps [1/m^3/J] at z = 0.
        """
        # T = 1.26 +- 0.09 K (nu/GHz)^-2.6 +- 0.04, see Holder 2012
        T = 1.26 * (nu/1e9)**-2.6
        return 8*np.pi / c0**3 / h**3 * eps**2 / (np.exp(eps/(kB*T)) - 1)

    def getEmin(self, z=0):
        """Minimum effective photon energy in [J]"""
        return 1e-10 * eV

    def getEmax(self, z=0):
        """Maximum effective photon energy in [J]"""
        return 0.1 * eV


if __name__ == '__main__':
    from pylab import *
    eps = logspace(-3, 1, 200) * eV
    x  = eps / eV
    c =  eps**2 / eV
    y1   = c * EBL_Kneiske04().getDensity(eps)
    y3   = c * EBL_Stecker05().getDensity(eps)
    y5   = c * EBL_Franceschini08().getDensity(eps)
    y6   = c * EBL_Finke10().getDensity(eps)
    y7   = c * EBL_Dominguez11().getDensity(eps)
    y8   = c * EBL_Gilmore12().getDensity(eps)
    y7up = c * EBL_Dominguez11('upper').getDensity(eps)
    y7lo = c * EBL_Dominguez11('lower').getDensity(eps)
    y9up = c * EBL_Stecker16('upper').getDensity(eps)
    y9lo = c * EBL_Stecker16('lower').getDensity(eps)

    figure()
    plot(x, y1, label='Kneiske 2004')
    plot(x, y3, label='Stecker 2005')
    plot(x, y5, label='Franceschini 2008')
    plot(x, y6, label='Finke 2010')
    plot(x, y7, label='Dominguez 2011')
    plot(x, y8, label='Gilmore 2012')
    fill_between(x, y7lo, y7up, facecolor='m', edgecolor='none', alpha=0.2, zorder=-1, label='Dominguez 2011 (limits)')
    fill_between(x, y9lo, y9up, facecolor='g', edgecolor='none', alpha=0.2, zorder=-1, label='Stecker 2016 (limits)')

    legend(loc='lower center', fontsize='x-small')
    loglog()
    grid()
    # ylim(1e1, 2e6)
    ylabel('$\epsilon^2 ~ dn/d\epsilon$ [eV/m$^3$]')
    xlabel('$\epsilon$ [eV]')
    savefig('figures/EBL.png')
    show()
