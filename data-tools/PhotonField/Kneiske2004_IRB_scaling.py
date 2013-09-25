from pylab import *
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.optimize import brentq
from crpropa import h_planck, c_light, k_boltzmann

# This script plots and calculates the integral of the cosmic infrared background (IRB) for different redshifts as given in Kneiske et al. 2004 (astro-ph/0309141)
# This integral can be used as an overall scaling factor for said background

def CMBSpectralRadiance(nu, z):
    """CMB spectral radiance [W/sr/m^2/Hz] at frequency nu [1/s], see http://en.wikipedia.org/wiki/Planck%27s_law"""
    return 2 * nu**3 * h_planck / c_light**2 / ( exp(nu * h_planck / k_boltzmann / 2.725 / (1 + z)) - 1 )

def CMBSpectralEnergyDensity(nu, z):
    """CMB spectral energy density [J/m^3/Hz] at frequency nu [1/s]"""
    return CMBSpectralRadiance(nu, z) * 4 * pi / c_light

def CMBSpectralNumberDensity(nu, z):
    """CMB spectral number density [1/m^3/Hz] at frequency nu [1/s]"""
    return CMBSpectralEnergyDensity(nu, z) / (nu * h_planck)

def process(fname, z):
    # read file
    # wavelength lambda [mu m] and wavelength scaled spectral radiance lambda * I_lambda [nW/m^2/sr]
    x, y = genfromtxt(fname, unpack=True)
    y[y < 1.01e-1] = 0 # cut off values that are artefacts from parsing the plots from the paper

    x *= 1e-6 # [mu m] -> [m]
    y *= 1e-9 # [nW/m^2/sr] -> [W/m^2/sr]
    
    # remove wavelength scaling 
    y /= x # [W/m^2/sr] -> [W/m^3/sr]
    
    # convert to spectral energy density dEdV_lambda
    y *= 4 * pi / c_light # [W/m^3/sr] -> [J/m^3/m]

    # convert to spectral number density: divide by photon energy
    y /= h_planck * c_light / x # [J/m^3/m] -> [1/m^3/m]

    # change to frequency as spectral variable: B_nu = lambda^2 / h * B_lambda
    nu = c_light / x # [m] -> [Hz]
    ynu = y * x**2 / c_light # [1/m^3/Hz] 

    # find frequency nu0 above which the IRB number density exceeds the CMB
    nu = nu[::-1]
    ynu = ynu[::-1]
    ynu_interp = interp1d(nu, ynu)
    def objective(nu):
        return CMBSpectralNumberDensity(nu, z) - ynu_interp(nu)
    nu0 = brentq(objective, 3.5e11, 1e13) # frequency of intersection

    # integral of spectral number density from nu0 to nuMax
    I1 = quad(ynu_interp, nu[0], nu[-1])[0]
    I2 = quad(ynu_interp, nu0, nu[-1])[0]

    # plot
    figure()
    nuCMB = logspace(9, 14)
    plot(nuCMB, CMBSpectralNumberDensity(nuCMB, z), label='CMB')
    plot(nu, ynu, label='IRB Kneiske')
    xlabel('Frequency [Hz]')
    ylabel('Number Density [1/m$^{3}$/Hz]')
    ylim(1e-16, 1)
    xlim(1e9, 1e16)
    text(5e14, 0.1, 'z = %.1f'%z)
    loglog()
    axvline(nu0, lw=1, ls='--', color='k')
    savefig('z%.1f.png'%z, bbox_inches='tight')

    print z, nu0, I1, I2
    return I1, I2


fnames = ['z0.txt', 'z0.2.txt', 'z0.4.txt', 'z0.6.txt', 'z1.txt', 'z2.txt', 'z3.txt', 'z4.txt']
z = array([0, 0.2, 0.4, 0.6, 1, 2, 3, 4, 5])

s1, s2 = zeros(9), zeros(9)
for i in range(8):
    s1[i], s2[i] = process('Kneiske2004_IRB/'+fnames[i], z[i])

s1 /= s1[0]
s2 /= s2[0]

figure()
plot(z, s1)
plot(z, s2)
xlim(0, 5)
xlabel('Redshift z')
ylabel('Overall Scaling')
savefig('Kneiske2004_IRB_scaling.png',bbox_inches='tight')
show()
