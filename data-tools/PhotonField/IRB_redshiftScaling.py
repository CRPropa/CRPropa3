from pylab import *
from scipy.interpolate import interp1d
from scipy.integrate import quad
from mpc import h_planck, c_light

# This script plots and calculates the integral of the cosmic infrared background (IRB) for different redshifts as given in Kneiske et al. 2004 (astro-ph/0309141)
# This integral can be used as an overall scaling factor for said background


def process1(fin):
    x,y = genfromtxt(fin, unpack=True)
    y[y < 1.01e-1] = 0
    return sum( (y[1:] + y[:-1]) / 2 * (x[1:] - x[:-1]) )


def numberDensity(fin='Kneiske2004_IRB/z0.txt'):
    """IRB (Kneiske et al. 2004) number density"""
    # read file
    # wavelength lambda [mu m] and wavelength scaled spectral radiance lambda * I_lambda [nW/m^2/sr]
    x, y = genfromtxt(fin, unpack=True)
    y[y < 1.01e-1] = 0 # cut off values that are artefacts from parsing the plots from the paper

    x *= 1e-6 # [mu m] -> [m]
    y *= 1e-9 # [nW/m^2/sr] -> [W/m^2/sr]
    
    # remove wavelength scaling 
    y /= x # [W/m^2/sr] -> [W/m^3/sr]
    
    # convert to spectral energy density dEdV_lambda
    y *= 4 * pi / c_light # [W/m^3/sr] -> [J/m^3/m]

    # convert to spectral number density: divide by photon energy
    y /= h_planck * c_light / x # [J/m^3/m] -> [1/m^3/m]

    # integrate over wavelength spectrum
    f = interp1d(x, y)
    n = quad(f, x[0], x[-1]) # numerical integral from x_min to x_max
    print n
    return n[0]


files = ['z0.txt', 'z0.2.txt', 'z0.4.txt', 'z0.6.txt', 'z1.txt', 'z2.txt', 'z3.txt', 'z4.txt']

s = zeros(8)
s2 = zeros(8)
for i,f in enumerate(files):
    s[i] = numberDensity('Kneiske2004_IRB/'+f)
    s2[i] = process1('Kneiske2004_IRB/'+f)
s /= s[0]
s2 /= s2[0]

z = array([0, 0.2, 0.4, 0.6, 1, 2, 3, 4])
s *= (1+z)**3
s2 *= (1+z)**3

figure()
plot(z, s, 'k')
plot(z, s2, 'r--')

xlim(0, 5)
xlabel('Redshift z')
ylabel('Overall Scaling')

#savefig('Kneiske2004_IRB_scaling.png',bbox_inches='tight')
show()
