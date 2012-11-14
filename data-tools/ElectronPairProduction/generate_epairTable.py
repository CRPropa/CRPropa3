from pylab import *
from mpc import *
from scipy import interpolate, integrate



def densityCMB(eps):
  # CMB spectral number density [1/m^3/J] at photon energy eps [J]
  return 8 * pi * eps**2 / h_planck**3 / c_light**3 / (exp(eps / k_boltzmann / 2.725) - 1)

xKneiske, yKneiske = genfromtxt('KneiskeIRB_z0.txt', unpack=1)
def densityIRB(eps):
  # IRB (Kneiske et al. 2004) spectral number density [1/m^3/J] at photon energy eps [J]
  x = h_planck * c_light / eps # wavelength [m]
  y = interp(x * 1e6, xKneiske, yKneiske, 0, 0) * 1e-9 # lambda * I_lambda [W/m^2/sr]
  y /= x # remove scaling with lambda; I_lambda [W/m^3/sr]
  y *= x**2 / h_planck / c_light # convert to I_eps (using I_eps = - dlambda / deps * I_lambda,   deps = - h*c / lambda^2 dlambda)
  y *= 4 * pi / c_light # convert to energy density
  y /= eps # convert to number density
  return y



### Calculate the energy loss -dE/dx due to electron pair production
# c.f. Blumenthal 1970, Phys. Rev. D Volume 1, Number 6, Page 1596-1602

# Tabulated values for equation (14) taken from figure 2
# xi = 2 * gamma * eps / (me*c^2), where eps is the photon energy, gamma the nucleus lorentz factor
# phi(xi) integral over differential cross sections
# This could be replaced with an explicit calculation of the integral
x, y = genfromtxt('Blumenthal_phi.txt', comments='#', unpack=1)
phi = interpolate.interp1d(10**x, 10**y)
xiMin = 10**x[0]
xiMax = 10**x[-1]

# prefactor of equation (13)
r0 = 2.817940e-15 # classical electron radius [m]
alpha = 7.297352e-3 # fine-structure constant
c = alpha * r0**2 * 1**2 * mass_electron**2 * c_light**4



### generate CMB table
# integrand of equation (13) with CMB
def f(xi, gamma):
  eps = xi * mass_electron * c_light**2 / 2 / gamma
  return densityCMB(eps) * phi(xi) / xi**2

N = 21
A1 = zeros((N, 2))

for i, g in enumerate(linspace(8.5, 12.5, N)):
  F, err = integrate.quad(f, xiMin, xiMax, 10**g)
  A1[i,1] = c * F
  A1[i,0] = g

savetxt('epair_CMB.txt', A1, header='log10(lorentzFactor)  dE/dx [J/m]', fmt=['%.1f','%.6e'], delimiter='\t')



### generate IRB table
def f(xi, gamma):
  eps = xi * mass_electron * c_light**2 / 2 / gamma
  return densityIRB(eps) * phi(xi) / xi**2

N = 31
A2 = zeros((N, 2))

for i, g in enumerate(linspace(6.5, 12.5, N)):
  F, err = integrate.quad(f, xiMin, xiMax, 10**g)
  A2[i,1] = c * F
  A2[i,0] = g

savetxt('epair_IRB.txt', A2, header='log10(lorentzFactor)  dE/dx [J/m]', fmt=['%.1f','%.6e'], delimiter='\t')



### generate combined table
A3 = A2.copy()
A3[10:,1] += A1[:,1]

savetxt('epair_CMB_IRB.txt', A3, header='log10(lorentzFactor)  dE/dx [J/m]', fmt=['%.1f','%.6e'], delimiter='\t')

