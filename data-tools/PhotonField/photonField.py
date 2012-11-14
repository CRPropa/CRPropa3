from pylab import *
from mpc import h_planck, c_light, k_boltzmann, eV


def CMBSpectralRadiance(nu):
	# CMB spectral radiance [W/sr/m^2/Hz] at frequency nu [1/s], see http://en.wikipedia.org/wiki/Planck%27s_law
	return 2 * nu**3 * h_planck / c_light**2 / ( exp(nu * h_planck / k_boltzmann / 2.725) - 1 )

def CMBSpectralEnergyDensity(nu):
	# CMB spectral energy density [J/m^3/Hz] at frequency nu [1/s]
	return CMBSpectralRadiance(nu) * 4 * pi / c_light

def CMBSpectralDensity(nu):
	# CMB spectral number density [1/m^3/Hz] at frequency nu [1/s]
	return CMBSpectralEnergyDensity(nu) / (nu * h_planck)


def IRBSpectralRadiance(nu):
	# IRB (Kneiske et al. 2004) spectral radiance [W/sr/m^2/Hz] at frequency nu [1/s]
	
	# calculate wavelength [m]
	x = c_light / nu
	
	# load data points
	xKneiske, yKneiske = genfromtxt('Kneiske2004_IRB/z0.txt', unpack=1)
	y = interp(x * 1e6, xKneiske, yKneiske, 0, 0) * 1e-9 # lambda * I_lambda [W/m^2/sr]
	
	# remove scaling with lambda; I_lambda [W/m^3/sr]
	y /= x
	
	# convert to I_nu (using I_nu = - dlambda / dnu * I_lambda,   dnu = - c / lambda^2 dlambda)
	y *= x**2 / c_light
	return y

def IRBSpectralEnergyDensity(nu):
	# IRB (Kneiske et al. 2004) spectral energy density [J/m^3/Hz] at frequency nu [1/s]
	return IRBSpectralRadiance(nu) * 4 * pi / c_light

def IRBSpectralDensity(nu):
	# IRB (Kneiske et al. 2004) spectral density [1/m^3/Hz] at frequency nu [1/s]
	return IRBSpectralEnergyDensity(nu) / (nu * h_planck)


if __name__ == "__main__":
	nu1 = logspace(-5, -2.3) * eV / h_planck
	sed1 = CMBSpectralEnergyDensity(nu1)
	snd1 = CMBSpectralDensity(nu1)

	nu2 = logspace(-5, 1.5) * eV / h_planck
	sed2 = IRBSpectralEnergyDensity(nu2)
	snd2 = IRBSpectralDensity(nu2)

	# number of photons per m^3 = integral over spectral number density
	# integration using midpoint rule
	N1 = sum( (snd1[1:] + snd1[:-1]) / 2 * (nu1[1:] - nu1[:-1]) )
	N2 = sum( (snd2[1:] + snd2[:-1]) / 2 * (nu2[1:] - nu2[:-1]) )
	print N1 / 1e6, 'CMB-photons per cubic centimeter'
	print N2 / 1e6, 'IRB-photons per cubic centimeter'

	# plotting
	figure()
	plot(nu1, sed1, label='CMB')
	plot(nu2, sed2, label='IRB Kneiske')
	xlabel('Frequency [Hz]')
	ylabel('Spectral Energy Density [J/m$^{3}$/Hz]')
	loglog()
	legend(frameon=0)
	savefig('spectralEnergyDensity.png', bbox_inches='tight')

	figure()
	plot(nu1, snd1, label='CMB')
	plot(nu2, snd2, label='IRB Kneiske')
	xlabel('Frequency [Hz]')
	ylabel('Spectral Density [1/m$^{3}$/Hz]')
	loglog()
	legend(frameon=0)
	savefig('spectralDensity.png', bbox_inches='tight')

