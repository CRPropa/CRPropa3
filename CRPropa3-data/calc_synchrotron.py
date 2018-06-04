from numpy import *
from scipy.special import kv
from scipy import integrate

def synchrotron_function(xval):
		"""
		Calculate cumulative synchrotron function based on J.D. Jackson (p. 785, formula 14.91)
		F(x) = (\int_{x}^{\infinity} K_{5/3}(x') dx') 
		for x = E_gamma / E_critical, where E_critical = hbar * 3/2 * c/rho * (E/(mc**2))**2
		E_gamma : Energy of synchrotron photon
		E       : Energy of particle
		rho     : gyroradius 
		Returns : cumulative synchrotron function
		"""
		F = zeros(len(xval))
		for i, value in enumerate(xval):
			a = xval[i:]
			F[i] = integrate.trapz(x=a, y=kv(5./3.,a))
		for i,value in enumerate(xval):
			b = integrate.cumtrapz(x=xval, y=xval*F, initial = 0)
		return b

# ----------------------------------------------------------------
# Cumulative differential synchrotron spectrum (comp. J.D. Jackson(p. 785, formula 14.91))
# ----------------------------------------------------------------

x = logspace(-6, 2, 801)
cdf = synchrotron_function(x)
lx = log10(x)
data = c_[lx, cdf]
fname  = 'data/Synchrotron/spectrum.txt'
header = 'x\t: fraction synchrotron photon frequency to critical frequency\nlog10(x)\tCDF\n'
fmt = '%.5e\t%.5e'
savetxt(fname, data, fmt=fmt, header=header)
