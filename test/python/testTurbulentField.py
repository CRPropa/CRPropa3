# evaluating functions for the turbulentField module
from pylab import *
from mpc import *
from scipy import optimize


def fftAutoCorrelation(a):
	'''
	calculate 2point-autocorrelation for an n-dimensional array of real input
	'''
	b = rfftn(a)
	return irfftn(b * b.conjugate()).real

class VectorFieldAutoCorrelation():
	'''
	2point-autocorrelation for a field of 3-vectors
	'''
	def __init__(self, Bx, By, Bz):
		self.n = shape(Bx)[0]
		self.Rx = fftAutoCorrelation(Bx)
		self.Ry = fftAutoCorrelation(By)
		self.Rz = fftAutoCorrelation(Bz)

	def getCorrelationCurve(self, step=(1, 1, 1)):
		'''
		returns the correlation curve in the direction of [step], summed over all field components
		'''
		step = array(step)
		# number of steps that can be walked accross the field
		nSteps = int(self.n / 2 / max(abs(step)))
		x = arange(0, nSteps) * norm(step)
		y = zeros(nSteps)
		for i in range(nSteps):
			y[i] += self.Rx[ i * step[0], i * step[1], i * step[2] ]
			y[i] += self.Ry[ i * step[0], i * step[1], i * step[2] ]
			y[i] += self.Rz[ i * step[0], i * step[1], i * step[2] ]
		y /= y[0]
		return x, y
	
	def plotCorrelationCurve(self, step=(1, 1, 1)):
		x, y = self.getCorrelationCurve(step)
		fig = figure()
		plot(x, y, label=str(step))
		legend()
		grid()
		xlabel(r'$|\vec{d}| \left[\rm{gridpoints}\right]$')
		ylabel(r'$R(\vec{d})/R(0)$')
		return fig

	def getIntegralLengthscale(self, step=(1, 1, 1)):
		'''
		returns the integral lengthscale in the direction of [step]
		'''
		x, y = self.getCorrelationCurve(step)
		# use symmetry: Lc = \int_{-inf}^{inf} R/R_0 dx = 2*\int_{0}^{inf} R/R_0 dx
		return (2 * sum(y[1:]) + y[0]) * norm(step)


def fftEnergySpectralDensity(a):
	'''
	calculate the energy spectral density for an n-dimensional array
	'''
	b = fftn(a)
	return (b * conjugate(b)).real


class VectorFieldEnergySpectralDensity():
	def __init__(self, Bx, By, Bz):
		# sum of the energy spectral density of each component
		E = fftEnergySpectralDensity(Bx)
		E += fftEnergySpectralDensity(By)
		E += fftEnergySpectralDensity(Bz)

		n = shape(Bx)[0]
		K = fftfreq(n)

		# project E(kx,ky,kz) -> E(k) in a dictionary
		d = {}
		for ix in range(n):
			for iy in range(n):
				for iz in range(n):
					magK = (K[ix]**2 + K[iy]**2 + K[iz]**2)**.5
					energy = E[ix,iy,iz]
					d.setdefault(magK,[energy]).append(energy)
		
		# make two arrays out of the dictionary
		keys = d.keys()
		keys.sort()
		self.k = zeros(len(keys))
		self.Ek = zeros(len(keys))
		for i,key in enumerate(keys):
			self.k[i] = key
			self.Ek[i] = mean(d[key])

		# normalize Ek
		self.Ek /= max(self.Ek)


	def findKmin(self, threshold=1e-6):
		for i in range(len(self.Ek)):
			if self.Ek[i] > threshold:
				return self.k[i]

	def findKmax(self, threshold=1e-6):
		for i in range(len(self.Ek)-1,-1,-1):
			if self.Ek[i] > threshold:
				return self.k[i]

	def powerLaw(self, p, x):
		return (p[1]*x)**(p[0])

	def costFunction(self, p, x, y):
		return self.powerLaw(p, x) - y

	def fitPowerLaw(self):
		# determine fitting range
		kmin=self.findKmin()
		kmax=self.findKmax()
		iMin = self.k.searchsorted(kmin)
		iMax = self.k.searchsorted(kmax)
		x = self.k[iMin:iMax]
		y = self.Ek[iMin:iMax]
		# fit
		p0 = [-11./3., 1.] # initial guess
		p1, success = optimize.leastsq(self.costFunction, p0, args=(x, y))
		return p1

	def plot(self):
		plot(self.k, self.Ek, label='Sim. Turbulence')
		p1 = self.fitPowerLaw()
		alabel = 'Fit: $k^{%.2f/3}$'%(p1[0]*3)
		plot(self.k, self.powerLaw(p1, self.k), label=alabel)
		kmin = self.findKmin()
		kmax = self.findKmax()
		axvline(kmin, color='r',linestyle='--',label='$k_{min}=1/'+str(1/kmin)+'$')
		axvline(kmax, color='r',linestyle='--',label='$k_{max}=1/'+str(1/kmax)+'$')
		loglog()
		legend(loc=0)
		xlabel('Wavenumber $k$')
		ylabel('Energy Spectral Density $E(k)$')
		grid()
		

def retrieveField(field):
	n = field.getGridSamples()
	Bx, By, Bz = zeros((3, n, n, n))
	for ix in range(n):
		for iy in range(n):
			for iz in range(n):
				b = field.getField(Vector3(ix, iy, iz))
				Bx[ix, iy, iz] = b.x()
				By[ix, iy, iz] = b.y()
				Bz[ix, iy, iz] = b.z()
	return Bx, By, Bz


#if __name__ == '__main__':
n = 64
field = TurbulentMagneticFieldGrid(Vector3(0, 0, 0), n, 1, 2, 16, 1, -11. / 3)
Bx, By, Bz = retrieveField(field)


### periodicity
figure()
subplot(111, aspect='equal')
A = zeros((3*n,3*n))
for i,j in ((0,1), (1,0), (1,1), (1,2), (2,1)):
	A[i*n:(i+1)*n, j*n:(j+1)*n] = Bx[:,:,32]
pc = pcolor(ma.masked_array(A, A == 0))
cbar = colorbar(pc)
cbar.set_label(r'$|\vec{B_x}|/B_{rms}$')
xlim(0,3*n)
ylim(0,3*n)
xlabel(r'$x$ [gridpoints]')
ylabel(r'$y$ [gridpoints]')
savefig('TurbulentField_periodicity.png', bbox_inches='tight')


### slice in configuration space
Bkx = fftshift(fftn(Bx))
Bky = fftshift(fftn(By))
Bkz = fftshift(fftn(Bz))
Bk = ((Bkx*Bkx.conjugate() + Bky*Bky.conjugate() + Bkz*Bkz.conjugate())**.2).real
del Bkx, Bky, Bkz
figure()
subplot(111, aspect='equal')
pc = pcolor(Bk[:,:,n/2])
cbar = colorbar(pc)
cbar.set_label(r'$|\vec{B}(\vec{k})|$ [a.u.]')
xlabel(r'$\vec{k}_x$')
ylabel(r'$\vec{k}_y$')
k = fftshift(fftfreq(n))
idx = arange(0,n,n/4)
xticks(idx, k[idx])
yticks(idx, k[idx])
xlim(0,n)
ylim(0,n)
savefig('TurbulentField_configurationSpace.png', bbox_inches='tight')


### correlation length + isotropy
figure()
corr = VectorFieldAutoCorrelation(Bx,By,Bz)
Lc = []
steps = []
for ix in arange(-2,3):
	for iy in arange(-2,3):
		for iz in arange(-2,3):
			if ix == 0 and iy == 0 and iz == 0:
				continue
			if (-ix,-iy,-iz) in steps:
				continue
			step = (ix,iy,iz)
			# steps.append(step)
			Lc.append( corr.getIntegralLengthscale(step) )
			x,y = corr.getCorrelationCurve(step)
			plot(x,y,label=str(step))
xlabel('Distance [gridpoints]')
ylabel('Normalized Autocorrelation')
xlim(0,32)
grid()
s = 'Correlation Length\n Nominal %.2f\n Simulated %.2f $\pm$ %.2f'%(field.getCorrelationLength(), mean(Lc), std(Lc)/(len(Lc))**.5)
text(0.5, 0.95, s, ha='left', va='top', transform=gca().transAxes)
savefig('TurbulentField_correlation.png', bbox_inches='tight')


### energy spectrum
esd = VectorFieldEnergySpectralDensity(Bx, By, Bz)
figure()
esd.plot()
savefig('TurbulentField_spectrum.png', bbox_inches='tight')


### field strength, mean and brms
Bx.resize(n**3)
By.resize(n**3)
Bz.resize(n**3)
figure()
hist(Bx, bins=40, range=(-3,3), histtype='step', normed=True, label='$B_x$', linewidth=2)
hist(By, bins=40, range=(-3,3), histtype='step', normed=True, label='$B_y$', linewidth=2)
hist(Bz, bins=40, range=(-3,3), histtype='step', normed=True, label='$B_z$', linewidth=2)
legend()
grid()
xlabel('Magnetic Field Amplitude$')
ylabel('Frequency')
Brms = (mean( Bx**2 + By**2 + Bz**2 ))**.5
text(1.45, 0.5, '$B_{RMS}$ = %.2f'%(Brms)) 
savefig('TurbulentField_amplitude.png', bbox_inches='tight')
