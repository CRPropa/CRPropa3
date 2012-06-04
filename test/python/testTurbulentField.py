# evaluating functions for the turbulentField module
from pylab import *
from mpc import *

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

def getVectorFieldEnergySpectralDensity(Bx, By, Bz):
	'''
	calculate the energy spectral density for an 3-dimensional Vector field
	'''
	E = fftEnergySpectralDensity(Bx)
	E += fftEnergySpectralDensity(By)
	E += fftEnergySpectralDensity(Bz)
	n = shape(Bx)[0]
	K = fftfreq(n)
	# project E(kx,ky,kz) -> E(k)
	d = {}
	for ix in range(n):
		for iy in range(n):
			for iz in range(n):
				magK = (K[ix]**2 + K[iy]**2 + K[iz]**2)**.5
				energy = E[ix,iy,iz]
				d.setdefault(magK, [energy]).append(energy)
	k = d.keys()
	k.sort()
	k = array(k)
	Ek = zeros(len(k))
	for i,key in enumerate(k):
		Ek[i] = mean(d[key])
	return k, Ek



### create field
n = 128
lMin, lMax = 2, 32
Brms = 1
alpha = -11./3
field = TurbulentMagneticField(Vector3d(0, 0, 0), n, n)
field.initialize(lMin, lMax, Brms, alpha)
Lc = field.getCorrelationLength()

### copy field grid to array(s)
Bx, By, Bz = zeros((3, n, n, n))
for ix in range(n):
	for iy in range(n):
		for iz in range(n):
			b = field.getField(Vector3d(ix, iy, iz))
			Bx[ix, iy, iz] = b.x
			By[ix, iy, iz] = b.y
			Bz[ix, iy, iz] = b.z

### plot slice in position space
figure()
A = ((Bx**2 + By**2 + Bz**2)**.5)[:,:,n/2]
im = imshow(A, origin='lower', extent=[0,n,0,n], vmin=0, vmax=3)
cbar = colorbar(im)
cbar.set_label(r'$|\vec{B}(\vec{x})| / B_{rms}$')
xlabel('$x$')
ylabel('$y$')
text(0.8, 1.05, '$z=%i$'%(n/2), transform=gca().transAxes)
savefig('TurbulentField_slicePositionSpace.png', bbox_inches='tight')

### plot slice in configuration space
figure()
Bkx = fftshift(fftn(Bx))
Bky = fftshift(fftn(By))
Bkz = fftshift(fftn(Bz))
Bk = ((Bkx*Bkx.conjugate() + Bky*Bky.conjugate() + Bkz*Bkz.conjugate()).real)**.5
A = log10(Bk[:,:,n/2])
im = imshow(A, origin='lower', vmin=0)
cbar = colorbar(im)
cbar.set_label(r'$log_{10}(|\vec{B}(\vec{k})| / B_{rms})$')
xlabel('$k_x$')
ylabel('$k_y$')
k = fftshift(fftfreq(n))
idx = arange(0,n,n/4)
xticks(idx, k[idx])
yticks(idx, k[idx])
xlim(0,n)
ylim(0,n)
text(0.8, 1.05, '$k_z=%.2f$'%k[n/2], transform=gca().transAxes)
savefig('TurbulentField_sliceConfigurationSpace.png', bbox_inches='tight')

### plot slice and periodical extension in position space
figure()
A = zeros((3*n,3*n))
for i,j in ((0,1), (1,0), (1,1), (1,2), (2,1)):
	A[i*n:(i+1)*n, j*n:(j+1)*n] = Bx[:,:,n/2]
im = imshow(ma.masked_array(A, A == 0), origin='lower', extent=[-n,2*n,-n,2*n], vmin=0, vmax=3)
cbar = colorbar(im)
cbar.set_label(r'$|\vec{B_x}|/B_{rms}$')
xticks([-n, 0, n, 2*n])
yticks([-n, 0, n, 2*n])
xlabel('$x$')
ylabel('$y$')
text(0.8, 1.05, '$z=%i$'%(n/2), transform=gca().transAxes)
savefig('TurbulentField_slicePeriodicity.png', bbox_inches='tight')

### plot (2pt-auto-)correlation curves for various directions
figure()
corr = VectorFieldAutoCorrelation(Bx,By,Bz)
Lcs = []
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
			Lcs.append( corr.getIntegralLengthscale(step) )
			x,y = corr.getCorrelationCurve(step)
			plot(x,y,label=str(step))
xlabel('Distance')
ylabel('Normalized Autocorrelation')
xlim(0,32)
grid()
s = 'Correlation Length\n Nominal %.2f\n Simulated %.2f $\pm$ %.2f'%(Lc, mean(Lcs), std(Lcs)/(len(Lcs))**.5)
text(0.5, 0.95, s, ha='left', va='top', transform=gca().transAxes)
savefig('TurbulentField_coherenceLength.png', bbox_inches='tight')

### plot energy spectrum
figure()
k, Ek = getVectorFieldEnergySpectralDensity(Bx, By, Bz)
plot(k, Ek, label='Turbulent Spectrum')
i = k.searchsorted(1./lMax) + 2
plot(k, Ek[i]/k[i]**alpha * k**alpha, label='Slope $k^{-11/3}$')
axvline(1./lMax, color='r', linestyle='--', label='$k_{min}=1/%.1f$'%lMax)
axvline(1./lMin, color='r', linestyle='--', label='$k_{max}=1/%.1f$'%lMin)
loglog()
legend(loc='center left')
xlabel('Wavenumber $k$')
ylabel('Energy Spectral Density $E(k) [a.u.]$')
grid()
savefig('TurbulentField_spectrum.png', bbox_inches='tight')

### plot histogram of field strengths
figure()
Bx.resize(n**3)
By.resize(n**3)
Bz.resize(n**3)
hist(Bx, bins=40, range=(-3,3), histtype='step', normed=True, label='$B_x$', linewidth=2)
hist(By, bins=40, range=(-3,3), histtype='step', normed=True, label='$B_y$', linewidth=2)
hist(Bz, bins=40, range=(-3,3), histtype='step', normed=True, label='$B_z$', linewidth=2)
legend()
grid()
xlabel('$B/B_{RMS}$')
ylabel('Frequency')
brms = (mean( Bx**2 + By**2 + Bz**2 ))**.5
bmean = abs(mean(Bx + By + Bz)) 
text(0.95, 0.7, '$RMS$ = %.2f\nMean = %.2f'%(brms, bmean), ha='right', va='top', transform=gca().transAxes) 
savefig('TurbulentField_amplitude.png', bbox_inches='tight')
