from pylab import *
from mpc import *


# 100 EeV proton
c = Candidate()
c.current.setId(getNucleusId(1,1))
c.current.setEnergy(100 * EeV)
c.current.setDirection(Vector3d(1, 0, 0))

# uniform perpendicular magnetic field of 10 nG
field = UniformMagneticField(Vector3d(0,0,1e-12))

# resulting gyroradius
R = c.current.getMomentum().mag() / c.current.getCharge() / 1e-12


def propagate(tolerance):
	prop = DeflectionCK(field, tolerance, 0.1 * kpc)
	maxLen = MaximumTrajectoryLength(2 * pi * R)

	c.setTrajectoryLength(0)
	c.setCurrentStep(0)
	c.setNextStep(1 * Mpc) # set a large initial step, so that an initial acceleration is uneccessary
	c.current.setPosition(Vector3d(0, R, 0))
	c.current.setDirection(Vector3d(1, 0, 0))
	c.setActive(True)

	posX, posY, dirX, dirY, dirDeviation, theta = [], [], [], [], [], []
	nSteps = 0
	while c.isActive():
		prop.process(c)
		maxLen.process(c)
		nSteps +=1

		posX.append(c.current.getPosition().x())
		posY.append(c.current.getPosition().y())

		dirX.append(c.current.getDirection().x())
		dirY.append(c.current.getDirection().y())

		t = c.getTrajectoryLength() / R
		theta.append(t)
		dirDeviation.append(c.current.getDirection().dot(Vector3d(cos(t), -sin(t), 0)))

	return array(posX), array(posY), array(dirX), array(dirY), array(dirDeviation), array(theta), nSteps



### Plot trajectory in x-space with 1e-4 tolerance
figure()
posX, posY, dirX, dirY, dirDev, theta, n = propagate(1e-4)
plot(posX/R, posY/R, "o", label='Simulated')
plot(sin(linspace(0,2*pi)), cos(linspace(0,2*pi)),'k--',label='True')
legend()
xlim(-1.5,1.5)
ylim(-1.5,1.5)
xlabel(r'$x / R_L$')
ylabel(r'$y / R_L$')
grid()
text(0.05, 0.95, 'Tolerance 1e-4', ha='left', transform=gca().transAxes)
savefig('DeflectionCK_xtrajectory', bbox_inches='tight')

### Plot trajectory in p-space with 1e-4 tolerance
figure()
posX, posY, dirX, dirY, dirDev, theta, n = propagate(1e-4)
plot(dirX, dirY, "o", label='Simulated')
plot(sin(linspace(0,2*pi)), cos(linspace(0,2*pi)),'k--',label='True')
legend()
xlim(-1.5,1.5)
ylim(-1.5,1.5)
xlabel(r'$p_x / |p_x|$')
ylabel(r'$p_y / |p_y|$')
grid()
text(0.05, 0.95, 'Tolerance 1e-4', ha='left', transform=gca().transAxes)
savefig('DeflectionCK_ptrajectory', bbox_inches='tight')

### Directional error as function of distance for different tolerances
figure()
for tolerance in [1e-3, 1e-4, 1e-5, 1e-6]:
	posX, posY, dirX, dirY, dirDev, theta, n = propagate(tolerance)
	plot(theta/2/pi, arccos(dirDev)*180/pi, label='%.0e, %i'%(tolerance, n))
legend(title='Tolerance, Steps', loc='upper left')
xlabel(r'Travelled Distance / $2 \pi R_L$')
ylabel(r'Direction Error [$^\circ$]')
grid()
semilogy()
savefig('DeflectionCK_pdeviation.png',bbox_inches='tight')

### Positional error as function of distance for different tolerances
figure()
for tolerance in [1e-3, 1e-4, 1e-5, 1e-6]:
	posX, posY, dirX, dirY, dirDev, theta, n = propagate(tolerance)
	displacement = ((posX/R - sin(theta))**2 + (posY/R - cos(theta))**2 )**.5
	plot(theta/2/pi, displacement, label='%.0e, %i'%(tolerance, n))
legend(title='Tolerance, Steps', loc='upper left')
xlabel(r'Travelled Distance / $2 \pi R_L$')
ylabel(r'Position Error / $R_L$')
grid()
semilogy()
ylim(1e-8,1e-1)
savefig('DeflectionCK_xdeviation.png',bbox_inches='tight')


