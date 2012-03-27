from pylab import *
from mpc import *


# 100 EeV proton
c = Candidate()
c.current.setId(getNucleusId(1,1))
c.current.setEnergy(100 * EeV)
c.current.setDirection(Vector3(1, 0, 0))

# uniform perpendicular magnetic field of 10 nG
field = UniformMagneticField(Vector3(0,0,1e-12))

# resulting gyroradius
R = c.current.getMomentum().mag() / c.current.getCharge() / 1e-12


def propagate(tolerance):
	prop = DeflectionCK(field, tolerance, DeflectionCK.WorstOffender, 0.1 * kpc)
	maxLen = MaximumTrajectoryLength(2 * pi * R)

	c.setTrajectoryLength(0)
	c.setCurrentStep(0)
	c.setNextStep(1 * kpc)
	c.current.setPosition(Vector3(0, R, 0))
	c.current.setDirection(Vector3(1, 0, 0))
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
		dirDeviation.append(c.current.getDirection().dot(Vector3(cos(t), -sin(t), 0)))

	return array(posX), array(posY), array(dirX), array(dirY), array(dirDeviation), array(theta), nSteps



### Plot trajectory in x-space with 1e-4 tolerance
figure()
posX, posY, dirX, dirY, dirDev, theta, n = propagate(1e-4)
plot(posX/R, posY/R, label='Simulated')
plot(sin(theta), cos(theta),'k--',label='True')
title('Trajectory using $10^{-4}$ tolerance')
legend()
xlim(-1.5,1.5)
ylim(-1.5,1.5)
xlabel(r'$x / R_L$')
ylabel(r'$y / R_L$')
grid()
savefig('DeflectionCK_xtrajectory', bbox_inches='tight')


### Plot trajectory in p-space with 1e-4 tolerance
figure()
posX, posY, dirX, dirY, dirDev, theta, n = propagate(1e-4)
plot(dirX, dirY, label='Simulated')
plot(sin(theta), cos(theta),'k--',label='True')
title(r'Trajectory in $\vec{p}$-space using $10^{-4}$ tolerance')
legend()
xlim(-1.5,1.5)
ylim(-1.5,1.5)
xlabel(r'$p_x / |p_x|$')
ylabel(r'$p_y / |p_y|$')
grid()
savefig('DeflectionCK_ptrajectory', bbox_inches='tight')


### Directional error as function of distance for different tolerances
figure()
for tolerance in [1e-3, 1e-4, 1e-5, 1e-6]:
	posX, posY, dirX, dirY, dirDev, theta, n = propagate(tolerance)
	plot(theta/2/pi, arccos(dirDev)*180/pi, label='%.0e, %i'%(tolerance, n))
legend(title='Tolerance, Steps', loc='upper left')
xlabel(r'Travelled Distance / $2 \pi R_L$')
ylabel(r'Direction Error [$^\circ$]')
xlim(0,1)
grid()
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
xlim(0,1)
grid()
savefig('DeflectionCK_xdeviation.png',bbox_inches='tight')


