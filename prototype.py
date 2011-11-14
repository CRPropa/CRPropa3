import numpy
import scipy.constants


### Object to handle simulation features
class FeatureHandler():
	def __init__(self):
		self.features = []

	def attach(self, f):
		if not f in self.features:
			self.features.append(f)

	def detach(self, f):
		if f in self.features:
			self.features.remove(f)

	def listFeatures(self):
		print self.features

	def propagate(self, p):
		for f in self.features:
			f.apply(p)


### Simulation features with common interface
class Feature():
	def __init__(self):
		pass
	def __repr__(self):
		return self.__class__.__name__
	def apply(self, p):
		pass

class OneDimPropagation(Feature):
	def apply(self, p):
		step =  p['nextStep']
		p['age'] += step
		p['step'] = step
		p['nextStep'] = 0.1

class EnergyLoss(Feature):
	def apply(self, p):
		p['energy'] *= .95

class EnergyBreakCondition(Feature):
	def __init__(self, Emin):
		self.Emin = Emin
	def apply(self, p):
		if p['energy'] <= self.Emin:
			p['state'] = -1

class TrajectoryLengthBreakCondition(Feature):
	def __init__(self, maxLength):
		self.maxLength = maxLength
	def apply(self, p):
		if p['age'] >= self.maxLength:
			p['state'] = -1

class NeutronDecay(Feature):
	def __init__():
		self.decayLength = 0

	def apply(self, p):
		if not((p['Z'], p['A']) == (0,1)):
			# not a neutron
			return

		if self.decayLength <= 0:
			# decayLength not yet set for this neutron
			self.decayLength = - tau * log(numpy.random.rand()) * getLorentzFactor(p) * scipy.constants.c

		# update decayLength
		self.decayLength -= p['step']

		if self.decayLength <= 0:
			# neutron decay
			(p['Z'], p['A']) = (1, 0)
			p['energy'] *= scipy.constants.m_p / scipy.constants.m_n
			return

		if p['nextStep'] > self.decayLength:
			# decay would happen before the next step, so reduce step size
			p['nextStep'] = self.decayLength
			return


### Particle
def newParticle(Z,A,m,E):
	p = {}
	p['state'] = 0
	p['step'] = 0.
	p['nextStep'] = 0.
	p['age'] = 0.
	p['Z'] = Z
	p['A'] = A
	p['mass'] = m
	p['energy'] = E
	return p

def getLorentzFactor(p):
	return p['energy'] / p['mass'] / c**2


### Simple simulation
def main():
	sim = FeatureHandler()
	sim.attach(EnergyBreakCondition(3))
	sim.attach(TrajectoryLengthBreakCondition(5))
	sim.attach(EnergyLoss())
	sim.attach(OneDimPropagation())

	print 'Forward tracking with the following features'
	sim.listFeatures()

	# create proton with 100 EeV
	p = newParticle(1,1,938e-12,100)
	p['step'] = 0.1
	p['nextStep'] = 0.1
	print 'propagating ..'
	while p['state'] <> -1:
#		print p
		sim.propagate(p)
	print 'done, final state:', p


if __name__ == '__main__':
	main()

