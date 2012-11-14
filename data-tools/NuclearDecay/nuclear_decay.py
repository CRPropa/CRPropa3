from numpy import *
import mpc

# Script to preprocess the nuclear decay data table from the BNL NuDat2 database
# See http://www.nndc.bnl.gov/nudat2/indx_sigma.jsp
# Query settings:
# Z: 0-36, output: formatted file

class Decay:
	def __init__(self):
		pass

	def __repr__(self):
		return 'Z=%i N=%i mode=%s tau=%.1e br=%.2f'%(self.Z, self.N, self.mode, self.tau, self.br)

	def load(self, s):
		l = s.split('\t')

		# Z, N, id
		self.Z = int(l[2])
		self.N = int(l[3])
		self.id = self.Z * 1000 + self.N

		# decay time
		s = l[9].strip()
		if s == 'infinity':
			self.tau = inf
		elif s == '':
			self.tau = 0
		else:
			self.tau = float(s) / log(2)

		# mode
		self.mode = l[12].strip()

		# branching ratio
		s = ''.join(c for c in l[13] if c not in ('>','<','=','~','%',' ','?','\n'))
		self.brString = s
		if s == '':
			self.br = 0.
		else:
			self.br = float(s) / 100.

	def isStable(self):
		return self.tau == inf

	def isBetaPlus(self):
		return self.mode.find('E') > -1

	def isBetaMinus(self):
		return self.mode.find('B') > -1


### parse data file
print '\nParsing data file'
print '-------------------------------------'
fin = open('NuDat2.txt','r')
lines = fin.readlines()
fin.close()

decayTable = [[[] for n in range(31)] for z in range(27)]

for line in lines[1:-3]:
	d = Decay()
	d.load(line)
	# skip if Z > 26 (Fe-56)
	if d.Z > 26:
		continue
	# skip if N > 30 (Fe-56)
	if d.N > 30:
		continue
	# skip if isomeric transition
	if d.mode == 'IT':
		print d, '<- skip (isomeric transition)'
		continue
	# skip if lifetime missing
	if d.tau == 0:
		print d, '<- skip (missing lifetime)'
		continue
	# skip if decay mode missing
	if d.mode == '':
		if not(d.isStable()):
			print d, '<- skip (missing decay mode)'
			continue
	# else store in decay table
	decayTable[d.Z][d.N].append(d)


### remove duplicate decays
print '\n\nRemoving duplicates'
print '-------------------------------------'
for z in range(27):
	for n in range(31):
		dList = decayTable[z][n]

		if len(dList) < 2:
			continue

		for i, d1 in enumerate(dList):
			for d2 in dList[i+1:]:
				if d1.mode == d2.mode:
					print d1
					print d2, ' <- remove \n'
					dList.remove(d2)


### explicitly edit some entries
print '\nExplicitly editing certain entries'
print '-------------------------------------'

# remove Li-5 alpha decay (equivalent to existing proton emission)
d0 = decayTable[3][2][0]
d1 = decayTable[3][2][1]
print d0
print d1, ' <- remove (equivalent to neutron emission)\n'
decayTable[3][2].remove(d1)

# remove He-5 alpha decay (equivalent to existing neutron emission)
d0 = decayTable[2][3][0]
d1 = decayTable[2][3][1]
print d1
print d0, ' <- remove (equivalent to neutron emission)\n'
decayTable[2][3].remove(d0)

# modify B-12 "B3A" decay to "B2A" as it would leave an empty nucleus
d = decayTable[5][7][1]
print d, ' <- change decay mode to B2A\n'
d.mode = 'B2A'

# Fe-45: to make beta+ decays exclusive
d = decayTable[26][19][0]
print d, ' <- set branching ratio to 0 (ratio equal to sum of following ratios)'
d.br = 0
brSum = 0
for d in decayTable[26][19][1:]:
	print d
	brSum += d.br
for d in decayTable[26][19][1:]:
	d.br /= brSum


### calculate exclusive mean life times
print '\n\nCalculating exclusive life times'
print '-------------------------------------'
for z in range(27):
	for n in range(31):
		dList = decayTable[z][n]

		# skip for 0 or 1 entry
		if len(dList) < 2:
			continue

		# get sum of branching ratios
		brSum = 0
		for d in dList:
			brSum += d.br

		# if sum is 0, set branching ratios to equal values
		if brSum == 0:
			for d in dList:
				d.br = 1. / len(dList)

		# else if sum not 1, search for an inclusive decay and/or normalize the branching ratios
		elif brSum != 1.:
			dInclusive = None
			brSumExclusive = 0
			for i,d in enumerate(dList):
				if d.br == 1.0:
					dInclusive = d # inclusive decay found
				else:
					brSumExclusive += d.br # add exclusive branching ratio

			if dInclusive != None:
				if dInclusive.br <= brSumExclusive:
					dList.remove(dInclusive) # remove if purely inclusive
				else:
					dInclusive.br -= brSumExclusive # else make exclusive

			# normalize all branching ratios
			for d in dList:
				d.br /= brSum

		# finally, calculate exclusive decay time by dividing with branching ratio, while removing zeros
		for d in dList:
			if d.br == 0:
				print d, ' <- remove (branching ratio 0)'
				dList.remove(d)
			else:
				d.tau /= d.br


### correct for electron capture contribution in beta+ decays
print '\nBeta+ correction'
print '-------------------------------------'
a0 = 5.29177e-11 # Bohr radius [m]
a0 /= mpc.c_light * (mpc.h_planck / 2 / pi) / mpc.eplus # convert to [1/eV]
Qe = mpc.mass_electron * mpc.c_squared / mpc.eV # [eV]

def I1(Q):
	x = Q / Qe
	return Qe**5 * ((2*x**4 - 9*x**2 - 8) * sqrt(x**2-1) + x*log(x+sqrt(x**2-1))) / 15

for z in range(27):
	for n in range(31):
		for d in decayTable[z][n]:
			if not(d.isBetaPlus()):
				continue

			m1 = mpc.getNucleusMass(mpc.getNucleusId(z+n, z))
			m2 = mpc.getNucleusMass(mpc.getNucleusId(z+n, z - 1))
			Q = (m1 - m2) * mpc.c_squared / mpc.eV

			Qbeta = (Q - Qe) # [eV]
			Qec = (Q + Qe) # [eV]

			# check if energetically possible
			if Qbeta < 0:
				print d, ' <- make stable (beta+ decay not possible)'
				d.tau = inf
				continue

			# see Basdevant, Fundamentals in Nuclear Physics, 4.3.2 and 4.3.3
			# ratio tau_beta+ / tau_ec
#			f = 2 / pi**2 * (z/a0)**3 * Qec**2 / Q**5 * 30
			f = 2 / pi**2 * (z/a0)**3 * Qec**2 / I1(Q)
			print d, ' <- beta+ correction %.1e'%f
			d.tau *= 1 + f


### set proton / neutron dripping with very short life time for all other isotopes
print '\n\nSet proton / neutron dripping for all other isotopes'
print '-------------------------------------'
for z in range(0,27):
	for n in range(0,31):
		if (z + n)==0:
			continue
		
		dList = decayTable[z][n]

		if len(dList) > 0:
			continue

		# else set p/n dripping
		d = Decay()
		d.Z = z
		d.N = n
		d.tau = 1e-99
		d.br = 1.
		if z > n: # proton dripping
			d.mode = 'P'
		else: # neutron dripping
			d.mode = 'N'
		dList.append(d)


### write to file
print '\n\nWrite to file'
print '-------------------------------------'
fout = open('nuclear_decay.txt','w')
fout.write('# Z, N, Mean Life Time [s], Decay Mode (#beta- #beta+ #alpha #p #n), dE\n')

# decay mode codes: #beta- #beta+ #alpha #p #n
modeDict = {'STABLE' : '0',
		'N'  : '00001',
		'2N' : '00002',
		'P'  : '00010',
		'2P' : '00020',
		'A'  : '00100',
		'2A' : '00200',
		'B-' : '10000',
		'2B-': '20000',
		'BN' : '10001',
		'B2N': '10002',
		'B3N': '10003',
		'B4N': '10004',
		'BNA': '10101',
		'BA' : '10100',
		'B2A': '10200',
		'B3A': '10300',
		'EC' : '01000',
		'2EC': '02000',
		'EA' : '01100',
		'EP' : '01010',
		'E2P': '01020',
		'E3P': '01030'}

for z in range(0,27):
	for n in range(0,31):
		if (z + n)==0:
			continue

		for d in decayTable[z][n]:
			# skip stable
			if d.tau == inf:
				continue
			fout.write('%i %i %s %e\n'%(d.Z, d.N, modeDict[d.mode], d.tau))

fout.close()
print 'done'

