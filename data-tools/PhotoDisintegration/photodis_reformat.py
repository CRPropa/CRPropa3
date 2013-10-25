from numpy import *

# This script reformats the exclusive mean free path data files from CRPropa2.
# - The four CRPropa2 files are merged into one,
# - nuclides with Z>26 or N>30 are omitted and
#
# The folders with the date files being used, are found in [CRPropa2 source]/TabulatedTALYSMeanFreePath/.
# Copy the content or execute the script from within this folder.
#
# The used files are
# 'PDInitNucleusId.cmt' : Particle ID (Z*1000+A) and corresponding line numbers in 'PDExclTabMnFrPthCrossId.cmt'
# 'PDExclTabMnFrPthCrossId.cmt' : Disintegration channel (#n #p #H2 #H3 #He3 #He4) and data index in 'PDExclTabMnFrPthCross.cmt'
# 'PDExclTabMnFrPthCross.cmt' : Inverse mean free path in [1/Mpc] for 200 equidistant Lorentz-factors gamma = 10^6 - 10^14
#
# CMB, IRB (Kneiske)


def getAZ(s):
	# returns A, Z from the particle id (Z*1000+A)
	pid = int(s)
	return pid % 1000, pid // 1000


def getDigit(number, d):
	# returns the d-th digit, counted from the back
	return number % (10**d) // (10**(d-1))


def isBogus(Z, A, channel):
	# checks if the disintegration channel is impossible or leaves an empty nucleus
	nN = getDigit(channel, 6)
	nP = getDigit(channel, 5)
	nH2 = getDigit(channel, 4)
	nH3 = getDigit(channel, 3)
	nHe3 = getDigit(channel, 2)
	nHe4 = getDigit(channel, 1)

	N = A - Z

	dZ = nP + nH2 + nH3 + 2*nHe3 + 2*nHe4
	dA = nN + nP + 2*nH2 + 3*nH3 + 3*nHe3 + 4*nHe4
	dN = dA - dZ

	if (Z - dZ < 0):
		print 'Z\' < 0: skipping line', j
		print 'Z=%i N=%i dZ=%i dN=%i channel=%i'%(Z,N,dZ,dN,channel)+'\n'
		return True
	if (N - dN < 0):
		print 'N\' < 0: skipping line', j
		print 'Z=%i N=%i dZ=%i dN=%i channel=%i'%(Z,N,dZ,dN,channel)+'\n'
		return True
	if (A - dA <= 0):
		print 'A\' <= 0: skipping line', j
		print 'Z=%i N=%i dZ=%i dN=%i channel=%i'%(Z,N,dZ,dN,channel)+'\n'
		return True


def reformat(photonField):
	# input files
	data1 = genfromtxt(photonField+'/PDInitNucleusId.cmt', skip_header=1)
	data2 = genfromtxt(photonField+'/PDExclTabMnFrPthCrossId.cmt')
	data3 = genfromtxt(photonField+'/PDExclTabMnFrPthCross.cmt')
	# output file
	fout = open('photodis_'+photonField+'.txt','w')
	fout.write('#Z, N, disintegration channel (#n#p#H2#H3#He3#He4), disintegration rate [1/Mpc] (200 samples for log10(gamma) = 6-14)')

	for i in range(size(data1, 0)):
		pid = data1[i, 0]

		A, Z = getAZ(pid) # get particle charge- and mass number
		N = A - Z

		if Z>26 or N>30: # skip isotopes heavier than Fe-56
			continue

		j0 = int(data1[i, 1]) # start index for channels in 'PDExclTabMnFrPthCrossId.cmt'
		j1 = int(data1[i, 2]) # end index

		for j in range(j0, j1):
			channel = int(data2[j, 0])

			if isBogus(Z,A,channel):
				print ' .. in line', j
				continue

			fout.write('\n')
			fout.write(str(Z)+'\t'+str(N)+'\t') # write the particle's Z, N
			fout.write(str(channel)) # write disintegration channel

			k0 = int(data2[j, 1]) # start index in 'PDExclTabMnFrPthCross.cmt'

			for rate in data3[k0:k0 + 200]:
				fout.write('\t'+str(rate))


reformat('CMB')
reformat('IRB')

