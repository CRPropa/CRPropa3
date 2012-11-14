# This script combines the photodisintegration rate tables for CMB and IRB.

from numpy import *


fout = open('photodis_CMB_IRB.txt', 'w')
fout.write('#Z, N, disintegration channel (#n#p#H2#H3#He3#He4), disintegration rate [1/Mpc] (200 samples for log10(gamma) = 6-14)\n')


data1 = genfromtxt('photodis_CMB.txt', skip_header=1, unpack=True)
Z1, N1, channel1 = data1[0], data1[1], data1[2]

data2 = genfromtxt('photodis_IRB.txt', skip_header=1, unpack=True)
Z2, N2, channel2 = data2[0], data2[1], data2[2]

# Create unique identifiers for each species + channel
ID1 = Z1 * 1e8 + N1 * 1e6 + channel1
ID2 = Z2 * 1e8 + N2 * 1e6 + channel2


# Create a dictionary (id -> index) of ID1 for lookups
d = {}
for i, id1 in enumerate(ID1):
	d[id1] = i


# loop over and save all IRB data while adding the corresponding CMB data
for i, id2 in enumerate(ID2):
	rate = data2[3:, i]

	# find the CMB rate for the same channel if available
	if d.has_key(id2):
		j = d.pop(id2) # get and remove the dictionary entry
		rate += data1[3:, j] # add up the inverse mean free paths

	fout.write('%i\t%i\t%i'%(Z2[i], N2[i], channel2[i]))
	for r in rate:
		fout.write('\t'+str(r))
	fout.write('\n')


# write the remaining CMB rates
keys = d.keys()
keys.sort()
for id1 in keys:
	i = d[id1]
	fout.write('%i\t%i\t%i'%(Z1[i], N1[i], channel1[i]))
	rate = data1[3:, i]
	for r in rate:
		fout.write('\t'+str(r))
	fout.write('\n')
