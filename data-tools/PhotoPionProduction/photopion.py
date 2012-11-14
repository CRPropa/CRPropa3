from numpy import *

# This script reformats the pion production rate data files from CRPropa2.
# Mean free pathes larger than the Hubble length are omitted
lHubble = 4300. # Mpc
#
# The data files used are found in CRPropa/src/Interactions/proton_sophia/
# Copy the content or execute the script from within this folder.
#
# The used files are
# 'pionprodrate_n' : pion production rates of neutrons against the CMB
# 'pionprodrate_n_ir_y0' : pion production rates of neutrons against the IRB (Kneiske)
# 'pionprodrate_p' : pion production rates of protons against the CMB
# 'pionprodrate_p_ir_y0' : pion production rates of protons against the IRB (Kneiske)


# pion production rate on the CMB
eCMB, nCMB = genfromtxt('pionprodrate_n.txt', comments='#', unpack=True)
pCMB = genfromtxt('pionprodrate_p.txt', usecols=(1), comments='#')

fout1 = open('photopion_CMB.txt', 'w')
fout1.write('# Pion production rate for protons and neutrons on the CMB\n')
fout1.write('# Energy [EeV], Rate [1/Mpc]\n')

for E, p, n in zip(eCMB, pCMB, nCMB):
	if (p < 1/lHubble) and (n < 1/lHubble):
		continue
	fout1.write('%.5f %.5f %.5f\n'%( E, p, n ))

fout1.close()


# pion production rate on the IRB (Kneiske 2002)
eIRB, nIRB = genfromtxt('pionprodrate_n_ir_y0.txt', comments='#', unpack=True)
pIRB = genfromtxt('pionprodrate_p_ir_y0.txt', usecols=(1), comments='#')

fout2 = open('photopion_IRB.txt', 'w')
fout2.write('# Pion production rate for protons and neutrons on IRB\n')
fout2.write('# Energy [EeV], Rate [1/Mpc]\n')

for E, p, n in zip(eIRB, pIRB, nIRB):
	if (p < 1/lHubble) and (n < 1/lHubble):
		continue
	fout2.write('%.5f %.5f %.5f\n'%( E, p, n ))

fout2.close()
