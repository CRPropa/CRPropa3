from numpy import *

# This script reformats the pion production rate data files from CRPropa2.
#
# The data files used are found in CRPropa/src/Interactions/proton_sophia/
# Copy the content or execute the script from within this folder.
# The used files are
# 'pionprodrate_n' : pion production rates of neutrons against the CMB
# 'pionprodrate_n_ir_y0' : pion production rates of neutrons against the IRB (Kneiske)
# 'pionprodrate_p' : pion production rates of protons against the CMB
# 'pionprodrate_p_ir_y0' : pion production rates of protons against the IRB (Kneiske)

# Interaction rates smaller than 1 / 1000 Gpc are not tabulated for the CMB.
# Even when accounting for cosmology these will not influence typical simulation results.
# For symmetry with the CMB rates, and for performance reasons the threshold is also applied to the IRB rates.
threshold = 1e-6

# pion production rate on the CMB
energies, nRate = genfromtxt('pionprodrate_n.txt', comments='#', unpack=True)
energies2, pRate = genfromtxt('pionprodrate_p.txt', comments='#', unpack=True)
if any(energies != energies2):
    print 'Mismatch of tabulated energies'

# save to one file, skip 0 entries
fout1 = open('photopion_CMB.txt', 'w')
fout1.write('# Pion production rate for protons and neutrons on the CMB\n')
fout1.write('# Energy [EeV], Rate [1/Mpc]\n')

for E, p, n in zip(energies, pRate, nRate):
    if (p < threshold) and (n < threshold):
        print 'CMB: skip', E, 'EeV:', p, n
        continue
    fout1.write('%7g %.6f %.6f\n'%( E, p, n ))

fout1.close()


# pion production rate on the IRB (Kneiske 2002)
energies, nRate = genfromtxt('pionprodrate_n_ir_y0.txt', comments='#', unpack=True)
energies2, pRate = genfromtxt('pionprodrate_p_ir_y0.txt', comments='#', unpack=True)
if any(energies != energies2):
    print 'Mismatch of tabulated energies'

# save to one file, skip 0 entries
fout2 = open('photopion_IRB.txt', 'w')
fout2.write('# Pion production rate for protons and neutrons on IRB\n')
fout2.write('# Energy [EeV], Rate [1/Mpc]\n')

for E, p, n in zip(energies, pRate, nRate):
    if (p < threshold) and (n < threshold):
        print 'IRB: skip', E, 'EeV:', p, n
        continue
    fout2.write('%8.6g %.6f %.6f\n'%( E, p, n ))

fout2.close()
