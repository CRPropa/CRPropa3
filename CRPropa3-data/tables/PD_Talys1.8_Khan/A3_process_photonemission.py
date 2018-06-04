from numpy import *
import os


isotopes = genfromtxt('isotopes.txt')
path = './'  # path to simulated TALYS output


# -----------------------------------------------------------------------------
# Step 1: Calculate total cross section for all photodisintegration channels
# from one mother to one daughter isotope, e.g.
# sigma_tot(N-14->C-12) = sigma(N-14->C-12+p+n) + sigma(N-14->C-12+d) + ...
# -----------------------------------------------------------------------------
def parseChannel(channel):
    """
    Parse TALYS channel and return number of emitted protons / neutron (dZ, dN).
    channel : int/float/string (#n, #p, #d, #t, #He-3, #He-4)
    """
    s = '%06d' % int(channel)
    dZ = int(s[1]) + int(s[2]) +   int(s[3]) + 2*int(s[4]) + 2*int(s[5])
    dN = int(s[0]) + int(s[2]) + 2*int(s[3]) +   int(s[4]) + 2*int(s[5])
    return dZ, dN

# obtain (Z, N, Zdaughter, Ndaughter) for all selected photodisintegration channels
data = genfromtxt('xs_pd_thin.txt', dtype=[('Z',int), ('N',int), ('channel',int), ('xs','301f8')])
Z, N, channel, xs = data['Z'], data['N'], data['channel'], data['xs']

# Calculate daughters from number of emitted Z, N
dZ, dN = array([parseChannel(c) for c in channel]).T
Zd = Z - dZ
Nd = N - dN

# find unique tuples (Z,N,Zd,Nd)
A = c_[Z, N, Zd, Nd]
A = A.view(A.dtype.descr * 4)  # create record array of A
Aunique = unique(A)

# calculate total cross section for from mother (Z,N) to daughter (Zd,Nd)
xs_sum_daughter = array([sum(xs[ravel(a == A)], axis=0) for a in Aunique])

savetxt('xs_photon_sum.txt',
    c_[Aunique.view(int).reshape(Aunique.shape + (-1,)), xs_sum_daughter],
    fmt='%i\t%i\t%i\t%i' + '\t%g'*301,
    header='Total disintegration cross section from (Z,N) to (Zd,Nd)\nZ\tN\tZ_daughter\tN_daughter\tEgamma [MeV]\txs [mb] for eps = 0.2 - 200 MeV in steps of logE = 0.01')


# -----------------------------------------------------------------------------
# Step 2: Collect the disintegration cross sections with photon emission --> xs_photon.txt
# - Parse all gam....tot files holding the cross sections for (Z,N) --> (Zd, Nd) + Egamma
# - Discard those for which the disintegration channel was previously thinned out
# -----------------------------------------------------------------------------
print ("Collect disintegration cross sections with photon emission")

# output file for cross sections (Z,N) --> (Zd, Nd) + Egamma
fout = open('xs_photon.txt', 'w')
fout.write('# Z\tN\tZ_daughter\tN_daughter\tEgamma [MeV]\txs [mb]\n')
fout.write('#cross sections [mb] for emission of photon with Egamma [MeV] for incident background photon energies eps = 0.2 - 200 MeV in steps of logE = 0.01\n')
fmt = '%i\t%i\t%i\t%i\t%.4f' + '\t%.4g'*301 + '\n'  # output format

# parse all gam....tot files for each isotope
for (Z,N,A) in isotopes:
    print (Z, N)
    folder = path + '%i-%i/' % (Z, A)

    for fname in os.listdir(folder):
        # consider the gam....tot files only
        if not( fname.startswith('gam') and fname.endswith('tot') ):
            continue

        # parse daughter nucleus from filename
        s = fname.strip('gam.tot')
        Zd = int(s[0:3])
        Nd = int(s[3:6]) - Zd

        # check if channel survived the previous thinning
        if not array((Z, N, Zd, Nd), dtype=Aunique.dtype) in Aunique:
            continue

        # save photon energy and cross section
        Egamma = float(open(folder + fname).readline().split()[-1])
        xs = genfromtxt(folder + fname, usecols=1)
        fout.write( fmt % ((Z, N, Zd, Nd, Egamma) + tuple(xs)) )

fout.close()


# -----------------------------------------------
# Step 3: Apply thinning to photon emission cross sections --> xs_photon_thin.txt
# - Determine branching ratios and only select photon energies with at least 5% branching ratio at any incident photon energy
# - Renormalize the remaining ratios
# -----------------------------------------------
data = genfromtxt('xs_photon.txt')
select = zeros(len(data), dtype=bool)

for i, (Z,N,Zd,Nd) in enumerate(Aunique):
    # find all cross sections for (Z,N,Zd,Nd)
    idx = (data[:,0] == Z) * (data[:,1] == N) * (data[:,2] == Zd) * (data[:,3] == Nd)
    if sum(idx) == 0:
        continue  # no photon emission for this (Z,N) -> (Zd,Nd)

    # total cross section
    xs = data[idx,5:]
    xs_sum = sum(xs, axis=0)

    # calculate branching ratios
    br = nan_to_num(xs / xs_sum)

    # find channels with branching ratios above the 5% threshold
    s = amax(br, axis=1) > 0.05

    # renormalize cross section of surviving channels
    xs_sum_thin = sum(xs[s], axis=0)
    data[idx,5:] *= nan_to_num(xs_sum / xs_sum_thin)

    # remember selected channels
    select[idx] = s

print ('Selected %i of %i channels' %  (sum(select) ,len(data)))
savetxt('xs_photon_thin.txt', data[select], fmt='%i\t%i\t%i\t%i\t%.4f' + '\t%g'*301)
