from numpy import *
import os


isotopes = genfromtxt('isotopes.txt')
path = './'  # path to simulated TALYS output


# -----------------------------------------------------------------------------
# Step 1: Collect all exclusive cross sections from TALYS --> xs_pd.txt
# -----------------------------------------------------------------------------
print ('Collecting TALYS output')

def parseChannel(channel):
    """
    Parse TALYS channel and return number of emitted protons / neutron (dZ, dN).
    channel : int/string (#n, #p, #d, #t, #He-3, #He-4)
    """
    s = '%06d' % int(channel)
    dZ = int(s[1]) + int(s[2]) +   int(s[3]) + 2*int(s[4]) + 2*int(s[5])
    dN = int(s[0]) + int(s[2]) + 2*int(s[3]) +   int(s[4]) + 2*int(s[5])
    return dZ, dN

def isBogus(Z, N, channel):
    """
    Check if disintegration channel leaves a valid residual nucleus.
    Note: We also mark channels with one residual proton or neutron as bogus
    because this occurs only for A < 12 and the corresponding cross sections
    from TALYS are clearly wrong.
    """
    dZ, dN = parseChannel(channel)
    if (dZ + dN) == 0:
        print '    no photo-disintegration:', channel
        return True
    if (dZ + dN) == (Z + N):
        print '    no nucleon left:', channel
        return True
    if dZ > Z:
        print '    too many protons lost:', channel
        return True
    if dN > N:
        print '    too many neutrons lost:', channel
        return True
    if ((N - dN) == 0) and ((Z - dZ) > 0):
        print '    residual nucleus contains only protons:', channel
        return True
    if ((Z - dZ) == 0) and ((N - dN) > 0):
        print '    residual nucleus contains only neutrons:', channel
        return True
    return False


fexcl = open('xs_pd.txt', 'w')
fexcl.write('# Z\tN\tchannel\txs\n')
fexcl.write('#cross sections [mb] for incident photon energies eps = 0.2 - 200 MeV in steps of logE = 0.01\n')
fmt = '%i\t%i\t%s' + '\t%.4g'*301 + '\n'  # output format

# parse all xs....tot files for each isotope
for (Z,N,A) in isotopes:
    print Z, N
    folder = path + '%i-%i/' % (Z, A)

    for f in os.listdir(folder):
        if not( f.startswith('xs') and f.endswith('tot') ):
            continue  # consider the xs....tot files only

        channel = f.strip('xs.tot')
        if isBogus(Z, N, channel):
            continue  # skip bogus channels

        xs = genfromtxt(folder + f, usecols=1)
        fexcl.write( fmt % ((Z, N, channel) + tuple(xs)) )

fexcl.close()


# -----------------------------------------------------------------------------
# Step 2:
# - Calculate the total (sum of exclusive) cross sections for each isotope --> xs_pd_sum.txt
# - Thin out channels: For each isotope select channels contributing at least 5% to the total cross section at any photon energy --> xs_pd_thin.txt
# -----------------------------------------------------------------------------
print ('Thinning exclusive channels')

data = genfromtxt('xs_pd.txt')
select = zeros(len(data), dtype=bool)

fsum = open('xs_pd_sum.txt', 'w')

for (Z, N, A) in isotopes:
    if Z <= 2:
        continue  # no TALYS output for H and He

    # find all channels for given isotope
    idx = (data[:,0] == Z) * (data[:,1] == N)
    xs = data[idx,3:]

    # save sum of exclusive cross-sections
    xs_sum = sum(xs, axis=0)
    fsum.write( ('%i\t%i' + '\t%g'*301 + '\n') % ((Z, N) + tuple(xs_sum)) )

    # calculate branching ratios, set to 0 when total cross-section = 0
    branching = nan_to_num(xs / xs_sum)

    # select channels with branching ratios above the threshold
    select[idx] = amax(branching, axis=1) > 0.05

fsum.close()

print ('Selected', sum(select), 'of' ,len(data), 'channels')
savetxt('xs_pd_thin.txt', data[select], fmt='%i\t%i\t%06d' + '\t%g'*301)
