from numpy import *

# -------------------------------------------------
# Thin out photo-disintegration channels: For each isotope, select only those
# channels, whose cross section is at least 5% of the total cross section
# for any photon energy
# -------------------------------------------------

isotopes = genfromtxt('isotopes.txt', dtype=int)
d = genfromtxt('xs_excl.txt')
s = ones(len(d), dtype=bool)
print 'done reading'

fsum = open('xs_sum.txt', 'w')
fmt = '%i\t%i' + '\t%g'*301 + '\n'  # output format

for (Z, N, A) in isotopes:
    print Z, N
    # find all channels for given isotope
    idx = (d[:,0] == Z) * (d[:,1] == N)
    xs = d[:,3:][idx]

    # save sum of exclusive cross-sections
    xs_sum = sum(xs, axis=0)
    fsum.write( fmt % ((Z, N) + tuple(xs_sum)) )

    # calculate branching ratios, set to 0 when total cross-section = 0
    br = xs / xs_sum
    br[isnan(br)] = 0

    # find channels with branching ratios above the threshold
    s[idx] = amax(br, axis=1) > 0.05

fsum.close()
print 'selected', sum(s), 'of' ,len(d), 'channels'
savetxt('xs_thin.txt', d[s], fmt='%i\t%i\t%06d' + '\t%g'*301)
