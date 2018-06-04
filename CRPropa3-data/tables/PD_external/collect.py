from numpy import *
import os

isotopes = genfromtxt('isotopes.txt', dtype=int)


# 1) collect all exclusive cross sections
f1 = open('xs_excl.txt', 'w')
f1.write('# Z\tN\tchannel\txs [mb]\n')
fmt1 = '%i\t%i\t%s' + '\t%.4g'*500 + '\n'  # output format

for (Z,N,A) in isotopes:
    print Z, N
    folder = '%i-%i/' % (Z, A)

    for f in os.listdir(folder):
        if not( f.startswith('xs') and f.endswith('tot') ):
            continue  # consider the xs....tot files only

        channel = f.strip('xs.tot')
        xs = genfromtxt(folder + f, usecols=1)
        f1.write( fmt1 % ((Z, N, channel) + tuple(xs)) )

f1.close()


# 2) sum up exclusive cross sections
d = genfromtxt('xs_excl.txt')

f2 = open('xs_sum.txt', 'w')
f2.write('# Z\tN\txs [mb]\n')
fmt2 = '%i\t%i' + '\t%g'*500 + '\n'

for (Z, N, A) in isotopes:
    print Z, N
    # find all channels for given isotope
    idx = (d[:,0] == Z) * (d[:,1] == N)
    xs = d[:,3:][idx]

    # save sum of exclusive cross-sections
    xs_sum = sum(xs, axis=0)
    f2.write( fmt2 % ((Z, N) + tuple(xs_sum)) )

f2.close()
