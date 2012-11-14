from pylab import *
from scipy.interpolate import interp1d
from cStringIO import StringIO
import struct
from mpc import *


profile = StringIO("""
# Profile of magnetic field strength vs baryonic density in the Minati box
# Lower extrapolation B = 206.404 + 12.9872 * rho + 0.193117 * rho**2
# Higher extrapolation B = -37.3625 -2.74343 * rho -0.0585562 * rho**2
#
# log10(rho/[g/cm^3]), log10(B/[G])
-32.00000 -11.43459
-31.97993 -11.42361
-31.91850 -11.37153
-31.85708 -11.33202
-31.79566 -11.28935
-31.73424 -11.25074
-31.67281 -11.20560
-31.61139 -11.16113
-31.54997 -11.11273
-31.48854 -11.06385
-31.42712 -11.01348
-31.36570 -10.96107
-31.30427 -10.90346
-31.24285 -10.83643
-31.18143 -10.76371
-31.12001 -10.67964
-31.05858 -10.58822
-30.99716 -10.48548
-30.93574 -10.37223
-30.87431 -10.24660
-30.81289 -10.11632
-30.75147 -9.96395
-30.69004 -9.80515
-30.62862 -9.63804
-30.56720 -9.46599
-30.50578 -9.28175
-30.44435 -9.09687
-30.38293 -8.91053
-30.32151 -8.73668
-30.26008 -8.56165
-30.19866 -8.39986
-30.13724 -8.24704
-30.07581 -8.10704
-30.01439 -7.98537
-29.95297 -7.86703
-29.89155 -7.76893
-29.83012 -7.67947
-29.76870 -7.59368
-29.70728 -7.51493
-29.64585 -7.45123
-29.58443 -7.39008
-29.52301 -7.33190
-29.46158 -7.28033
-29.40016 -7.22940
-29.33874 -7.18511
-29.27732 -7.14020
-29.21589 -7.09864
-29.15447 -7.05942
-29.09305 -7.01471
-29.03162 -6.96864
-28.97020 -6.93470
-28.90878 -6.90071
-28.84735 -6.86472
-28.78593 -6.81636
-28.72451 -6.80306
-28.66309 -6.75334
-28.60166 -6.74123
-28.54024 -6.70163
-28.47882 -6.65852
-28.41739 -6.64615
-28.35597 -6.62567
-28.29455 -6.57452
-28.23312 -6.56336
-28.17170 -6.51866
-28.11028 -6.48613
-28.04886 -6.47825
-28.00000 -6.45452
""")

tablgRho, tablgB = genfromtxt(profile, unpack=True, comments='#')
f = interp1d(tablgRho, tablgB)

def modulation(x):
	if (x <= -32):
		return 206.404 + 12.9872 * x + 0.193117 * x * x
	elif (x >= -28):
		return -37.3625 - 2.74343 * x - 0.0585562 * x * x
	else:
		return f(x)

#gadget mass -db mhd_z.db -ox 54000 -oy 54000 -oz 54000 -size 132000 -bins 257 -o density_54-186Mpc_257bins.raw -dump
#gadget mass -db mhd_z.db -ox 54000 -oy 54000 -oz 54000 -size 132000 -bins 513 -o density_54-186Mpc_513bins.raw -dump

# takes a N+1 grid, the last entry is skipped to be conform with a periodic grid
fin = open('density_54-186Mpc_513bins.raw', 'rb')
fout = open('mod_MiniatiProfile_DolagDensity_54-186Mpc_512bins.raw', 'wb')
N = 513

sumb2 = 0
for ix in range(N):
	print ix
	for iy in range(N):
		for iz in range(N):
			rho = struct.unpack('f',fin.read(4))[0]
			rho *= 1.98892e40 * 0.7**2 / kpc**3 # density [kg/m^3]
			x = log10(rho) - 3 # log10(density [g/cm^3])
			b = 10**modulation(x) # field strength [G]
			b *= gauss # field strength [T]
			sumb2 += b**2
			fout.write(struct.pack('f', b))
		fin.read(4)
	fin.read(4 * (N + 1))

print 'Brms =', (sumb2 / N**3)**.5, "T"

fin.close()
fout.close()
