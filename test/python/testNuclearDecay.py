# CRPRopa test script
# Visualizes the nuclear decay tables in CRPropa
#
import matplotlib
matplotlib.use('Agg')
import numpy as np

from crpropa import *
from pylab import *
from matplotlib.colors import LogNorm


# color values of decay modes
vS = 0.5 # stable
vP = 0.9 # proton emission
vN = 0.125 # neutron emission
vA = 0.99 # alpha decay
vBp = 0.25 # beta+ decay
vBm = 0.75 # beta- decay
modeDict = {1:vN, 2:vN, 10:vP, 20:vP, 100:vA, 200:vA, 10000:vBp, 20000:vBp, 10001:vBp, 10002:vBp, 10003:vBp, 10004:vBp, 10101:vBp, 10100:vBp, 10200:vBp, 10300:vBp, 1000:vBm, 2000:vBm, 1100:vBm, 1010:vBm, 1020:vBm, 1030:vBm}


def readtable(pathtofile, comments):
    data = []
    with open(pathtofile, 'r') as f:
        for line in f:
            if line.startswith(comments):
                continue
            items = line.split()
            converted_line = list(map(int, items[0:3]))
            converted_line += list(map(float, items[3:]))
            data.append(converted_line)
    return data

### read and evaluate decay tables
data = readtable(getDataPath('nuclear_decay.txt'), comments='#')
modes = zeros((27, 31))
times = zeros((27, 31))
multi = zeros((27, 31))
inclu = zeros((27, 31))

for el in data:
    z = el[0]
    n = el[1]
    mode = el[2]
    time = el[3]
    inclu[z,n] += 1/time
    multi[z,n] += 1
    if modes[z, n] != 0:
        if time >= times[z, n]:
            continue
    modes[z,n] = modeDict[mode]
    times[z,n] = time

### stable nuclei
inclu[np.where(inclu == 0)] = 1e-99

inclu = 1 / inclu

modes[1,0] = vS
for z in range(1,27):
    for n in range(1,31):
        if modes[z, n] == 0:
            modes[z, n] = vS


### plot data
fig = figure()
ax = fig.add_subplot(111)
ax.imshow(ma.masked_array(modes, modes==0), aspect='equal', interpolation='nearest', origin='lower')
ax.set_xlabel('Neutrons')
ax.set_ylabel('Protons')
ax.grid()
ax.text(0.41, 0.59, 'beta+ decay', ha='center', va='center', rotation=40, transform=ax.transAxes)
ax.text(0.48, 0.52, 'stable', ha='center', va='center', rotation=40, transform=ax.transAxes)
ax.text(0.56, 0.44, 'beta- decay', ha='center', va='center', rotation=40, transform=ax.transAxes)
ax.text(0.78, 0.22, 'neutron dripping', ha='center', va='center', rotation=40, transform=ax.transAxes)
ax.text(0.2, 0.8, 'proton dripping', ha='center', va='center', rotation=40, transform=ax.transAxes)
fig.savefig('NuclearDecay_mode.png',bbox_inches='tight')

fig = figure()
ax = fig.add_subplot(111)
im = ax.imshow(ma.masked_array(inclu, inclu==np.inf), aspect='equal', interpolation='nearest', origin='lower', norm=LogNorm(vmin=1e-3, vmax=1e3))
cbar = fig.colorbar(im, orientation='horizontal')
cbar.set_ticks(logspace(-3,3,7))
cbar.set_ticklabels(['0.001','0.01','0.1','1','10','100','1000'])
cbar.set_label('Lifetime [s]')
ax.set_xlabel('Neutrons')
ax.set_ylabel('Protons')
ax.grid()
fig.savefig('NuclearDecay_lifetime.png',bbox_inches='tight')

fig = figure()
ax = fig.add_subplot(111)
mmulti = ma.masked_array(multi, multi==0)
im = ax.imshow(mmulti, aspect='equal', interpolation='nearest', origin='lower')
cbar = fig.colorbar(im, orientation='horizontal')
cbar.set_label('\# Decay Channels')
ax.set_xlabel('Neutrons')
ax.set_ylabel('Protons')
ax.grid()
fig.savefig('NuclearDecay_multiplicity.png',bbox_inches='tight')

show()
