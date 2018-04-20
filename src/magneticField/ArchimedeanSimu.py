from crpropa import *

import time

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pylab import *

#######PARAMETER##########

N_max = 1000 # maximale Anzahl an Plot Punkten

energy = 10 * TeV     # Gyro
charge = 1

x = Vector3d(500 *kpc,0,10*pc)  # position
p = Vector3d(1.,0,0.0)  # direction

step = 2. * pc 
lenght = 100 * Mpc

m = energy / c_light/c_light
B_max = 4* muG # siehe Laborbuch
R_Gmin = m*c_light / eplus/B_max/pc

print "%s %e" % ("Anzahl Schritte = " , lenght/step)
print "%s %e" % ("minimaler Gyro-Radius [pc] = " , R_Gmin)
print "%s %e" % ("step size [pc] = " , step / pc )
print "%s %e" % ("steps/Gyro = " , R_Gmin /(step / pc))

##########################

Bfield = ArchimedeanSpiralField(1., 1., 1., 1.)

################MODUL########################
sim = ModuleList()

#~ sim.add(PropagationCK(Bfield, 1., step, step))
sim.add(PropagationBP(Bfield, step))

sim.add(MaximumTrajectoryLength(lenght))  
########OUTPUT###########

out = TextOutput('Halo_schlecht_BP.txt')

out.disableAll()
out.enable(Output.CurrentPositionColumn) 

#########################

#######OBSERVER############

obs = Observer()
obs.add(ObserverTimeEvolution(lenght/N_max, lenght/N_max, N_max) )  
obs.onDetection(out)
obs.setDeactivateOnDetection(False)
sim.add(obs) 

###########################

c = Candidate(nucleusId(1, charge), energy , x, p)

Zeit = time.time()

sim.run(c, True)
out.close()

################GRAPHIK########################

print "%s %e" % ("Rechenzeit = " , time.time() - Zeit)
print "%s %e" % ("Rechenzeit pro Schritt = " , (time.time() - Zeit)/(lenght/step))


data = pd.read_csv("Halo_schlecht_BP.txt", names = ['X', 'Y', 'Z'], delimiter='\t', comment='#')

x, y, z = np.array(data['X']), np.array(data['Y']), np.array(data['Z']*1000000.)

fig = plt.figure()
#~ ax = fig.add_subplot(111) # 2D
#~ ax.scatter(x, y)
ax = fig.gca(projection='3d') #3D
ax.scatter(x[::1], y[::1], z[::1])

ax.set_xlabel('x [Mpc]')
ax.set_ylabel('y [Mpc]')
ax.set_zlabel('z [pc]')

plt.show()
