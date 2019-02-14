#!/usr/bin/env python
# coding: utf-8

# # Comparison of Propagation Modules (BP - CK)

# 
# 
# Numerical simulations of the propagation of charged particles through magnetic fields, solving the equation of motion can be achieved in principle with many different algorithms. There are, however, an increasing number of studies that have found that there are two algorithms, which work in general best for propagating charged particles within a magnetic field. These two algorithms are compared and evaluated.
# 
# Both the Boris push and the Cash-Karp algorithms solve the Lorentz equation and thus enable the propagation of charged particles through magnetic fields. For the simplest case of a homogeneous background field, the particle trajectory can be derived analytically. This enables comparison with numerical integrators and provides information on the errors of the algorithms for the corresponding parameters being used. The aim is to find a relationship between the error and the parameters applied. The pitch angle and the step length are suitable parameters to determine the influence. 
# 
# Consequently, we want to understand which algorithm is suitable for each simulation setup.

# ## Analytic Solution

# The trajectory of a charged particle in a homogeneous background field can be solved analytically. We assume that the particle starts off with non-zero component of the momentum ($p_z > 0$) along the background magnetic field. The particle, then moves in a helix whose gyration radius is given by: $r_g = \frac{E}{B\cdot q \cdot c}$
# 
# The velocity component parallel to the background field remains constant. The analytical solution for the $direction = Vector3d(p_x, p_y, p_z)$ with $p_x^2+p_y^2+p_z^2 = 1$ and the $position = Vector3d(0, 0, 0)$ yields:

# In[1]:


def analytical_solution(max_trajectory, p_z, r_g_0, number_steps):
    # calculate the time stamps similar to that used in the numerical simulation 
    t = np.linspace(0, max_trajectory/pc, number_steps+1)
    # shift the phase so that the analytical solution 
    # also starts at (0,0,0) with in the direction (p_x,p_y,p_z)
    d = t[1:]/r_g_0-3*math.pi/4.
    
    # for parallel motion corrected gyro radius
    r_g = r_g_0*(1-p_z**2)**0.5 

    # at these trajectory lengths, the numerical solutions are known
    x_ana = r_g*np.cos(d)
    y_ana = -r_g*np.sin(d)
    z_ana = p_z*t[1:]
    return x_ana, y_ana, z_ana


# ## Helper Functions

# We define functions with which we can calculate the gyration radius and the circumference. In addition we only want to define how many steps we require per circumference and how many gyrations should be performed. The last helping function should then calculate the maximum trajectory length and the step size, so that our requirements are fulfilled.

# In[3]:


import numpy as np
from crpropa import *
import math

# Gyro radius should be R_g = 10.810076 parsecs for p_z = 0 and B = 10nG and E = 100 TeV
def larmor_radius(c, field):
    p = c.current.getMomentum()
    try:
        B = field.getRegularField(c.current.getPosition())
    except:
        B = field.getField(c.current.getPosition())
        q = c.current.getCharge()
        p_perp = p.getPerpendicularTo(B)
    try:
        r =  abs(p_perp.getR() / (B.getR()*q))
    except ZeroDivisionError:
        r = 1000000
    return r

# Calculate gyro radius
def larmor_circumference(c, field):
    r_g_0 = larmor_radius(c, field)/pc # gyro radius for p_z = 0
    l_gyration = 2*math.pi*r_g_0*pc    # trajectory lenght of one gyration
    return r_g_0, l_gyration

# Trajectory length of particle in number of gyrations
def maximum_trajectory(steps_per_gyrations, number_of_gyrations, c, field, p_z):
    r_g_0, l_gyration = larmor_circumference(c, field)
    max_trajectory = l_gyration*number_of_gyrations
    steplength = l_gyration/np.array(steps_per_gyrations)
    return max_trajectory, steplength, r_g_0


# ## Run Simulation 

# We can now compare the analytical solution with both particle integrators in CRPropa, namely the Boris push and the Cash-Karp.
# 
# First we have to add our propagation module with the above specialized background magnetic field to our module list. Afterwards we can specify that the particle information is collected at each step along the particle trajectory. With the output module, we can specify where these information will be saved. 
# 
# Now that we have our module list ready we can fire up our simulation with both propagation algorithms and hope that something visually interesting is going to happen. 

# In[4]:


import time as Time

# We use only a Background magnetic field in the z-direction. 
# We could add more complex magentic fields to our MagneticFieldList.
B = 10*nG
direction_B = Vector3d(0, 0, 1)
  
const_mag_vec = direction_B * B
reg_field = UniformMagneticField(const_mag_vec)


### Running the simulation with either CK or BP
def run_simulation(module, steps_per_gyrations, number_of_gyrations, p_z):
    # Initial condition of candidate
    p_x = (1-p_z**2)**(1/2.)/2**0.5
    p_y = p_x
    E = 100 * TeV
    direction = Vector3d(p_x, p_y, p_z)
    position = Vector3d(0, 0, 0)
    
    c = Candidate(nucleusId(1, 1), E, position, direction)
    
    max_trajectory, steplength, r_g_0 = maximum_trajectory(steps_per_gyrations, number_of_gyrations, c, reg_field, p_z)
    
    sim = ModuleList()
    if module == 'CK':
        sim.add(PropagationCK(reg_field,1e-4,steplength, steplength))
        output = TextOutput('trajectory_CK.txt', Output.Trajectory3D)
    elif module == 'BP':
        sim.add(PropagationBP(reg_field, steplength))
        output = TextOutput('trajectory_BP.txt', Output.Trajectory3D)
    else:
        print('no module found. Use either BP or CK.')
        return
    
    # we only want to simulate a certain trajectory length
    sim.add(MaximumTrajectoryLength(max_trajectory))
    # the output information will be saved in pc instead of the default which is Mpc
    output.setLengthScale(pc)
    # each particle position will be saved in the above specified text field.
    sim.add(output)
    # compare the simulation time of both propagation methods
    t0 = Time.time()
    # run the simulation
    sim.run(c, True)
    t1 = Time.time()
    output.close()
    print('Simulation time with module '+ str(module)+' is '+str(t1-t0)+'s.')
    Time.sleep(0.5)
    return max_trajectory, p_z, r_g_0


# ## Load Simulation Data

# There are several ways to load the simulation data. Pandas is helpful to load files. To illustrate this, we will load and process the data with pandas.

# In[5]:


import pandas as pd

def load_data(text, r_g):
    data = pd.read_csv(text, 
                 names=['D','ID','E','X','Y','Z','Px','Py','Pz'], delimiter='\t', comment='#',
                 usecols=["D", "X", "Y", "Z","Px","Py","Pz"])

    ### distances are saved in units of pc
    ### transform so that the center of the gyromotion is at (0,0)
    data.X = data.X.values-r_g/2**0.5
    data.Y = data.Y.values+r_g/2**0.5
    ### convert disctance in kpc
    data.D = data.D.values/1000.
    
    ### calcualte gyro radius
    data['R'] = (data.X**2+data.Y**2)**0.5

    return data


# ## Compare the Propagation in the Perpendicular Plane

# Here, we can compare the numerical results of both modules with the analytical solution of the particle trajectory. First, we want to show the trajectories in the $xy$-plane (left plot). The black squares present the positions of the analytical solution for the chosen step size. In the middle plot, we can display the relative error of the gyro radius as a function of time. Finally, we can compute the distance of the numerical solutions to the analytical solution in the $xy$-plane (right plot).

# In[6]:


import matplotlib.pyplot as plt

steps_per_gyrations = 10
number_gyrations = 40000
number_of_steps = steps_per_gyrations * number_gyrations

def plot_subplots(ax1, ax2, ax3, data, x_ana, y_ana, r_g, module, color):
    # numerical calculated positions
    
    ax1.scatter(data.X,data.Y, s=1,color = color, label = module)
    
    # analytical solution shwon in black squares!
    if module == 'CK':
        ax1.scatter(x_ana, y_ana, color = 'k', marker = 's', s = 15)
        # for the legend
        ax1.scatter(x_ana, y_ana, color = 'k', marker = 's', s = 1, label = 'Analytical solution')
    
    ax1.legend(markerscale=5)
    
    # numerical solutions
    ax2.scatter(data.D,(data.R-r_g)/r_g*100., s=1,color = color, label= module)
    ax2.legend(markerscale=5)
    
    ax3.scatter(data.D,((x_ana-data.X)**2+(y_ana-data.Y)**2)**0.5, s=1,color = color, label=module)
    ax3.legend(markerscale=5)
   
 # We use this function to plot the whole figure for the particle motion in the xy-plane
def plot_figure_perp(max_trajectory, p_z, r_g_0, number_of_steps):
    fig, ((ax1, ax2, ax3)) = plt.subplots(1, 3,figsize=(12,4))

    # for parallel motion corrected gyro radius
    r_g = r_g_0*(1-p_z**2)**0.5 
    
    # Initial condition of candidate
    p_x = (1-p_z**2)**(1/2.)/2**0.5
    p_y = p_x
    
    x_ana, y_ana, z_ana = analytical_solution(max_trajectory, p_z, r_g_0, number_of_steps)
    
    data = load_data('trajectory_BP.txt', r_g)
    plot_subplots(ax1, ax2, ax3, data, x_ana, y_ana, r_g, 'BP', 'brown')

    data = load_data('trajectory_CK.txt', r_g)
    plot_subplots(ax1, ax2, ax3, data, x_ana, y_ana, r_g, 'CK', 'dodgerblue')

    ax1.set_xlabel('$x$ [pc]')
    ax1.set_ylabel('$y$ [pc]')
    ax1.set_title('$p_z/p$ = '+str(p_z), fontsize=18)

    ax2.set_xlabel('distance [kpc]')
    ax2.set_ylabel('relative error in $r_\mathrm{g}$ [%]')
    ax2.set_title('$p_z/p$ = '+str(p_z), fontsize=18)

    ax3.set_xlabel('distance [kpc]')
    ax3.set_ylabel('$xy-$deviation from ana. solution [pc]')
    ax3.set_title('$p_z/p$ = '+str(p_z), fontsize=18)
    
    fig.tight_layout()
    plt.show()


# We can study $\textbf{different pitch angles}$ of the particle with respect to the background magnetic field and start with $p_z/p = 0.01$ which represents a strong perpendicular component:

# In[7]:


p_z = 0.01
max_trajectory, p_z, r_g_0 = run_simulation('CK', steps_per_gyrations, number_gyrations/2, p_z)
run_simulation('BP', steps_per_gyrations, number_gyrations/2, p_z)
plot_figure_perp(max_trajectory, p_z, r_g_0, number_of_steps/2)


# We increase the component of the momentum along the parallel component:

# In[8]:


p_z = 0.5
max_trajectory, p_z, r_g_0 = run_simulation('CK', steps_per_gyrations, number_gyrations, p_z)
run_simulation('BP', steps_per_gyrations, number_gyrations, p_z)
plot_figure_perp(max_trajectory, p_z, r_g_0, number_of_steps)


# In[9]:


p_z = 0.99
max_trajectory, p_z, r_g_0 = run_simulation('CK', steps_per_gyrations, number_gyrations, p_z)
run_simulation('BP', steps_per_gyrations, number_gyrations, p_z)
plot_figure_perp(max_trajectory, p_z, r_g_0, number_of_steps)


# In[10]:


p_z = 0.9999
max_trajectory, p_z, r_g_0 = run_simulation('CK', steps_per_gyrations, number_gyrations, p_z)
run_simulation('BP', steps_per_gyrations, number_gyrations, p_z)
plot_figure_perp(max_trajectory, p_z, r_g_0, number_of_steps)


# $\textbf{Conclusions:}$ 
# - the $\textbf{Boris push}$ propagates the particle in the perpendicular plane, as expected, on a cirlce with a radius that is in great agreement with the gyro radius. The problem is that the frequency deviates from the analytical solution. Consequently, even though we can see an agreement between the numerical and analytical gyro radius, the deviation in the $xy$-plane is between 0 and $2r_g$. The $\textbf{global error is constrained}$.
# - the $\textbf{Cash-Karp}$ algorithm results in a relative $\textbf{small local error}$ in the perpendicular plane for small components of the momentum along the parallel component of the magnetic field. The problem is that the small $\textbf{local error accumulates to a large global error}$. In addition, the $\textbf{error depends on the pitch angle}$ (angle between particle direction and background magnetic field) of the particle and thus may artificially pollute the results if we use a background field in combination with an isotropic source of particles.
# 
# + Note that other implementations of the CashKarp algorithm lead to different results, since CRPropa prevents the candidate from losing energy due to propagation alone. This is not true in general for other implementations. 
