This list gives an overview over the simulation modules together with a short description.
For a more detailed explanation refer to the Doxygen documentation.

### Propagation modules
Propagation modules are responsible for proposing a step size, evaluating the bids for the step size of the previous round and spatially moving the particle according to this step. Every simulation needs exactly one propagation module, that is usually put at the beginning of the module list.

* **SimplePropagation** 
  * Rectalinear propagation, i.e. propagates a candidate from a given starting position (e.g. source) with an initial direction onwards until one of the breaking conditions is fulfilled
* **PropagationCK** 
  * Deflections of charged particles in magnetic fields using the Cash-Karp algorithm (Runge-Kutta of order 4/5)

### Interaction modules
Interaction modules implement physical interactions which modify the particle and possibly produce secondaries

* **ElectronPairProduction**
  * Electron pair production energy losses for charged nuclei using the continuous energy loss approximation
* **PhotoPionProduction**
  * Photo meson production for nuclei using free nucleon approximation
  * Uses SOPHIA to calculate the outcome of a photopion interaction
* **PhotoDisintegration**
  * Photo disintegration using TALYS tables
* **NuclearDecay**
  * Nuclear decay
  * Electron and neutrino secondaries
* **Redshift**
  * Redshift calculation
  * Adiabatic energy loss

### Boundary modules
Boundaray modules flag a particle once a certain condition is fulfilled. Optionally they can signal the termination of propagation

#### Break conditions
* **MaximumTrajectoryLength**
* **MinimumRedshift**
* **MinimumEnergy**

#### Boundaries
Cubic-, Spherical- and EllipsoidalBoundary flag (and optionally inactivate) particles that exit the volume they define. They can be used to limit the simulation volume.
* **CubicBoundary**
* **SphericalBoundary**
* **EllipsoidalBoundary**
Periodic- and ReflectiveBox implement boundary conditions for the particles. They are useful for a 3D setup where an initial volume is to be repeated (periodically or reflectively). When using them, a couple of things need to be considered. Observers will shadow the volume behind if they are set to inactivate particles. Also Observers should be placed at a distance to the boundaries that is larger than the maximum step size of the propagator, since step size limitation does not work beyond periodic/reflective boundaries.
* **PeriodicBox**
  * Periodic boundary conditions for the particle
  * If a particle leaves the box it will enter from the opposite side and the initial position will be changed as if it had come from that side.
* **ReflectiveBox**
  * Reflective boundaray conditions for the particle
  * If a particle leaves the box it will be reflected (mirrored) and the initial position will be changed as if it had come from that side.

#### Observers
Observers can be defined using a collection of ObserverFeatures.
The names of ObserverFeatures all start with "Observer" so you can discover the available options from an interactive python sesion by typing "Observer" and pressing "tab". The list includes
* **ObserverSmallSphere**
  * Detects particle when they enter the sphere
* **ObserverLargeSphere**
  * Detects particles when they leave the sphere
* **ObserverTracking**
  * For recording the tracks of particles inside a small observer sphere
* **ObserverPoint**
  * Observer for 1D simulations
* **ObserverDetectAll**
  * Detects all particles
* **ObserverRedshiftWindow**
  * Detect particles within a given redshift interval around z=0
* **ObserverInactiveVeto**
  * Prevent inactive particles from being detected
* **ObserverPhotonVeto**, **ObserverElectronVeto**, **ObserverNeutrinoVeto**, **ObserverNucleusVeto**

### Output modules
Main output modules
* **ShellOutput** Output to the shell
* **TextOutput** Plain text output
  * Customizable with the presets Event1D, Event3D, Trajectory1D, Trajectory3D, Everything.
  * If the filename ends with '.gz' the output is compressed 
* **HDF5Output** Output in the HDF5 format

Legacy output modules (CRPropa 2 format)
* **ROOTEventOutput1D**
* **ROOTEventOutput3D**
* **ROOTTrajectoryOutput1D**
* **ROOTTrajectoryOutput3D**
* **CRPropa2EventOutput1D**
* **CRPropa2EventOutput3D**
* **CRPropa2TrajectoryOutput1D**
* **CRPropa2TrajectoryOutput3D**

### Other modules
* **PerformanceModule**
  * Measure execution time for a number of modules