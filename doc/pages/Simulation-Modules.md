## Modules

This list gives an overview over the simulation modules together with a short description.
For a more detailed explanation refer to the Doxygen documentation.

### Propagation modules
Propagation modules are responsible for proposing a step size, evaluating the bids for the step size of the previous round and spatially moving the particle according to this step. Every simulation needs exactly one propagation module, that is usually put at the beginning of the module list.

* **SimplePropagation** - Simple rectlinear propagation
* **PropagationBP** - Deflections of charged particles in magnetic fields using the Boris Push algorithm with dynamic step size control
* **PropagationCK** - Deflections of charged particles in magnetic fields using the Cash-Karp algorithm (Runge-Kutta of order 4/5) with dynamic step size control
* **DiffusionSDE** - Solves the Fokker-Planck transport equation using stochastic differential equations (SDEs).

### Interaction modules
Interaction modules implement physical interactions which modify the particle and eventually produce secondary particles. Hadronic secondaries are always generated, non-hadronic secondaries are optionally generated.
Currently, only interactions with extragalactic background photon fields (CMB, EBL, CRB) are implemented.
Hadronic interactions with matter distributions is highly subdominant except for high density regions, and is currently not implemented.

Interactions of protons, neutrons, and nuclei (Z = 1 - 26, N = 1 - 30)

* **ElectronPairProduction** - Electron pair production (Bethe-Heitler) for charged nuclei using the continuous energy loss approximation, optional secondaries: electrons/positrons
* **PhotoPionProduction** - photo-meson production for protons, neutrinos and nuclei, uses SOPHIA as event generator, secondaries: protons/neutrons, optional secondaries: antiprotons/antineutrons, photons, electrons/positrons and neutrinos
* **PhotoDisintegration** - photodisintegration using TALYS cross sections (alternatively, PSB and Kossov models are available), secondaries: protons, neutrons, deuterons, tritons, alpha-3, alpha-4, optional secondaries: photons
* **NuclearDecay** - decay of neutrons and nuclei up to iron, optional secondaries: photons, electrons/positrons and neutrinos

Interactions of photons, electrons and positrons

* **EMPairProduction** - electron pair production (Breit-Wheeler process), optional secondaries: electrons/positrons
* **EMDoublePairProduction** - double electron pair production, optional secondaries: electrons/positrons
* **EMTripletPairProduction** - triplet pair production, optional secondaries: electrons/positrons
* **EMInverseComptonScattering** - inverse compton scattering, optional secondaries: photons

General interactions/processes

* **Redshift** - updates the redshift and calculates the adiabatic energy loss
* **FutureRedshift** - same as Redshift, but allows for negative redshifts (for symmetric window around observer)
* **SynchrotronRadiation** - synchrotron radiation of charged particles in magnetic fields, optional secondaries: photons
* **AdiabaticCooling** - takes adiabatic cooling (or heating) of the particles due to expansion (or compression) of the plasma into account

### Conditional modules
Conditional modules implement certain conditions for stopping propagation.
They provide interfaces to act onReject or onAccept of a cosmic ray.
Boundary modules can be used to limit the simulation volume.
Periodic-/ReflectiveBox and ReflectiveShell implement boundary conditions for the particles. They are useful for a 3D setup where an initial volume is to be repeated (periodically or reflectively). When using them, a couple of things need to be considered. Observers will shadow the volume behind if they are set to inactivate particles. Also Observers should be placed at a distance to the boundaries that is larger than the maximum step size of the propagator, since step size limitation does not work beyond periodic/reflective boundaries.

* **MaximumTrajectoryLength** - Stop after reaching maximum trajectory length
* **MinimumEnergy** - Stop after reaching a minimum energy
* **MinimumEnergyPerParticleId** - Same as `MinimumEnergy` but allows for different particle types to have individual minimum energies.
* **MinimumRedshift** - Stop after reaching a minimum redshift
* **CubicBoundary** - Cubic simulation volume
* **SphericalBoundary** - Spherical simulation volume
* **EllipsoidalBoundary** - Ellipsoidal simulation volume
* **CylindricalBoundary** - Cylindric simulation volume
* **PeriodicBox** - Periodic boundary conditions for the particle: If a particle leaves the box it will enter from the opposite side and the initial position will be changed as if it had come from that side.
* **ReflectiveBox** - Reflective boundary conditions for the particle: If a particle leaves the box it will be reflected (mirrored) and the initial position will be changed as if it had come from that side.
* **ReflectiveShell** - Reflective spherical boundary condition for the particle; If a particle crosses the shell during the current step it will be reflected back at the boundary according to the law of reflection. Particles coming from the outside are only reflected if the current step ends inside the shell.
* **DetectionLength** - Detects the candidate at a given trajectory length.

### Observers
Observers can be defined using a collection of ObserverFeatures.
The names of ObserverFeatures all start with "Observer" so you can discover the available options from an interactive python session by typing "Observer" and pressing "tab". The list includes
* **ObserverSurface** - Detects particles crossing the boundaries of a defined surface (see, e.g., `Geometry` module)
* **ObserverTracking** - For recording the tracks of particles inside an observer sphere
* **Observer1D** - Observer for 1D simulations that detects particles when reaching x = 0
* **ObserverDetectAll** - Detects all particles
* **ObserverRedshiftWindow** - Detect particles within a given redshift interval
* **ObserverInactiveVeto** - Veto for inactive particles
* **ObserverPhotonVeto** - Veto for photons
* **ObserverElectronVeto** - Veto for electrons/positrons
* **ObserverNeutrinoVeto** - Veto for neutrinos
* **ObserverNucleusVeto** - Veto for protons/neutrons and nuclei
* **ObserverTimeEvolution** - Records all candidates along their trajectory using linear or logarithmic steps

### Output modules
Main output modules
* **ShellOutput** - Output to the shell
* **TextOutput** - Plain text output, customizable with the presets Event1D, Event3D, Trajectory1D, Trajectory3D, Everything, or more fine grained control. If the filename ends with '.gz' the output is compressed.
* **HDF5Output** - Output in the HDF5 format
* **ParticleCollector** - A temporary container for storing candidates in memory (use with care due to memory limitations, e.g. 1e6 candidates ~ 500MB of RAM)

### Other modules
* **PerformanceModule** - Measure execution time for a number of modules
