Cosmic ray sources can be used to generate particles and feed them to a simulation.
The source is modular in design. You can add properties to a source, which define the particles that the source generates.

Source objects do not have a logic. The order in which the source properties modify a particle is the order they were added to the source. If several source properties modify the same particle property, previous modifications are overwritten.

Several sources can be combined in a SourceList, with individual (relative) luminosities.

* **SourceParticleType**
  * Type (particle ID= of the emitted particles
* **SourceMultipleParticleTypes**
  * Multiple emitted particle types with individual total abundances
* **SourceEnergy**
  * Discrete energy
* **SourcePowerLawSpectrum**
  * Energy drawn from a single power law spectrum
* **SourceComposition**
  * Particle type and energy for nuclei with a power law spectrum and individual abundances at equal energy per nucleon
* **SourcePosition**
  * Discrete initial position
* **SourceMultiplePositions**
  * Multiple discrete initial positions with inividual luminosities
* **SourceUniformDistribution**
  * Initial position drawn from a uniform distribution
* **SourceDensityGrid**
  * Initial position drawn from a 3D grid
* **SourceDensityGrid1D**
  * Initial position drawn from a 1D grid
* **SourceIsotropicEmision**
  * Isotropic initial direction
* **SourceDirection**
  * Discrete initial direction
* **SourceEmissionCone**
  * Initial direction drawn uniformly from within cone
* **SourceRedshift**
  * Discrete initial redshift
* **SourceUniformRedshift**
  * Initial redshift drawn from uniform distribution
* **SourceRedshift1D**
  * Initial redshift according to distance to observer