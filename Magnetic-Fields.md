* **MagneticFieldList**
  * Container for the superposition of several fields
* **UniformMagneticField**
  * Uniform magnetic field
* **MagneticFieldGrid**
  * Regular grid of with trilinear interpolation
  * The grid is either periodically or reflectively repeated
* **ModulatedMagneticFieldGrid**
  * Vector grid scaled with a scalar grid
  * The grids can have different sizes
* **SPHMagneticField**
  * Magnetic field from large scale structure simulations using the SPH formalism
* **SPHMagneticFieldGrid**
  * Magnetic field from large scale structure simulations using the SPH formalism, precalculated on a regular grid
* **TurbulentMagneticField**
  * Turbulent field from a superposition of a number of random linear waves
  * Implementation currently not correct! For a turbulent use MagneticFieldGrid and use initTurbulence
* **JF2012Field**
  * Galactic magnetic field from the JF2012 model
  * Regular, striated and turbulent
