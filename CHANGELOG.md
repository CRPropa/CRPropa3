# CRPropa NEXT

### Bug fixes:


### New features:
* new candidate property tagOrigin to trace back which source or which interaction created the candidate

### Interface changes:


### Features that are deprecated and will be removed after this release


### New plugins and resources linked on the webpages:



## CRPropa 3.2

### Bug fixes:
* Fix of reflective boundary condition for scalar- and vectorgrids
  that showed asymmetry and discontinuities (See issue [#361]).
* Fix in EMTripletPairProduction
* Fix of the data files of the Hackstein EGMF models as well as the 
  corresponding example notebook.
* Fix of axis normalization of getRotated in Vector3.h.
* Fix of secondary spectra in electromagnetic interactions
  (EM*-modules), issue [#334] and pull request [#15] in crpropa-data.
* Fix weight inheritance for secondary particles; they are created with 
  their parents weights as intial weights now.

### New features:
* Add modules for first and second order Fermi acceleration
* Make custom photon fields available in the PhotoPionProduction module.
* Add tricubic- and nearest neighbour interpolation routines for scalar-
  and vectorgrids.
* Add the new PolarizedSingleModeMagneticField class for polarized/
  helical single mode magnetic field models.
* Add a source feature for targeted emission, following the 
  von-Mises-Fisher distribution
* Updates in SNR and pulsar source distributions 

### Interface changes:
* Plane wave and grid turbulence models use same parameter convention now 

### Features that are deprecated and will be removed after this release
* External extensions DINT and Eleca, which can be replaced with the 
  EM*-modules combined with the thinning option for reasonable computation
  times.

### New plugins and resources linked on the webpages:
* Updated version of the CLUES EGMF models


## CRPropa 3.1.7

### Bug fixes:

* Re-add URB that was temporarily not available in the (unreleased) master branch.
* Turbulent fields generated on a grid were limited up to 2048 grid-size due to
  an integer overflow (i.e. 2048^3 index > signed int); solved by replacing
  int with size_t

### New features:

* New model for the radio background added: Nitu et al. 2021 measurements.
* New model for the radio background added: Fixsen et al. 2011 (ARCADE-2) measurements.
* Weighted sampling thinning of electromagnetic processes  (EMPairProduction,
  EMInverseComptonScattering, EMDoublePairProduction, EMTripletPairProduction).
* Planck JF12b variant of the JF12Field. See arXiv:1601.00546. Thanks to
	Mikhail Zotov for contributing.
* ParticleCollector can provide Candidates directly to ModuleList::run
* Basic file versioning of the data archive in CMakeLists.txt
* Python docstrings are generated automatically from doxygen documentation
* Photon field classes replaced the photon field enumerator,
  consequently, the new photon fields implementation follows the same logic as
  of the other modules making it is easier to introduce custom ones
* New class-based interface for turbulent fields introduced
* New turbulence modules implemented:
  - GridTurbulence (with the bendover scale) which should in general be used
    instead of initTurbulence (before, it was implemented as
    initTurbulenceWithBendover);
  - SimpleGridTurbulence which provides the exact field as initTurbulence;
  - HelicalGridTurbulence which provides the exact field as initHelicalTurbulence
  - PlaneWaveTurbulence a new algorithm based on Giacalone & Jokipii, 1999 and
    Tautz & Dosch, 2013.
* New CMake option: `BUILD_DOC` for building Doxygen & Sphinx docs


### Interface change:

* ParticleCollector::getAll() -> ParticleCollector::getContainer()
* Photon fields are no longer items of the PhotonField enumerator but independent
  classes that share the same interface, so instead of `CMB` one should use `CMB()`,
  instead of `IRB_Kneiske04` - `IRB_Kneiske04()`, etc.
* initTurbulenceWithBendover() removed (as it was just briefly present in the code)
  and replaced with GridTurbulence

### Features that are deprecated and will be removed after this release:

* Turbulence-related functions: initTurbulence, turbulentCorrelationLength,
  initHelicalTurbulence


## CRPropa v3.1.6

### Bug fixes:

* Fix of (#254): Redshift evolution in PhotoPionProduction
  The reshift evolution was always handled with simple scaling and never with
  the more accurate 2 dimensional interpolation.

### New features:

* New source feature SourceLambertDistributionOnSphere
  Simplifies simulations of scenarios with an initially isotropic and homogeneous
  distribution of cosmic rays on a sphere, e.g. for investigation of propagation
  of extagalactic cosmic rays in the Milky Way  (see issue #246 and pull
  request #247)
* For PhotoPionProduction now a two dimensional interpolation of the redshift
  evolution is available (see pull request  #255)
* A new break condition MinimumChargeNumber is added (see pull request #256)
* Vector3 now has index based access to its components. This improves
  interoperability with third party libraries, in particular with numpy arrays.
  (see pull request #262)
* Random seeds can be accessed from python (see pull request #263)
* New galactic magnetic field model by Terral & Ferriere (2017) (see pull request #258)
* Reimplementation of SOPHIA's photon field sampling used in
PhotoPionProduction in c++, leading to a factor 2-3 speed up of the
module (see pull request #260).
* Introducing a method, sophiaEvent(onProton, Eprimary, Ephoton), to
directly call SOPHIA's event generator from python (see pull request #260).


## CRPropa v3.1.5

### Bug fixes:
 * Fixed issue with secondaries in the PhotoPionProduction potentially relevant
   for primary energies above approx. 10**21 eV (See issue [#225]).
 * Fix of azimuthal component of ArchimedeanSpiral Field (See commit 1f79e2c).
   Thanks to Leander Schlegel for reporting and providing a patch.

### New features:
 * The PropagatorBP implements the Boris-Push algorithm as as alternative to
   the Cash-Karp integrator PropagatorCK. Thanks to Patrick Reichherzer for
   this contribution.
 * A basic geometry system was added to e.g. restrict modules to a spatial region.
 * Candidates now have weights.
 * User can define custom candidate properties identified with a string (slow but sometimes handy).
 * Random seeds can be stored in output files. Please note that perfect
   reproducibility of simulations is still not guaranteed when using multi
   cores.
 * The 'ParticleCollector' allows in-memory storage of Candidates for
   collective-processing during a simulation.
 * The SDE solver supports advection and adiabatic cooling.
 * A selection of models for mass distributions in our Galaxy as first step of
   the development of a hadron-hadron interaction module was added.
   Thanks to Julien Dörner for implementing these models.
 * Source distributions following the Galactic pulsar distribution and the
   Galactic SNR distributions are available.
 * Turbulent magnetic fields can now be optionally helical
   (initHelicalTurbulence) or obey the homogeneous turbulence theory
   requirement (initTurbulenceWithBendover).
 * Modified JF12 Galactic magnetic field model with a truly solenoidal disk
   spiral field and a smooth X field with parabolic field lines according to
   Kleimann et al. (arXiv:1809.07528). Thanks to Timo Schorlepp for sharing
   this code.

### Features that are deprecated and will be removed after this release:
 * PhotonOutput1D is obsolete as the functionality is provided by the regular
   output modules (See #211).
 * Future releases will use C++11 features. Outdated compilers will not be
   supported anymore.

### New plugins and resources linked on the webpages:
 * ROOTOutput. This code has been moved from CRPropa to a plugin and is no
   longer maintained by us. If you want to maintain this code, please contact
   us.
 * Data for the constrained 'Hackstein' models of the local Universe using
   initial conditions from the CLUES project. Thanks to Stefan Hackstein for
   sharing.
