## CRPropa vNEXT

### Bug fixes:

### New features:

* Planck JF12b variant of the JF12Field. See arXiv:1601.00546. Thanks to
	Mikhail Zotov for contributing.
* ParticleCollector can provide Candidates directly to ModuleList::run
* Basic file versioning of the data archive in CMakeLists.txt

### Interface change:

* ParticleCollector::getAll() -> ParticleCollector::getContainer()

### Features that are deprecated and will be removed after this release:

### New plugins and resources linked on the webpages:

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
