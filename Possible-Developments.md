* **Photodisintegration cross-sections** The photodisintegration cross-sections currently pose the largest uncertainty for the propagation of cosmic ray nuclei. A dedicated study bringing together the available evaluated (ENDF) and unevaluated (EXFOR) data is highly advisable. TALYS should be tuned to maximize the global agreement. **TOP-PRIORITY**

* **On-the-fly calculation of interaction rates**:
Currently all interactions use lookup tables.
Taking into account the evolution of photon fields requires the tabulation of interaction rates for different redshifts.
An alternative would be the on-the-fly calculation of interaction rates using a sufficiently fast method such as Romberg integration.
The advantages compared to the lookup tables would be an increased accuracy (no interpolation necessary in energy or redshift) and a reduced memory requirement.
Even if the calculations are not done on-the-fly, they could performed initially, which would allow to ship CRPropa with just the tabulated cross-sections and the description of photon fields.

* **Photons from photo-disintegration**:
Currently, photons from photo-disintegration of cosmic ray nuclei are neglected.
Taking these photons into account would require a suitable implementation of predicted photon energies from e.g. TALYS **(in progress)**

* **Galactic magnetic field models**: Implementation of new models for the galactic magnetic field such as the revised random field of the JF12 model described in [arXiv:1409.5120](http://arxiv.org/abs/arXiv:1409.5120).

* **Inverse propagation**: Inverting the stochastic interactions and continuous energy losses would allow back-propagating observed cosmic rays through extragalactic space.

* **Synergies with Geant4** Geant4 basically provides the same physics as CRPropa, allowing for comparisons and possible technical/physical improvements. Additionally, event biasing techniques that are used in Geant could allow for more efficient MC sampling.

* **Replacement of SOPHIA** SOPHIA (PhotoPionProduction module) is generally the slowest part of CRPropa so it should be replaced with tabular data or some other idea.