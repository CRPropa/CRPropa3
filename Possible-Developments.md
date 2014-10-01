* **Galactic magnetic field models**: Implementation of new models for the galactic magnetic field such as the revised random field of the JF12 model described in [arXiv:1409.5120](http://arxiv.org/abs/arXiv:1409.5120).

* **On-the-fly calculation of interaction rates**:
Currently all interactions use lookup tables.
Taking into account the evolution of photon fields requires the tabulation of interaction rates for different redshifts.
An alternative would be the on-the-fly calculation of interaction rates using a sufficiently fast method such as Romberg integration.
The advantages compared to the lookup tables would be an increased accuracy (no interpolation necessary in energy or redshift) and a reduced memory requirement.

* **Photons from photo-disintegration**:
Currently, photons from photo-disintegration of cosmic ray nuclei are neglected.
Taking these photons into account would require a suitable implementation of predicted photon energies from e.g. TALYS

* **Direct implementation of EleCa**:
EleCa is currently used after the actual CRPropa simulation.
Since EleCa is a MC code for individual photons, the functionality could be directly implemented in CRPropa in the form of simulation modules. This would allow propagating photons in the same way as nuclei and neutrinos, using the potentially more efficient CRPropa framework and considerably reducing code complexity.

* **Inverse propagation**: Inverting the stochastic interactions and continuous energy losses would allow back-propagating observed cosmic rays through extragalactic space.

* **Synergies with Geant4** Geant4 basically provides the same physics as CRPropa, allowing for comparisons and possible technical/physical improvements. Additionally, event biasing techniques that are used in Geant could allow for more efficient MC sampling.