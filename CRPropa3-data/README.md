CRPropa3-data
=============

Tools to generate the interaction and nuclear data for [CRPropa 3](https://github.com/CRPropa/CRPropa3).

Interactions between cosmic ray nuclei and background photons
 - calc_epp.py : electron pair production
 - calc_ppp.py : photo-pion production
 - calc_photodisintegration.py : photodisintegration
 - calc_elasticscattering.py   : elastic scattering on the nuclear structure

Interactions between cosmic ray photons / electrons and background photons
 - calc_electromagnetic.py
    - photon    : pair and double-pair production
    - electrons : triplet pair production and inverse Compton scattering

Other processes:
 - calc_decay.py : nuclear decays
 - calc_synchrotron.py : synchrotron radiation of charged particles
 - calc_scaling : global redshift scaling of cosmic photon fields
 - calc_mass : table of nuclear masses

Helper modules
 - photonField.py     : collection of background photon fields (CMB, EBL, URB)
 - interactionRate.py : functions to calculate interaction rates with isotropic photon fields
