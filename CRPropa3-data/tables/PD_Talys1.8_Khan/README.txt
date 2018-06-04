Overview
-----------

TALYS 1.8 is used together with the GDR parameter files in the folder TALYS-patch1 to
 - calculate all photodisintegration cross sections
 - obtaining the energies and emission probabilities of photons during these interactions
 - calculating the cross section for elastic scattering of the incident photon

The photodisintegration cross sections are processed with
 - A1_talys_photodisintegration.py
 - A2_process-photodisintegration.py
 - A3_process-photonemission.py

The elastic scattering cross sections are generated and processed with
 - B1_talys_elasticscattering.py
 - B2_process_elasticscattering.py

Note: The raw TALYS output is too large but the processed cross section files are tracked in the repository.

Theses files are used in
 - CRPropa3-data/calc_photodisintegration.py
 - CRPropa3-data/calc_elasticscattering.py
to calculate the interaction rates with various photon backgrounds.
