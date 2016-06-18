CRPropa needs a number of data files to run. These include databases for nuclear mass and decay rates, interaction rates for photodisintegration, photo-pion production and electron-pair production, as well as the data files for DINT and EleCa.
The scripts and files that are used to prepare the data files are tracked by the [CRPropa3-data repository](https://github.com/CRPropa/CRPropa3-data).  

Since a git repository is not well suited to store these large files, they are downloaded automatically during cmake configuration. In case you want to use a different photodisintegration model use the link to the alternative tarballs below and extract them to your install folder. 

**Default data files:** https://crpropa.desy.de/images/8/83/Data_2016_05_23.tar.gz  
Contains interaction rates nuclei, electrons and photons and spectra of secondary particles from these interactions, as well as data on nuclear masses and decay rates.  
Uses photodisintegration cross sections (for A > 12) from TALYS 1.6 with a partially modified set of GDR parameters.

**Alternative 1:** https://crpropa.desy.de/images/d/d3/Data_2016_03_21_Kossov.tar.gz  
Photodisintegration interaction rates for the Kossov model.

**Alternative 2:** https://crpropa.desy.de/images/0/0a/Data_2016_03_21_PSB.tar.gz  
Photodisintegration interaction rates for the PSB model.