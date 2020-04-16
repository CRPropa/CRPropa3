CRPropa needs a number of data files to run. These include databases for nuclear mass and decay rates, interaction rates for photodisintegration, photo-pion production and electron-pair production, as well as the data files for DINT and EleCa.
The scripts and files that are used to prepare the data files are tracked by the [CRPropa3-data repository](https://github.com/CRPropa/CRPropa3-data).  

Since a git repository is not well suited to store these large files, they are downloaded automatically during cmake configuration. In case you want to use a different photodisintegration model or want to download the default data file manually use the following [link](https://www.desy.de/~crpropa/data/interaction_data/) to the tarballs files and extract them to your install folder. 

**Default data files:** (data.tar.gz)
Contains interaction rates nuclei, electrons and photons and spectra of secondary particles from these interactions, as well as data on nuclear masses and decay rates.  
Uses photodisintegration cross sections (for A > 12) from TALYS 1.8 with a partially modified set of GDR parameters.

**Alternative 1:** (data_2016_03_21_Kossov.tar.gz)
Photodisintegration interaction rates for the Kossov model.

**Alternative 2:** (data_2016_03_21_PSB.tar.gz)
Photodisintegration interaction rates for the PSB model.

### Verify the data files integrity

Every data file should have a corresponding -CHECKSUM file in which appropriate MD5 sums are stored. A checksum file can be generated with, for example:
```
$ ls data.tar.gz | xargs -I{} sh -c 'md5sum "$1" > "$1-CHECKSUM"' -- {}
```
To verify the integrity of the data files, a it is enough to download the checksum file in the same directory with the data file:
```
$ md5sum -c *-CHECKSUM
```
Files that are automatically downloaded with CMake are also automatically verified by CMake.