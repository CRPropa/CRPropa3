CRPropa needs a number of data files to run. These include databases for nuclear mass and decay rates, interaction rates for photodisintegration, photo-pion production and electron-pair production, as well as the data files for DINT and EleCa.

Since a git repository is not well suited to store these large files, they are currently stored on https://crpropa.desy.de/Special:UncategorizedFiles and downloaded automatically during cmake configuration.

The scripts and files that are used to prepare the data files are tracked by the [CRPropa3-data repository](https://github.com/CRPropa/CRPropa3-data).
When changes are made, a new tarball of all data files should be created and uploaded to e.g. https://crpropa.desy.de/.
The new URL and MDF5 checksum then have to be set accordingly in CMakeLists.txt.