
.. _file__Users_rab_Dropbox_softwares_CRPropa_CRPropa3-dev_include_crpropa_GridTools.h:

File GridTools.h
================


Grid related functions: load, dump, save, retrieve grid properties ... 



.. contents:: Contents
   :local:
   :backlinks: none

Definition (``/Users/rab/Dropbox/softwares/CRPropa/CRPropa3-dev/include/crpropa/GridTools.h``)
----------------------------------------------------------------------------------------------




Detailed Description
--------------------

This file contains a number of functions related to scalar and vector grids (Grid.h).

Dump/load functions are available for saving/loading grids to/from and binary and plain text files. In the files the grid points are stored from (0, 0, 0) to (Nx, Ny, Nz) with the z-index changing the fastest. Vector components are stored per grid point in xyz-order. In case of plain-text files the vector components are separated by a blank or tab and grid points are stored one per line. All functions offer a conversion factor that is multiplied to all values. 




Includes
--------


- ``array``

- ``crpropa/Grid.h`` (:ref:`file__Users_rab_Dropbox_softwares_CRPropa_CRPropa3-dev_include_crpropa_Grid.h`)

- ``crpropa/magneticField/MagneticField.h`` (:ref:`file__Users_rab_Dropbox_softwares_CRPropa_CRPropa3-dev_include_crpropa_magneticField_MagneticField.h`)

- ``string``



Included By
-----------


- :ref:`file__Users_rab_Dropbox_softwares_CRPropa_CRPropa3-dev_include_CRPropa.h`




Namespaces
----------


- :ref:`namespace_crpropa`


Functions
---------


- :ref:`exhale_function_group__Core_1gad5b5ba4c8f79b6b797c57d7ee658c1a6`

- :ref:`exhale_function_group__Core_1ga2be301c04cc62be7b1efdc3d8472dcda`

- :ref:`exhale_function_group__Core_1gad230c08955c7a505eefe413366f4717a`

- :ref:`exhale_function_group__Core_1ga25cabf4d5be84ea94f9bd90d50b17a13`

- :ref:`exhale_function_group__Core_1gaf89e81d370fe811cc9386261f02166a5`

- :ref:`exhale_function_group__Core_1ga0501e5b31dd2cdb5159459d874f16d30`

- :ref:`exhale_function_group__Core_1gab7c740cf66f560083029b8a6196a46a4`

- :ref:`exhale_function_group__Core_1ga091a17cca41dba2b79ec383f7f9e8f37`

- :ref:`exhale_function_group__Core_1gab66dd49f6fc38ee6834eff297e730ae3`

- :ref:`exhale_function_group__Core_1gaad9520b837fa3ff060ed1786da8c7f52`

- :ref:`exhale_function_group__Core_1ga35a4fc74b7deadd58ff6209df14140c0`

- :ref:`exhale_function_group__Core_1ga095a3cc70d79ca2a58df236d823fc7f2`

- :ref:`exhale_function_group__Core_1ga14e9170907096634f33d71ec5aae75d3`

- :ref:`exhale_function_group__Core_1ga788ff1002bda23d428e0e43ac5b7cf3a`

- :ref:`exhale_function_group__Core_1ga42f4cdca40b5d3f1632ac93a886a8781`

- :ref:`exhale_function_group__Core_1gad823c8cf8555f1b45432468adef2d404`

- :ref:`exhale_function_group__Core_1ga0d66e0a4632d21364c4227784074f53e`

- :ref:`exhale_function_group__Core_1gad86fe3acb119896b39adb6d751b57700`

