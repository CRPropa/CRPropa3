Diffusion of Cosmic Rays
---------------------------------

Diffusion of cosmic rays can be modeled in an ensemble averaged way
using transport equations. In CRPropa this is done based on stochastic
differential equations (SDEs) which are mathematically equivalent to 
the more familar transport equations.

The first notebooks give an overview how to set up a simulation using
a user-defined diffusion coefficient.

.. toctree::

   example_notebooks/Diffusion/DiffusionValidationI.v4.ipynb
   example_notebooks/Diffusion/DiffusionValidationII.v4.ipynb

Advection and Adiabatic Energy Changes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Advection in non-divergence-free velocity fields causes adiabatic energy
changes, which can be modeled following the next examples.

.. toctree::
   example_notebooks/Diffusion/AdiabaticCooling.ipynb


Example of diffusion in the Milky way
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   example_notebooks/Diffusion/GalacticDiffusion.ipynb

Gas Densities
^^^^^^^^^^^^^

The last two examples show how to use the gas density fields. They can be 
used to, e.g., model the CR source distribution and will later be used as 
target fields for hadron-hadron interation.

.. toctree::
   example_notebooks/density/density_easy_example.ipynb
   example_notebooks/density/density_model_example.ipynb
