Galactic Cosmic Rays
--------------------

Propagation of Galactic cosmic rays can be done with the single-particle 
approach (solving the equation of motion) or with an ensemble averaged 
ansatz (solving the transport equation for a particle distribution). The 
latter one is explained in section "Diffusion of Cosmic Rays". The single-particle
approach is demonstrated in the following example.

.. toctree::
   example_notebooks/galactic_trajectories/galactic_trajectories.ipynb


Diffusion of Cosmic Rays
------------------------

Diffusion of cosmic rays can be modeled in an ensemble averaged way
using transport equations. In CRPropa this is done based on stochastic
differential equations (SDEs) which are mathematically equivalent to 
the more familar transport equations.

The first notebooks give an overview how to set up a simulation using
a user-defined diffusion coefficient.

.. toctree::

   example_notebooks/Diffusion/DiffusionValidationI.ipynb
   example_notebooks/Diffusion/DiffusionValidationII.ipynb

Advection and Adiabatic Energy Changes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Advection in non-divergence-free velocity fields causes adiabatic energy
changes, which can be modeled following the next examples.

.. toctree::
   example_notebooks/Diffusion/AdiabaticCooling.ipynb

Diffusive Shock Acceleration
^^^^^^^^^^^^^^^^^^

Diffusive shock acceleration or first order Fermi acceleration can be modeled
in the diffusive picture as an interplay between diffusion, advection and 
adiabatic cooling. Simple shock configurations are shown in the following 
example notebook.

.. toctree::
   example_notebooks/acceleration/diffusive_shock_acceleration.ipynb

In the next example time-dependent advection fields are introduced and acceleration
at a moving shock and stationary shock are compared.

.. toctree::
    example_notebooks/acceleration/DSA-moving-shock.ipynb

Momentum Diffusion 
^^^^^^^^^^^^^^^^^^

Momentum diffusion or second order Fermi acceleration is explained in the
following example notebook.

.. toctree::
   example_notebooks/Diffusion/MomentumDiffusion.ipynb

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
   example_notebooks/density/density_grid_sampling.ipynb
