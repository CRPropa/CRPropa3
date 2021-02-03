Acceleration
-------------------------

The acceleration modules work conceptually similar to energy losses in photon
fields - instead of the energy loss from an interaction with a background
photon, the particle energy is changed by interacting with a scattering center.

In the scattering process, the particle is transformed into the restframe of the
scattering center. In this frame the particle is scattered in a random direction
and transformed back. The properties of the velocity field of scattering center
define if this is a second or first order Fermi process.

.. toctree::

   example_notebooks/acceleration/second_order_fermi.ipynb
   example_notebooks/acceleration/first_order_fermi_acceleration.ipynb




The basic simulation of acceleration as sketched above is agnostic about the
nature of the scattering centers. For a random walk analogue to diffusion in a
magnetic field the step length is proportional to the diffusion coefficient
which depends on the energy of the particle and the properties of the turbulent
field. The corresponding modification of the step length is implemented as
:cpp:class:`crpropa::StepLengthModifier` that can be added to the corresponding simulation
regimes.

The exact nature of the dependency is subject of ongoing research. Here,
so far only quasi-linear theory, a simplistic approach to diffusion of charged
particles in magnetic fields is implemented in the corresponding modifier :cpp:class:`crpropa::QuasiLinearTheory`.





