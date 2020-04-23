## Basic Concepts

### Modular design
The basic concept of CRPropa 3 is having all aspects of cosmic ray propagation,
such as photodisintegration or the maximum trajectory length, split up into
independent simulation '''Module'''s. Their task is to modify a cosmic ray one
after the other and in small steps. The propagation then consists of repeatedly
looping a cosmic ray through a list of active modules ('''ModuleList''') until
a module signals that the propagation is finished.
The modules are independent in that they do not require each other and can basically be used in any combination.

The object that contains the initial and current state on the cosmic ray
particle is called '''Candidate'''. It also contains information on the state
in the previous simulation step, as well as the current redshift, the current
and next (preliminary) step size, status flags (e.g. "Detected") and a record
of interactions that are projected to happen to the cosmic ray particle.

### Parallelization/Multiprocessing
CPRopa 3 enables shared memory multiprocessing (using OpenMP), thus the
distribution of workload among several computing units ("threads"), which share
the same memory.
Distributing the workload is trivial in cosmic ray propagation as there is no
interaction between the cosmic rays. To propagate e.g. 20 cosmic rays with two
threads, thread A could propagate the first 10, and thread B the second 10
simultaneously.

To make this work in practice, "thread-safety" has to be ensured. Objects and
variables which are shared among the threads are not allowed to change during
the parallel sections of the code. If they do, a change by thread A is likely
to mess with the correct calculation in thread B. Most importantly, the
simulation modules which are in the current parallelization layout shared among
the threads cannot store information on the cosmic ray. Instead the information
has be stored on the cosmic ray itself.

Program parts in which "thread-safety" cannot be ensured need to be put inside
a "critical region". If a thread reaches a critical region that is occupied by
another thread, it waits for the thread to finish. Critical regions limit the
benefit of parallelization and should not contain computationally expensive
code.
Examples for critical region are in the output modules or the call to the
Fortran code SOPHIA for the photopion-production (although not checked, SOPHIA
is very likely not thread-safe).

Other paradigms for parallelization are distributed memory multiprocessing with
MPI, and GPGPU computing. These require a specific program layout and are
currently not supported.

### Unit system
In UHECR propagation many physical domains meet, each with their specific unit
system. As any choice is arbitrary, the '''SI-system''' is used internally
throughout the program. Symbols are provided for the most important units (eV,
kpc, Mpc, ...). This enforces expressive code, e.g.
MaximumTrajectoryLength(1000 * Mpc).

### Particle ID
The PDG Monte Carlo numbering scheme (see
http://pdg.lbl.gov/2012/mcdata/mc_particle_id_contents.html) is used for
compatibility with other Monte Carlo programs.
Important elementary particle IDs are
* 11 electron
* 12 electron neutrino
* 13 muon
* 14 muon neutrino
* 22 photon

Antiparticles have negative IDs, e.g. -11 for anti-electron. <br>
Nuclear codes are given as ±10LZZZAAAI with
* charge number Z,
* mass number A,
* ismoer level I (ground state I = 0, I > 0 excitation states) and
* total number of strange quarks L

### Reference Counting
For memory management CRPropa3 uses reference counting with the template class '''Referenced'''.
This is necessary for compatibility with Python's memory management.

