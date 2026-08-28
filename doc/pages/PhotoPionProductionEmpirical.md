# Empirical photomeson data tables

`PhotoPionProductionEmpirical` implements the empirical photomeson model of

> L. Morejón, A. Fedynitch, D. Boncioli, D. Biehl, W. Winter,
> *JCAP* **11** (2019) 007, [arXiv:1904.07999](https://arxiv.org/abs/1904.07999)

as published in [AstroPhoMes](https://github.com/mohller/AstroPhoMes). It is an
alternative to `PhotoPionProduction`; this page covers how to use it, the data it
needs, and how the port was validated.

## The three files

The module relies on three data files contained in `data/PhotoPionProductionEmpirical/`:

| File | Contents | Size |
|---|---|---|
| `basis.txt` | cross-section basis vectors on the 159-point `eps_r` grid | 19 kB |
| `redistribution.txt` | SOPHIA x-spectra of the secondaries | 2.1 MB |
| `fragments.txt` | exclusive fragmentation channels | 6.7 MB |

**All three are independent of the photon field.** Unlike `PhotoPionProduction`
there is no `rate_<field>.txt`: the module integrates the field itself at
construction. So nothing has to be regenerated when a photon field is added or
changed, and photon fields defined by the user in Python (having no table
defined) work directly.

Computing the rates for the bais cross section takes about **~8 ms**, parsing the 8.8 MB of tables takes **~220 ms**, so a module is built in **~230 ms** total.

## Usage

The module is an alternative to `PhotoPionProduction` without calls to SOHPIA. **Do not add both: they model the same process and would double the interaction rate.**

```python
import crpropa as crp

cmb = crp.CMB()
ppe = crp.PhotoPionProductionEmpirical(cmb, True, True, True)

sim = crp.ModuleList()
sim.add(crp.SimplePropagation(1 * crp.kpc, 10 * crp.Mpc))
sim.add(ppe)
sim.add(crp.PhotoDisintegration(cmb))
sim.add(crp.NuclearDecay())
sim.add(crp.MinimumEnergy(1 * crp.EeV))
sim.add(crp.MaximumTrajectoryLength(100 * crp.Mpc))

source = crp.Source()
source.add(crp.SourceParticleType(crp.nucleusId(56, 26)))
source.add(crp.SourceEnergy(2e4 * crp.EeV))
source.add(crp.SourcePosition(crp.Vector3d(0)))
source.add(crp.SourceDirection(crp.Vector3d(1, 0, 0)))

out = crp.ParticleCollector()
sim.add(out)
sim.run(source, 50, True)
```

Photopion production needs `log10(gamma) >~ 11`, and for a nucleus
`gamma ~ E / (A m_N)`, so iron has to carry several thousand EeV before it
interacts at all, hence the unphysically high source energy above.

Unlike in the superposition approach of `PhotoPionProduction`, the Empirical Model's  fragmentation table produces a larger range of masses, and many short-lived residuals that need to decay, so `NuclearDecay` is always required.

### Runtime options

All are optional. The defaults are what the validation below refers to; the

| Option | Default | Alternative | Set the alternative when |
|---|---|---|---|
| `setMultiplicitySampling()` | `FloorBernoulli`: `floor(m)` secondaries plus one more with probability `m - floor(m)`; exact mean at minimum variance | `PoissonSampling`: Poisson around the mean | You specifically want Poisson variance in the secondary counts. Be aware it manufactures multi-pion states near the Delta resonance, where the mean pion yield is ~0.4 and the kinematics allow exactly one: ~11% of events there get two or more pi0. |
| `setDecayMode()` | `PromptDecay`: daughters created at the interaction vertex | `DisplacedDecay`: daughters offset by `gamma c tau` along the parent direction | You need decay *positions* rather than vertex positions. The offset is geometric only: no energy loss, deflection or cooling in flight. It is negligible for extragalactic propagation (a 1 EeV muon travels 2e-10 Mpc) and misleading in compact sources, where cooling overtakes decay above `E_crit ~ 5.9 EeV / B[G]`. |
| `setHaveKaons()` | `True`: K+- produced and decayed via a renormalised two-mode approximation (75.5% mu nu, 24.5% pi pi0) | `False`: no kaons at all | You want to isolate how much of your neutrino flux comes through kaons. Not an accuracy improvement: kaons are ~10% of the charged-pion yield at the highest energies, so disabling them removes that contribution entirely. |
| `setConserveCharge()` | `True`: net meson charge is repaired against the ejected nucleons, then the residual, then the mesons | `False`: no repair; events may carry a net charge error | Essentially never, outside debugging. Use `getChargeViolationCount()` to see how often the repair could not close, which should be zero. |
| `setLimit()` | `0.1`: step limited to this fraction of the mean free path | any `double` in (0, 1] | The usual step-size/accuracy trade-off, as for any interaction module. |
| `setInteractionTag()` | `"PPPE"` | any string | You want to trace these secondaries separately in the output, e.g. to separate them from `"PPP"` when comparing the two photopion modules in one run. |

The two enumerated options are nested in the class:

```python
E = crp.PhotoPionProductionEmpirical
ppe.setMultiplicitySampling(E.PoissonSampling)   # E.FloorBernoulli is the default
ppe.setDecayMode(E.DisplacedDecay)               # E.PromptDecay is the default
```

### Diagnostics

Useful for checking the model rather than running it:

```python
eps = ppe.getEnergyBin(140)                  # rest-frame photon energy [J]
ppe.crossection(eps, 56, 26)                 # nonelastic cross section [m^2]
ppe.pionScaling(eps, 56, 26)                 # empirical pion scaling, 1 for A=1
ppe.getRate(56, 26, 1e12, 0)                 # interaction rate [1/m]
ppe.nucleusMFP(56, 26, 1e12, 0)              # mean free path [m]
ppe.getChargeViolationCount()                # should stay 0
ppe.getMeanFragmentBudget(56, 26)            # (<dA>, <dZ>) per interaction
```

`getMeanFragmentBudget(A, Z)` returns the pair `(<dA>, <dZ>)`: the mean mass
and charge lost per interaction. For Fe-56 it gives 5.79, against exactly 1 for
the superposition treatment of `PhotoPionProduction`.

## How are the data files produced?

Like every other CRPropa data file, the generation takes place in
[CRPropa3-data](https://github.com/CRPropa/CRPropa3-data):

```
dataGeneration/PhotoPionProductionEmpirical/export_tables.py
```

which reads the model from AstroPhoMes, as dependency submodule under
`lib/AstroPhoMes`. That script producing the tables is the source to
consult for file formats, unit conventions and acceptance ranges.

Two things to know before regenerating:

- **The 159-point `eps_r` grid must not be re-gridded.** AstroPhoMes blends its
  data-based pion normalisation in and out with sigmoids whose windows are
  hardcoded to grid *indices* (0–31 and 55–104), not to energies. A different
  grid would change the pion yield. The grid is also not log-equidistant.
- `pion_spl` is in **microbarn**, and is deliberately **not** clipped even where
  its extrapolation goes negative (see the notes in `initBasis`).

## Coverage

Mass coverage follows whatever nuclide range `fragments.txt` carries; the module is
written to be agnostic to it. It accepts up to the limits `NUCLEAR_ZMAX` (82) and `NUCLEAR_NMAX`
(132), and missing species contain empty table rows, so extending the model is a data swap with no C++
change. Rows for nuclei outside CRPropa's range are skipped at load.

The current tables reach lead. Nuclei the table does not cover simply do not
interact through this module.

## Possible cause of failure

```
PhotoPionProductionEmpirical: could not open .../basis.txt
```

means the data release predates these tables. Either update to one that contains
`PhotoPionProductionEmpirical/`, or point `CRPROPA_DATA_PATH` at a CRPropa3-data
checkout in which the generator has been run:

```sh
export CRPROPA_DATA_PATH=/path/to/CRPropa3-data/data
```

## Validation

The module is checked against three independent references, all automated in
`testInteraction.cpp` and `testPhotoMesonEmpirical.py`:

- `sigma_nonel` reproduces AstroPhoMes to 5e-10, and the pion scaling to 7.5e-3,
  for Fe-56, O-16 and Pb-208;
- the runtime-built proton rate reproduces CRPropa's shipped SOPHIA rate table to
  within 4%, the residual being the two codes discretising the same integral
  differently;
- for A = 1 the empirical model *is* SOPHIA, and both modules agree on the final
  state (leading nucleon, neutrino and photon spectra) to better than 1%.

Mass, charge and energy are conserved in every event by construction.

## Citation

If you use this module, please cite the source paper and the code it is ported
from, alongside the usual CRPropa reference. The empirical relations, the mass
scalings and the fragmentation table are the content of that paper, citing it
lets a reader trace the physical assumptions behind your results back to
where they were derived and validated.

> L. Morejón, A. Fedynitch, D. Boncioli, D. Biehl, W. Winter,
> *Improved photomeson model for interactions of cosmic ray nuclei*,
> JCAP **11** (2019) 007, [doi:10.1088/1475-7516/2019/11/007](https://doi.org/10.1088/1475-7516/2019/11/007);
> [arXiv:1904.07999](https://arxiv.org/abs/1904.07999)

> L. Morejón, *AstroPhoMes: Photomeson models for astrophysical nuclei* (2019),
> [DOI:10.5281/zenodo.2600177](https://doi.org/10.5281/zenodo.2600177)

The secondary spectra are SOPHIA redistribution tables, so the SOPHIA reference
applies as it does for `PhotoPionProduction`. The empirical relations themselves
are reproduced from Terranova and Tavares,
[Phys. Scr. 49 (1994) 267](https://doi.org/10.1088/0031-8949/49/3/004).
