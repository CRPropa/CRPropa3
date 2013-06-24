#ifndef CRPROPA_PARTICLE_ID_H
#define CRPROPA_PARTICLE_ID_H

#include <cstddef>
#include <vector>

namespace crpropa {

std::vector<int> neutrinos();
std::vector<int> leptons();

/** This implements the 2012 Monte Carlo nuclear code scheme.
 * Ion numbers are +/- 10LZZZAAAI.
 * AAA is A - total baryon number
 * ZZZ is Z - total charge
 * L is the total number of strange quarks.
 * I is the isomer number, with I=0 corresponding to the ground state.
 */
int nucleusId(int a, int z);
int chargeNumberFromNucleusId(int id);
int massNumberFromNucleusId(int id);

/* CRPropa2.0 code scheme */
int convertFromCRPropaId(int id);
int convertToCRPropaId(int id);

} // namespace crpropa

#endif /* CRPROPA_PARTICLE_ID_H */
