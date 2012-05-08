#ifndef NUCLEUS_H_
#define NUCLEUS_H_

#include "kiss/convert.h"

#include <stdexcept>

namespace mpc {

// This implements the 2006 Monte Carlo nuclear code scheme.
// Ion numbers are +/- 10LZZZAAAI.
// AAA is A - total baryon number
// ZZZ is Z - total charge
// L is the total number of strange quarks.
// I is the isomer number, with I=0 corresponding to the ground state.
inline int getNucleusId(int a, int z) {
	if (z < 0)
		throw std::runtime_error(
				"mpc::Nucleus: no nucleus with Z < 0, A=" + kiss::str(a) + " Z="
						+ kiss::str(z));
	if (a < 1)
		throw std::runtime_error(
				"mpc::Nucleus: no nucleus with A < 1, A=" + kiss::str(a) + " Z="
						+ kiss::str(z));
	if (a < z)
		throw std::runtime_error(
				"mpc::Nucleus: no nucleus with A < Z, A=" + kiss::str(a) + " Z="
						+ kiss::str(z));
	return 1000000000 + z * 10000 + a * 10;
}

void initNuclearMassTable();
double getNucleusMass(int id);

} // namespace mpc

#endif /* NUCLEUS_H_ */
