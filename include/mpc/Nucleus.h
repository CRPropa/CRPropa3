#ifndef NUCLEUS_H_
#define NUCLEUS_H_

#include "mpc/AssocVector.h"
#include <stdexcept>

namespace mpc {

// HepPID codes for ions 1AAAZZZ00nj
// see http://cepa.fnal.gov/psm/HepPID
// nj = 2J+1 is the nucleus total spin and is currently neglected
inline int getNucleusId(int a, int z) {
	if (z < 0)
		throw std::runtime_error("mpc::Nucleus: no nucleus with Z < 0");
	if (a < 0)
		throw std::runtime_error("mpc::Nucleus: no nucleus with A < 0");
	if (a < z)
		throw std::runtime_error("mpc::Nucleus: no nucleus with A < Z");
	return 1000000000 + a * 1000000 + z * 1000;
}

inline int getChargeNumberFromNucleusId(int id) {
	return ((id - 1000000000) % 1000000) / 1000;
}

inline int getMassNumberFromNucleusId(int id) {
	return (id - 1000000000) / 1000000;
}

static Loki::AssocVector<int, double> nuclearMassTable;
void initNuclearMassTable();
double getNucleusMass(int id);

} // namespace mpc

#endif /* NUCLEUS_H_ */
