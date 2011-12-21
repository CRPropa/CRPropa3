#ifndef HEPPID_H_
#define HEPPID_H_

namespace mpc {

// HepPID codes for ions 1AAAZZZ00n, where n contains the total spin as n = 2J+1
// see http://cepa.fnal.gov/psm/HepPID/
inline int getNucleusId(int a, int z) {
	return 1000000000 + z * 10000 + a * 10;
}

inline int getChargeNumberFromNucleusId(int id) {
	return ((id - 1000000000) % 1000000) / 1000;
}

inline int getMassNumberFromNucleusId(int id) {
	return (id - 1000000000) / 1000000;
}

}

#endif /* HEPPID_H_ */
