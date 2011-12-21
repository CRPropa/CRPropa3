#ifndef HEPPID_H_
#define HEPPID_H_

namespace mpc {

// HepPID codes for ions 1AAAZZZ00nj
// see http://cepa.fnal.gov/psm/HepPID
// nj = 2J+1 is the nucleus total spin and is currently neglected
inline int getNucleusId(int a, int z) {
	return 1e9 + a * 1e6 + z * 1e3;
}

inline int getChargeNumberFromNucleusId(int id) {
	return ((id - 1000000000) % 1000000) / 1000;
}

inline int getMassNumberFromNucleusId(int id) {
	return (id - 1000000000) / 1000000;
}

}

#endif /* HEPPID_H_ */
