#include "crpropa/ParticleID.h"

#include <HepPID/ParticleIDMethods.hh>
#include <kiss/convert.h>

namespace crpropa {

std::vector<int> neutrinos() {
	int a[] = { 12, -12, 14, -14, 16, -16 };
	std::vector<int> v(a, a + sizeof(a) / sizeof(int));
	return v;
}

std::vector<int> leptons() {
	int a[] = { 11, -11, 12, -12, 13, -13, 14, -14, 15, -15, 16, -16 };
	std::vector<int> v(a, a + sizeof(a) / sizeof(int));
	return v;
}

int nucleusId(int a, int z) {
	if (z < 0)
		throw std::runtime_error(
				"crpropa::Nucleus: no nucleus with Z < 0, A=" + kiss::str(a) + " Z="
						+ kiss::str(z));
	if (a < 1)
		throw std::runtime_error(
				"crpropa::Nucleus: no nucleus with A < 1, A=" + kiss::str(a) + " Z="
						+ kiss::str(z));
	if (a < z)
		throw std::runtime_error(
				"crpropa::Nucleus: no nucleus with A < Z, A=" + kiss::str(a) + " Z="
						+ kiss::str(z));
	return 1000000000 + z * 10000 + a * 10;
}

int chargeNumberFromNucleusId(int id) {
	return HepPID::Z(id);
}

int massNumberFromNucleusId(int id) {
	return HepPID::A(id);
}

bool isNucleus(int id) {
	return HepPID::isNucleus(id);
}

int convertFromCRPropa2NucleusId(int crp_id) {
	int Z = crp_id / 1000;
	int A = crp_id % 1000;
	return nucleusId(A, Z);
}

int convertToCRPropa2NucleusId(int id) {
	if (not HepPID::isNucleus(id)) // if not nucleus, just return the id
		return id;
	int Z = chargeNumberFromNucleusId(id);
	int A = massNumberFromNucleusId(id);
	return Z * 1000 + A;
}

}
