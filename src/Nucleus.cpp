#include "mpc/Nucleus.h"
#include "mpc/Common.h"

#include <kiss/convert.h>
#include <HepPID/ParticleIDMethods.hh>

#include <vector>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace mpc {

int getNucleusId(int a, int z) {
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

int getChargeNumberFromNucleusId(int id) {
	return HepPID::Z(id);
}

int getMassNumberFromNucleusId(int id) {
	return HepPID::A(id);
}

int convertFromCRPropaId(int crp_id) {
	int Z = crp_id / 1000;
	int A = crp_id % 1000;
	return getNucleusId(A, Z);
}

int convertToCRPropaId(int id) {
	int Z = getChargeNumberFromNucleusId(id);
	int A = getMassNumberFromNucleusId(id);
	return Z * 1000 + A;
}

struct NuclearMassTable {
	std::vector<double> table;

	NuclearMassTable() {
		std::string filename = getDataPath("/NuclearMass/nuclearMassTable.txt");
		std::ifstream infile(filename.c_str());

		if (!infile.good())
			throw std::runtime_error("mpc: could not open file " + filename);

		int Z, N;
		double mass;
		while (infile.good()) {
			if (infile.peek() != '#') {
				infile >> Z >> N >> mass;
				table.push_back(mass);
			}
			infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}

		infile.close();
	}
};

static NuclearMassTable nuclearMassTable;

double getNucleusMass(int id) {
	int Z = HepPID::Z(id);
	int N = HepPID::A(id) - Z;
	double mass = nuclearMassTable.table[Z * 31 + N];
	if (mass == 0)
		throw std::runtime_error("mpc: nucleus not found " + kiss::str(id));
	return mass;
}

} // namespace mpc
