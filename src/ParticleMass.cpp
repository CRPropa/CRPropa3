#include "mpc/ParticleMass.h"
#include "mpc/ParticleID.h"
#include "mpc/Common.h"

#include <kiss/convert.h>

#include <vector>
#include <fstream>
#include <stdexcept>
#include <limits>

namespace mpc {

struct NuclearMassTable {
	std::vector<double> table;

	NuclearMassTable() {
		std::string filename = getDataPath("nuclear_mass.txt");
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

double nucleusMass(int id) {
	int Z = chargeNumberFromNucleusId(id);
	int N = massNumberFromNucleusId(id) - Z;
	double mass = nuclearMassTable.table[Z * 31 + N];
	if (mass == 0)
		throw std::runtime_error("getNucleusMass: nucleus not found " + kiss::str(id));
	return mass;
}

} // namespace mpc
