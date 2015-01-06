#include "crpropa/ParticleMass.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Common.h"

#include <kiss/convert.h>

#include <vector>
#include <fstream>
#include <stdexcept>
#include <limits>

namespace crpropa {

struct NuclearMassTable {
	std::vector<double> table;

	NuclearMassTable() {
		std::string filename = getDataPath("nuclear_mass.txt");
		std::ifstream infile(filename.c_str());

		if (!infile.good())
			throw std::runtime_error("crpropa: could not open file " + filename);

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

double nuclearMass(int id) {
	int A = massNumber(id);
	int Z = chargeNumber(id);
	return nuclearMass(A, Z);
}

double nuclearMass(int A, int Z) {
	int N = A - Z;
	double mass = nuclearMassTable.table[Z * 31 + N];
	if (mass == 0)
		throw std::runtime_error("nuclearMass: nucleus not found " + kiss::str(A) + ", " + kiss::str(Z));
	return mass;
}

} // namespace crpropa
