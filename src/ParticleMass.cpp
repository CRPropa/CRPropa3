#include "crpropa/ParticleMass.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Common.h"
#include "crpropa/Units.h"

#include "kiss/convert.h"
#include "kiss/logger.h"

#include <vector>
#include <fstream>
#include <stdexcept>
#include <limits>

namespace crpropa {

struct NuclearMassTable {
	bool initialized;
	std::vector<double> table;

	NuclearMassTable() {
		initialized = false;
	}

	void init() {
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
		initialized = true;
	}

	double getMass(std::size_t idx) {
		if (!initialized) {
#pragma omp critical(init)
			init();
		}
		return table[idx];
	}
};

static NuclearMassTable nuclearMassTable;

double nuclearMass(int id) {
	int A = massNumber(id);
	int Z = chargeNumber(id);
	return nuclearMass(A, Z);
}

double nuclearMass(int A, int Z) {
	if ((A < 1) or (A > 56) or (Z < 0) or (Z > 26) or (Z > A)) {
		KISS_LOG_WARNING <<
		"nuclearMass: nuclear mass not found in the mass table for " <<
	        "A = " << A << ", Z = " << Z << ". " <<
		"Approximated value used A * amu - Z * m_e instead.";
		return A * amu - Z * mass_electron;
	}
	int N = A - Z;
	return nuclearMassTable.getMass(Z * 31 + N);
}

} // namespace crpropa
