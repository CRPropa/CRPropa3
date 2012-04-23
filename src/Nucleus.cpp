#include "mpc/Nucleus.h"
#include "mpc/Common.h"

#include <kiss/convert.h>
#include <HepPID/ParticleIDMethods.hh>

#include <vector>
#include <fstream>
#include <limits>

namespace mpc {

static std::vector<double> nuclearMassTable;

void initNuclearMassTable() {
	std::string filename = getDataPath("/NuclearMass/nuclearMassTable.txt");
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error("mpc: could not open file " + filename);

	int Z, N;
	double mass;
	while (infile.good()) {
		if (infile.peek() != '#') {
			infile >> Z >> N >> mass;
			nuclearMassTable.push_back(mass);
		}
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}

	infile.close();
}

double getNucleusMass(int id) {
	if (nuclearMassTable.size() == 0)
		initNuclearMassTable();

	int Z = HepPID::Z(id);
	int N = HepPID::A(id) - Z;
	double mass = nuclearMassTable[Z * 31 + N];
	if (mass == 0)
		throw std::runtime_error("mpc: nucleus not found " + kiss::str(id));
	return mass;
}

} // namespace mpc
