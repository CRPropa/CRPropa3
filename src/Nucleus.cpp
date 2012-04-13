#include "mpc/Nucleus.h"
#include "mpc/Common.h"

#include <kiss/convert.h>

#include <fstream>
#include <limits>

namespace mpc {

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
			nuclearMassTable[getNucleusId(Z + N, Z)] = mass;
		}
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}

	infile.close();
}

double getNucleusMass(int id) {
	if (nuclearMassTable.size() == 0)
		initNuclearMassTable();

	Loki::AssocVector<int, double>::const_iterator i = nuclearMassTable.find(id);
	if (i == nuclearMassTable.end())
		throw std::runtime_error("mpc: nuclear mass not found " + kiss::str(id));
	return i->second;
}

} // namespace mpc
