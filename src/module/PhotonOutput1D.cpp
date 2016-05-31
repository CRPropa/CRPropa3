#include "crpropa/module/PhotonOutput1D.h"
#include "crpropa/Units.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

namespace crpropa {

PhotonOutput1D::PhotonOutput1D(const string &filename) :
		filename(filename), output(filename.c_str()) {
	output << "#ID\tE\tD\tpID\tpE\tiID\tiE\n";
	output << "#\n";
	output << "# ID          Id of particle (photon, electron, positron)\n";
	output << "# E           Energy [EeV]\n";
	output << "# D           Comoving distance to origin [Mpc]\n";
	output << "# pID         Id of parent particle\n";
	output << "# pE          Energy [EeV] of parent particle\n";
	output << "# iID         Id of source particle\n";
	output << "# iE          Energy [EeV] of source particle\n";
	output << "#\n";
}

PhotonOutput1D::~PhotonOutput1D() {
}

void PhotonOutput1D::process(Candidate *candidate) const {
	int pid = candidate->current.getId();
	if ((pid != 22) and (abs(pid) != 11))
		return;

	char buffer[1024];
	size_t p = 0;

	p += sprintf(buffer + p, "%4i\t", pid);
	p += sprintf(buffer + p, "%g\t", candidate->current.getEnergy() / EeV);
	p += sprintf(buffer + p, "%8.4f\t", candidate->current.getPosition().getR() / Mpc);

	p += sprintf(buffer + p, "%10i\t", candidate->created.getId());
	p += sprintf(buffer + p, "%8.4f\t", candidate->created.getEnergy() / EeV);

	p += sprintf(buffer + p, "%10i\t", candidate->source.getId());
	p += sprintf(buffer + p, "%8.4f\n", candidate->source.getEnergy() / EeV);


#pragma omp critical
	{
		output.write(buffer, p);
	}

	candidate->setActive(false);
}

void PhotonOutput1D::close() {
	output.flush();
}

string PhotonOutput1D::getDescription() const {
	std::stringstream s;
	s << "PhotonOutput1D: Output file = " << filename;
	return s.str();
}

} // namespace crpropa
