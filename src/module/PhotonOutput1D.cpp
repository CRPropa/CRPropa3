#include "crpropa/module/PhotonOutput1D.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

namespace crpropa {

PhotonOutput1D::PhotonOutput1D(const string &filename) :
		filename(filename), output(filename.c_str()) {
	output << "# sE\tsID\tpE\tpID\tE\tR\tz\n";
	output << "#\n";
	output << "# iE          Energy [EeV] of source particle\n";
	output << "# iID         Id of source particle\n";
	output << "# pE          Energy [EeV] of parent particle\n";
	output << "# pID         Id of parent particle\n";
	output << "# E           Energy [EeV]\n";
	output << "# R           Distance point of creation [Mpc]\n";
	output << "# z           Redshift\n";
	output << "#\n";
}

PhotonOutput1D::~PhotonOutput1D() {
}

void PhotonOutput1D::process(Candidate *candidate) const {
	if (candidate->current.getId() != 22)
		return;

	char buffer[1024];
	size_t p = 0;

	p += sprintf(buffer + p, "%8.4f\t", candidate->source.getEnergy() / EeV);
	p += sprintf(buffer + p, "%10i\t", candidate->source.getId());

	p += sprintf(buffer + p, "%8.4f\t", candidate->created.getEnergy() / EeV);
	p += sprintf(buffer + p, "%10i\t", candidate->created.getId());

	p += sprintf(buffer + p, "%8.4f\t", candidate->current.getEnergy() / EeV);
	p += sprintf(buffer + p, "%8.4f\t",
			candidate->current.getPosition().getR() / Mpc);
	p += sprintf(buffer + p, "%8.4f\n", candidate->getRedshift());

#pragma omp critical
	{
		output.write(buffer, p);
		//output.flush();
	}

	candidate->setActive(false);
}

string PhotonOutput1D::getDescription() const {
	std::stringstream s;
	s << "PhotonOutput1D: Output file = " << filename;
	return s.str();
}

} // namespace crpropa
