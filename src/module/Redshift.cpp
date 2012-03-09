#include "mpc/module/Redshift.h"
#include "mpc/Units.h"

namespace mpc {

SimpleRedshift::SimpleRedshift(Vector3 center) {
	this->center = center;
}

void SimpleRedshift::process(Candidate *candidate) const {
	double d = (candidate->current.getPosition() - center).mag();
	double z = 0.024 * d / 100 / Mpc;
	double dz = candidate->getRedshift() - z;

	// redshift can only decrease
	if (dz < 0)
		return;

	// using dE/dz=E/(1+z)
	candidate->current.setEnergy(
			candidate->current.getEnergy() * (1 - dz) / (1 + z));
	candidate->setRedshift(z);
}

std::string SimpleRedshift::getDescription() const {
	return "Simple redshift";
}

} // namespace mpc

// from CRPropa
//double omegaM = 0.3;
//double omegaL = 0.7;
//double H0 = 7e4 * meter / second;
//size_t nBins = 1000;
//double zMin = 0.0001;
//double zMax = 100;
//
//std::vector<double> z;
//std::vector<double> R;
//std::vector<double> D;
//z.push_back(0);
//R.push_back(1);
//D.push_back(0);
//
//for (size_t i=0; i<(nBins-1); i++) {
//	z.push_back( zMin * pow( zMax/zMin, i / double(nBins-2) ) );
//	R.push_back( 1. / ( (1.+z[i]) * sqrt( omegaL + omegaM * pow(1.+z[i],3) ) ) );
//	D.push_back( D[i-1] + (z[i] - z[i-1]) / 2. * (R[i-1] + R[i]) * c_light/H0 );
//}
