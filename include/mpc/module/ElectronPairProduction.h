#ifndef ELECTRONPAIRPRODUCTION_H_
#define ELECTRONPAIRPRODUCTION_H_

#include "mpc/Module.h"
#include "mpc/Common.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

namespace mpc {

/**
 @class ElectronPairProduction
 @brief Electron-pair production of charged nuclei with background photons.

 This module simulates electron-pair production as a continuous energy loss.\n
 Several photon fields can be selected. They are considered as homogeneous and evolving as the CMB.
 */
class ElectronPairProduction: public Module {
private:
	int photonField;
	gsl_interp_accel *acc;
	gsl_spline *lossRate; // energy loss rate in [J/m]
	double xMin, xMax, yMax;

public:
	ElectronPairProduction(int photonField = CMB_IRB);
	void init(int photonField);
	void init(std::string filename);
	void process(Candidate *candidate) const;
};

} // namespace mpc

#endif /* ELECTRONPAIRPRODUCTION_H_ */
