#ifndef PHOTODISINTEGRATION_H_
#define PHOTODISINTEGRATION_H_

#include "mpc/module/StochasticInteraction.h"
#include "mpc/Random.h"

#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

namespace mpc {

/**
 @class PhotoDisintegration
 @brief Photo-disintegration of nuclei with background photons.

 This module simulates photo-disintegration of nuclei with background photons.\n
 Background photon fields are considered as homogeneous and evolving as the CMB.\n
 */
class PhotoDisintegration: public StochasticInteraction {
private:
	int photonField;
	struct PDMode {
		int channel; // number of emitted (n, p, H2, H3, He3, He4)
		gsl_spline *rate; // disintegration rate [1/m]
	};
	std::vector<std::vector<PDMode> > pdTable; // pdTable[Z * 31 + N] = vector<PDmode>
	gsl_interp_accel *acc; // accelerator for the interpolation

public:
	PhotoDisintegration(int photonField = CMB_IRB);
	~PhotoDisintegration();
	void init(int photonField);
	void init(std::string filename);
	bool setNextInteraction(Candidate *candidate,
			InteractionState &interaction) const;
	void performInteraction(Candidate *candidate) const;
};

} // namespace mpc

#endif /* PHOTODISINTEGRATION_H_ */
