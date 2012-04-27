#ifndef PHOTOPIONPRODUCTION_H_
#define PHOTOPIONPRODUCTION_H_

#include "mpc/module/StochasticInteraction.h"
#include "mpc/Random.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

namespace mpc {

/**
 @class PhotoPionProduction
 @brief Photo-pion interactions of nuclei with background photons.

 This module simulates photo-hadronic interactions of nuclei with background photons.\n
 Several photon fields can be selected. They are considered as homogeneous and evolving as the CMB.\n
 */
class PhotoPionProduction: public StochasticInteraction {
protected:
	int photonField;
	gsl_interp_accel *acc;
	gsl_spline *pRate; // interaction rate in [1/m] for protons
	gsl_spline *nRate; // interaction rate in [1/m] for neutrons
	double Emin, Emax;

public:
	PhotoPionProduction(int photonField = CMB);
	~PhotoPionProduction();
	void init(int photonField);
	void init(std::string filename);
	bool setNextInteraction(Candidate *candidate,
			InteractionState &interaction) const;
	void performInteraction(Candidate *candidate) const;
};

/**
 @class SophiaPhotoPionProduction
 @brief Photo-pion interactions of nuclei with background photons using SOPHIA.
 */
class SophiaPhotoPionProduction: public PhotoPionProduction {
protected:
	bool havePhotonsElectrons;
	bool haveNeutrinos;
	bool haveAntiNucleons;
public:
	SophiaPhotoPionProduction(int photonField = CMB, bool photonsElectrons =
			false, bool neutrinos = false, bool antiNucleons = false);
	void performInteraction(Candidate *candidate) const;
};

} // namespace mpc

#endif /* PHOTOPIONPRODUCTION_H_ */
