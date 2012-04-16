#ifndef SOPHIAPHOTOPIONPRODUCTION_H_
#define SOPHIAPHOTOPIONPRODUCTION_H_

#include "mpc/Module.h"
#include "mpc/Random.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

namespace mpc {

/**
 @class SophiaPhotoPionProduction
 @brief Photo-pion interactions of nuclei with background photons using SOPHIA.

 This module simulates photo-hadronic interactions of nuclei with background photons.\n
 Several photon fields can be selected. They are considered as homogeneous and evolving as the CMB.\n
 */
class SophiaPhotoPionProduction: public Module {
public:
	enum PhotonField {
		CMB, IR
	};

	SophiaPhotoPionProduction(PhotonField photonField);
	SophiaPhotoPionProduction();
	~SophiaPhotoPionProduction();
	void init(PhotonField photonField);
	void init(std::string filename);
	void process(Candidate *candidate) const;
	bool setNextInteraction(Candidate *candidate) const;
	void performInteraction(Candidate *candidate) const;
	std::string getDescription() const;

private:
	std::string name;
	PhotonField photonField;
	gsl_interp_accel *acc;
	gsl_spline *pRate; // interaction rate in [1/m] for protons
	gsl_spline *nRate; // interaction rate in [1/m] for neutrons
	double Emin, Emax;
};

} // namespace mpc

#endif /* SOPHIAPHOTOPIONPRODUCTION_H_ */
