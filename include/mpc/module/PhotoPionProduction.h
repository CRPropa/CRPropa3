#ifndef PHOTODISINTEGRATION_H_
#define PHOTODISINTEGRATION_H_

#include "mpc/Module.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>

namespace mpc {

class PhotoPionProduction: public Module {
public:
	enum PhotonField {
		CMB, IR, CMBIR
	};

	PhotoPionProduction(PhotonField photonField);
	PhotoPionProduction();
	~PhotoPionProduction();
	void init(PhotonField photonField);
	void init(std::string filename);
	std::string getDescription() const;
	void process(Candidate *candidate, std::vector<Candidate *> &secondaries);
	void setNextInteraction(Candidate *candidate);
	void performInteraction(Candidate *candidate,
			std::vector<Candidate *> &secondaries);

private:
	PhotonField photonField;
	gsl_rng *rand;
	gsl_interp_accel *acc;
	gsl_spline *protonFreePath;
	gsl_spline *neutronFreePath;
	int cached_id;
	int cached_interaction;
	double cached_distance;
};

} // namespace mpc

#endif /* PHOTODISINTEGRATION_H_ */
