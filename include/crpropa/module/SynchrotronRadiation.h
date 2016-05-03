#ifndef CRPROPA_SYNCHROTRONRADIATION_H
#define CRPROPA_SYNCHROTRONRADIATION_H

#include "crpropa/Module.h"
#include "crpropa/magneticField/MagneticField.h"

namespace crpropa {

/**
 @class SynchrotronRadiation
 @brief Energy loss by synchrotron photons which are radiated off by charged particles deflected in magnetic fields.

 This module simulates the production of synchrotron photons as a continuous energy loss.\n
 The synchrotron spectrum is calculated from J.D. Jackson p. 785 (14.91).\n
 The radiated power is taken from J.D. Jackson p. 770 (14.31).\n
 The production of secondary photons can be activated.
 */
class SynchrotronRadiation: public Module {
private:
	ref_ptr<MagneticField> field;
	double Brms;

	std::vector<double> tabx; /*< tabulated fraction of omega_{photon} to omega_{critical} from 1e-6 to 1e2 in 801 steps logspaced*/
	std::vector<double> tabCDF; /*< tabulated cumulative synchrotron spectrum*/

	double limit; ///< fraction of energy loss length to limit the next step
	bool havePhotons;

public:
	SynchrotronRadiation(ref_ptr<MagneticField> field, bool havePhotons =
			false, double limit = 0.1);
	SynchrotronRadiation(double Brms = 0., bool havePhotons =
			false, double limit = 0.1);

	void setField(ref_ptr<MagneticField> field);
	void setField(double Brms);
	ref_ptr<MagneticField> getField();
	double getBrms();
	void setHavePhotons(bool havePhotons);
	void setLimit(double limit);
	void initSpectrum();
	void process(Candidate *candidate) const;
	void addPhotons(Candidate *candidate, double loss) const;
	std::string getDescription() const;
};

} // namespace crpropa

#endif // CRPROPA_SYNCHROTRONRADIATION_H
