#ifndef CRPROPA_PHOTODISINTEGRATION_H
#define CRPROPA_PHOTODISINTEGRATION_H

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"

#include <vector>

namespace crpropa {

/**
 @class PhotoDisintegration
 @brief Photo-disintegration of nuclei with background photons.
 */
class PhotoDisintegration: public Module {
private:
	PhotonField photonField;
	double limit; // fraction of mean free path for limiting the next step
	struct PDMode {
		int channel; // number of emitted (n, p, H2, H3, He3, He4)
		std::vector<double> rate; // disintegration rate [1/m]
	};
	std::vector<std::vector<PDMode> > pdTable; // pdTable[Z * 31 + N] = vector<PDmode>
	double lgmin; // minimum log10(Lorentz-factor)
	double lgmax; // maximum log10(Lorentz-factor)

public:
	PhotoDisintegration(PhotonField photonField = CMB, double limit = 0.1);
	void setLimit(double l);
	void init(PhotonField photonField);
	void init(std::string filename);
	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate, int channel) const;
	double interactionRate(int id, double E, double z);
};

} // namespace crpropa

#endif // CRPROPA_PHOTODISINTEGRATION_H
