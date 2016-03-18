#ifndef CRPROPA_PHOTON_OUTPUT_THRESHOLD_1D_H
#define CRPROPA_PHOTON_OUTPUT_THRESHOLD_1D_H

#include "crpropa/Module.h"
#include "crpropa/Units.h"

#include <memory>
#include <fstream>

namespace crpropa {

class PhotonOutputThreshold1D: public Module {
private:
	std::string filename;
	mutable std::ofstream output;
	mutable double Ethreshold;
public:
	PhotonOutputThreshold1D(const std::string &filename);
	PhotonOutputThreshold1D(const std::string &filename, const double Ethr);
	~PhotonOutputThreshold1D();
	void process(Candidate *candidate) const;
	std::string getDescription() const;
	void endRun();
};

} // namespace crpropa

#endif // CRPROPA_PHOTON_OUTPUT_THRESHOLD_1D_H
