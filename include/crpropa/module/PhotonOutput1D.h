#ifndef CRPROPA_PHOTON_OUTPUT_H
#define CRPROPA_PHOTON_OUTPUT_H

#include "crpropa/Module.h"

#include <memory>
#include <fstream>

namespace crpropa {

class PhotonOutput1D: public Module {
private:
	std::string filename;
	mutable std::ofstream output;
public:
	PhotonOutput1D(const std::string &filename);
	~PhotonOutput1D();
	void process(Candidate *candidate) const;
	std::string getDescription() const;
	void endRun();
};

} // namespace crpropa

#endif // CRPROPA_PHOTON_OUTPUT_H
