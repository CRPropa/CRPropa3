#ifndef CRPROPA_PHOTON_OUTPUT_H
#define CRPROPA_PHOTON_OUTPUT_H

#include "crpropa/Module.h"

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
	void close();
};

} // namespace crpropa

#endif // CRPROPA_PHOTON_OUTPUT_H
