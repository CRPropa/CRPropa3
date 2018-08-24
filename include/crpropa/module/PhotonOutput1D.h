#ifndef CRPROPA_PHOTON_OUTPUT_H
#define CRPROPA_PHOTON_OUTPUT_H

#include "crpropa/Module.h"

#include <fstream>

namespace crpropa {
/**
 * \addtogroup Output
 * @{
 */

class PhotonOutput1D: public Module {
private:
	std::ostream *out;
	std::string filename;
	mutable std::ofstream outfile;

public:
	PhotonOutput1D();
	PhotonOutput1D(std::ostream &out);
	PhotonOutput1D(const std::string &filename);
	~PhotonOutput1D();
	void process(Candidate *candidate) const;
	std::string getDescription() const;
	void close();
	void gzip();
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_PHOTON_OUTPUT_H
