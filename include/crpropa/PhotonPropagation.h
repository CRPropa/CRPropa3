#ifndef CRPROPA_PHOTON_PROPAGATION_H
#define CRPROPA_PHOTON_PROPAGATION_H

#include <string>
#include <vector>

namespace crpropa {

void EleCaPropagation(const std::string &background,
		const std::string &inputfile, std::vector<double> &energy,
		std::vector<double> &spectrum);
void EleCaPropagation(const std::string &background,
		const std::string &inputfile, const std::string &outputfile);

} // namespace crpropa

#endif // CRPROPA_PHOTON_PROPAGATION_H
