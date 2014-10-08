#ifndef CRPROPA_PHOTON_PROPAGATION_H
#define CRPROPA_PHOTON_PROPAGATION_H

#include <string>
#include <vector>

namespace crpropa {


void EleCaPropagation(const std::string &inputfile,
		const std::string &outputfile,
		bool showProgress=true,
		double lowerEnergyThreshold = 1E16,
    const std::string &background = "ALL");


void DintPropagation(const std::string &inputfile,
		const std::string &outputfile, int IRFlag = 2, int RadioFlag = 2,
		double Zmax = 5);

} // namespace crpropa

#endif // CRPROPA_PHOTON_PROPAGATION_H
