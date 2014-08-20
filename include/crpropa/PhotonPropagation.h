#ifndef CRPROPA_PHOTON_PROPAGATION_H
#define CRPROPA_PHOTON_PROPAGATION_H

#include <string>
#include <vector>

namespace crpropa {

void EleCaPropagation(const std::string &inputfile,
		const std::string &background, std::vector<double> &energy,
		std::vector<double> &spectrum, 
		double lowerEnergyThreshold = 1E16
    );
void EleCaPropagation(const std::string &inputfile,
		const std::string &outputfile,
		double lowerEnergyThreshold = 1E16,
    const std::string &background = "ALL");
void DintPropagation(const std::string &inputfile,
		const std::string &outputfile, int IRFlag = 2, int RadioFlag = 2,
		double Zmax = 5);

} // namespace crpropa

#endif // CRPROPA_PHOTON_PROPAGATION_H
