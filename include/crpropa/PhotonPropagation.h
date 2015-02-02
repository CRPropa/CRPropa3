#ifndef CRPROPA_PHOTON_PROPAGATION_H
#define CRPROPA_PHOTON_PROPAGATION_H

#include <string>
#include <vector>

namespace crpropa {


/** Use EleCa to propagate the photons recorded in inputfile.
 *  Propagation is stopped if the particles reach the observer or their energy
 *  dropps below lowerEnergyThreshold [J].
 * **/
void EleCaPropagation(const std::string &inputfile,
		const std::string &outputfile,
		bool showProgress=true,
		double lowerEnergyThreshold = 0.001602176487,
    const std::string &background = "ALL");


/** Use Dint to calculate the observed sprectrum fom the photons recorded in inputfile.
 * **/
void DintPropagation(const std::string &inputfile,
		const std::string &outputfile, int IRFlag = 2, int RadioFlag = 2,
		double Zmax = 5);

} // namespace crpropa

#endif // CRPROPA_PHOTON_PROPAGATION_H
