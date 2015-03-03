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
		double lowerEnergyThreshold = 0.8010882435, // 5 EeV
		double magneticFieldStrength = 1E-13,				// 1 nG
    const std::string &background = "ALL");


/** Use Dint to calculate the observed sprectrum fom the photons recorded in inputfile.
 *  Magneticfield [T], default 1E-13 T = 1 nG
 *
 *  Implemented backgrounds are selected via the flags
 *		IRFlag = 0 (High IR), 1 (Low IR), 2 (Primack IR), 3 (Not defined), 4 (EleCa)
 *  RadioFlag = 0 (High Radio), 1 (Low Radio), 2 (Obs Radio), 3 (Null Radio), 4 (EleCa)
 * **/
void DintPropagation(const std::string &inputfile,
		const std::string &outputfile, 
		double magneticFieldStrength = 1E-13, // 1 nG
		int IRFlag = 4, int RadioFlag = 4,		
		double Zmax = 5);


/** Propagate high energy photons using eleca until the crossoverenergy is
 * reached. Then the cascade is continued using dint.
 **/
void DintElcaPropagation(const std::string &inputfile,
	const std::string &outputfile, 
	bool showProgress = true,
	double crossOverEnergy = 0.08010882435,  // in Joule! 5E17 eV = 0.080 J
	double magneticFieldStrength = 1E-13);

} // namespace crpropa

#endif // CRPROPA_PHOTON_PROPAGATION_H
