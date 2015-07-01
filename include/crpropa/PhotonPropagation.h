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
 *  IRFlag 0: High IR Background , 1: Low IR BAckground , 2: Primack Background, 4: Eleca Background (Stecker, Malkan, Scully (2006) arXiv:astro-ph/0510449v4)
 *  RadioFlag 0: High Radio Background , 1: MEd Radio BAckground , 2: Obs. Radio background , 3: No Radio, 4: Eleca Background (Protheroe, Bierman 1996.)
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
	double magneticFieldStrength = 1E-13);	 // 1 nG

} // namespace crpropa

#endif // CRPROPA_PHOTON_PROPAGATION_H
