#ifndef CRPROPA_PHOTON_PROPAGATION_H
#define CRPROPA_PHOTON_PROPAGATION_H

#include <string>
#include <vector>

/**
 @file
 @brief Calculation of the electromagnetic cascade outside of CRPropa's module list
 */

namespace crpropa {

/**
 Propagate photons, electrons and positrons using the EleCa code.
 The propagation is stopped when the particles reach the observer or their energy drops below the threshold energy.
 */
void ElecaPropagation(
	const std::string &inputfile,               //!< input in PhotonOutput1D format
	const std::string &outputfile,              //!< output in PhotonOutput1D format
	bool showProgress = true,                   //!< show a progress bar
	double lowerEnergyThreshold = 0.8010882435, //!< threshold energy [J], default = 5 EeV
	double magneticFieldStrength = 1E-13,       //!< magnetic field strength [T], default = 1 nG
	const std::string &background = "ALL"       //!< photon background string
	);

/**
 Calculate the electromagnetic cascade with DINT
 */
void DintPropagation(
	const std::string &inputfile,         //!< input in PhotonOutput1D format
	const std::string &outputfile,        //!< output spectrum (photons, electrons, positrons)
	int IRFlag = 4,                       //!< EBL background 0: high, 1: low, 2: Primack, 4: Stecker'06
	int RadioFlag = 4,                    //!< radio background 0: high, 1: medium, 2: obs, 3: none, 4: Protheroe'96
	double magneticFieldStrength = 1E-13, //!< magnetic field strength [T], default = 1 nG
	double aCutcascade_Magfield = 0       //!< a-parameter, see CRPropa 2 paper
	);

/**
 Propagate photons using EleCa for energies above the crossover energy and DINT below
 */
void DintElecaPropagation(
	const std::string &inputfile,           //!< input in PhotonOutput1D format
	const std::string &outputfile,          //!< output spectrum (photons, electrons, positrons)
	bool showProgress = true,               //!< show a progress bar
	double crossOverEnergy = 0.08010882435, //!< crossover energy [J] between EleCa and DINT, default = 0.5 EeV
	double magneticFieldStrength = 1E-13,   //!< magnetic field strength [T], default = 1 nG
	double aCutcascade_Magfield = 0         //!< a-parameter, see CRPropa 2 paper
	);

} // namespace crpropa

#endif // CRPROPA_PHOTON_PROPAGATION_H
