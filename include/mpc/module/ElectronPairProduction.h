#ifndef ELECTRONPAIRPRODUCTION_H_
#define ELECTRONPAIRPRODUCTION_H_

#include "mpc/Module.h"

#include <vector>
#include <string>

namespace mpc {

/**
 @class ElectronPairProduction
 @brief Electron-pair production of charged nuclei with background photons.

 This module simulates electron-pair production as a continuous energy loss.\n
 Several photon fields can be selected. They are considered as homogeneous and evolving as the CMB.
 */
class ElectronPairProduction: public Module {
protected:
	int photonField;
	std::vector<double> y; // energy loss rate table for protons in [J/m]
	std::vector<double> x; // energy table in [J]

public:
	ElectronPairProduction(int photonField);
	ElectronPairProduction();
	void init(int photonField);
	void init(std::string filename);
	std::string getDescription() const;
	void process(Candidate *candidate) const;
};

} // namespace mpc

#endif /* ELECTRONPAIRPRODUCTION_H_ */
