#ifndef ELECTRONPAIRPRODUCTION_H_
#define ELECTRONPAIRPRODUCTION_H_

#include "mpc/Module.h"
#include "mpc/ParticleState.h"

#include <vector>

namespace mpc {

class ElectronPairProduction: public Module {

public:
	enum PhotonField {
		CMB, IR, CMBIR
	};

	ElectronPairProduction(PhotonField photonField);
	ElectronPairProduction();
	void init(PhotonField photonField);
	void init(std::string filename);
	std::string getDescription() const;
	void process(Candidate *candidate, std::vector<Candidate *> &secondaries);

private:
	PhotonField photonField;
	std::vector<double> y; // energy loss rate table for protons in [J/m]
	std::vector<double> x; // energy table in [J]

};

} // namespace mpc

#endif /* ELECTRONPAIRPRODUCTION_H_ */
