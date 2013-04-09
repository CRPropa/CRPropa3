#ifndef MPC_MODULE_PHOTON_OUTPUT_H_
#define MPC_MODULE_PHOTON_OUTPUT_H_

#include "mpc/Module.h"
#include "mpc/magneticField/MagneticField.h"

#include <fstream>

namespace mpc {

class PhotonOutput: public Module {
private:
	std::string filename, dataPath;
	mutable std::ofstream fout;
	ref_ptr<MagneticField> field;

	int IRFlag, RadioFlag;
	double H0, OmegaLambda, OmegaM;
	double Zmax, Cutcascade_Magfield;

public:
	PhotonOutput(const std::string &filename, ref_ptr<MagneticField> field);
	~PhotonOutput();
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

} // namespace mpc

#endif /* MPC_MODULE_PHOTON_OUTPUT_H_ */
