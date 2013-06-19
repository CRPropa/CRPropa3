#ifndef MPC_MODULE_PHOTON_DINT_H_
#define MPC_MODULE_PHOTON_DINT_H_

#include "mpc/Module.h"
#include "mpc/magneticField/MagneticField.h"

#include <fstream>

namespace mpc {

class PhotonDINT: public Module {
private:
	std::string filename, dataPath;
	mutable std::ofstream fout;
	ref_ptr<MagneticField> field;

	int IRFlag, RadioFlag;
	double Zmax, Cutcascade_Magfield;

public:
	PhotonDINT(const std::string &filename, ref_ptr<MagneticField> field);
	~PhotonDINT();
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

} // namespace mpc

#endif /* MPC_MODULE_PHOTON_DINT_H_ */
