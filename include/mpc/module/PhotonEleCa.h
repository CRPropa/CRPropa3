#ifndef MPC_PHOTON_ELECA_H_
#define MPC_PHOTON_ELECA_H_

#include "mpc/Module.h"
#include "mpc/magneticField/MagneticField.h"

namespace mpc {

class PhotonEleCa: public Module {
private:
	ref_ptr<MagneticField> field;

public:
	PhotonEleCa(ref_ptr<MagneticField> field);
	~PhotonEleCa();
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

} // namespace mpc

#endif /* MPC_PHOTON_ELECA_H_ */
