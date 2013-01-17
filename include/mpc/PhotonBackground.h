#ifndef MPC_PHOTONBACKGROUND_H_
#define MPC_PHOTONBACKGROUND_H_

namespace mpc {

// Photon fields
enum PhotonField {
	CMB, IRB, CMB_IRB
};

// Returns overall photon field scaling factor at redshift z
double photonFieldScaling(int photonField, double z);

} // namespace mpc

#endif // MPC_PHOTONBACKGROUND_H_
