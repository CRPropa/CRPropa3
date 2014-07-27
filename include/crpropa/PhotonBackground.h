#ifndef CRPROPA_PHOTONBACKGROUND_H
#define CRPROPA_PHOTONBACKGROUND_H

namespace crpropa {

// Photon fields
// The default IRB model is that of Kneiske et al. 2004
enum PhotonField {
	CMB, IRB, IRB_Kneiske04, IRB_Kneiske10, IRB_Franceschini08, IRB_Stecker05
};

// Returns overall comoving scaling factor
double photonFieldScaling(PhotonField photonField, double z);

} // namespace crpropa

#endif // CRPROPA_PHOTONBACKGROUND_H
