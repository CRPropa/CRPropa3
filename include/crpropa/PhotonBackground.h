#ifndef CRPROPA_PHOTONBACKGROUND_H
#define CRPROPA_PHOTONBACKGROUND_H

namespace crpropa {

// Photon fields
enum PhotonField {
	CMB, IRB
};

// Returns overall comoving scaling factor
double photonFieldScaling(PhotonField photonField, double z);

} // namespace crpropa

#endif // CRPROPA_PHOTONBACKGROUND_H
