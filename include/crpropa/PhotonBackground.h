#ifndef CRPROPA_PHOTONBACKGROUND_H
#define CRPROPA_PHOTONBACKGROUND_H

namespace crpropa {

// Photon fields
enum PhotonField {
	CMB, IRB, CMB_IRB
};

// Returns overall photon field scaling factor at redshift z
double photonDensityScaling(int photonField, double z);

// Returns overall scaling factor for loss rates dE/dt at redshift z
double lossRateScaling(PhotonField photonField, double z);

} // namespace crpropa

#endif // CRPROPA_PHOTONBACKGROUND_H
