#ifndef CRPROPA_PHOTONBACKGROUND_H
#define CRPROPA_PHOTONBACKGROUND_H

#include <crpropa/Vector3.h>
#include <crpropa/Grid.h>

namespace crpropa {

/**
 * \addtogroup EnergyLosses
 * @{
 */
// Photon fields
// The default IRB model is that of Kneiske et al. 2004
// The slots PF1 to PF8 may be used for custom photon fields
enum PhotonField {
	CMB,
	IRB,  // same as IRB_Kneiske04
	IRB_Kneiske04,
	IRB_Stecker05,
	IRB_Franceschini08,
	IRB_Finke10,
	IRB_Dominguez11,
	IRB_Gilmore12,
	IRB_Stecker16_upper,
	IRB_Stecker16_lower,
	URB_Protheroe96,
	PF1, PF2, PF3, PF4,  // customizable
	PF5, PF6, PF7, PF8,  // field slots
};

// Returns overall comoving scaling factor
double photonFieldScaling(PhotonField photonField, double z);

// Returns a string representation of the field
std::string photonFieldName(PhotonField photonField);

/** 
 @class CustomPhotonField
 @brief Handler class for photon fields. Provides the sampleEps method.

 sampleEps draws a photon from a given photon background. This method 
 and all methods it depends on have been inspired by the SOPHIA code.
 */
class CustomPhotonField {
public:
	/** Constructor for photon field data
	 @param fieldPath  path/to/photonField.txt
	 */
	explicit CustomPhotonField(std::string fieldPath);

	/* Empty constructor to ease initialization in some modules
	 */
	CustomPhotonField();

	/** Draws a photon from the photon background
	 @param onProton  true=proton, false=neutron
	 @param Ein       energy of primary
	 @param zIn       redshift of primary
	 */
	double sampleEps(bool onProton, double Ein, double zIn) const;

	/** Returns the photon field density in 1/(Jm³).
		Multiply by h*nu for physical photon density.
	 @param eps   photon energy in eV
	 @param zIn   redshift
	*/
	double getPhotonDensity(double eps, double z) const;

	/** Returns the crossection of p-gamma interaction
	 @param eps       photon energy
	 @param onProton  true=proton, false=neutron
	*/
	double SOPHIA_crossection(double eps, bool onProton) const;

	std::vector<double> photonEnergy;
	std::vector<double> photonRedshift;
	std::vector<double> photonDensity;
protected:
	void init(std::string fieldPath);
	double SOPHIA_probEps(double eps, bool onProton, double Ein, double zIn) const;
	double SOPHIA_pl(double x, double xth, double xmax, double alpha) const;
	double SOPHIA_ef(double x, double th, double w) const;
	double SOPHIA_breitwigner(double sigma_0, double Gamma, double DMM, double epsPrime, bool onProton) const;
	double SOPHIA_functs(double s, bool onProton) const;
};

/** @}*/
} // namespace crpropa

#endif // CRPROPA_PHOTONBACKGROUND_H
