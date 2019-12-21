#ifndef CRPROPA_PHOTONBACKGROUND_H
#define CRPROPA_PHOTONBACKGROUND_H

#include <string>

namespace crpropa {

/**
 * \addtogroup EnergyLosses
 * @{
 */
// Photon fields
// The default IRB model is that of Kneiske et al. 2004
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
	URB_Protheroe96
};

// Returns overall comoving scaling factor
double photonFieldScaling(PhotonField photonField, double z);

// Returns a string representation of the field
std::string photonFieldName(PhotonField photonField);


/**
 @class photonFieldSampling
 @brief Reimplementation of SOPHIA photon sampling. Naming and unit conventions are taken from SOPHIA to ease comparisions.
 */
class PhotonFieldSampling {
public:
	PhotonFieldSampling();

	/**
	 Constructor to mimic SOPHIA structure.
	  @param bgFlag		1: CMB | 2: IRB_Kneiske04
	 */
	explicit PhotonFieldSampling(int bgFlag);

	/**
	 SOPHIA's photon sampling method. Returns energy [J] of a photon of the photon field.
	 @param onProton	particle type: proton or neutron
	 @param E_in		energy of incoming nucleon
	 @param z_in		redshift of incoming nucleon
	 */
	double sample_eps(bool onProton, double E_in, double z_in) const;
protected:
	int bgFlag;

	// called by: sample_eps
	// - input: photon energy [eV], redshift
	// - output: photon density per unit energy [#/(eVcm^3)]
	double getPhotonDensity(double eps, double z_in) const;

	// called by: sample_eps, prob_eps
	// - input:  type: specifier of function over which to integrate, integration limits A and B, arguents for functions over which to integrate
	// - output: 8-points Gauß-Legendre integral
	double gaussInt(std::string type, double lowerLimit, double upperLimit, bool onProton, double E_in, double z_in) const;

	// called by: sample_eps  
	// - input: s [GeV^2] 
	// - output: (s-p^2) * sigma_(nucleon/gamma) [GeV^2 * mubarn]
	double functs(double s, bool onProton) const;
	
	// called by: sample_eps, gaussInt
	// - input: photon energy eps [eV], E_in [GeV]
	// - output: probability to encounter photon of energy eps
	double prob_eps(double eps, bool onProton, double E_in, double z_in) const;

	// called by: functs
	// - input: photon energy [eV]
	// - output: crossection of nucleon-photon-interaction [mubarn]
	double crossection(double eps, bool onProton) const;

	// called by: crossection
	// - input: photon energy [eV], threshold [eV], max [eV], unknown [no unit]
	// - output: unknown [no unit]
	double Pl(double x, double xth, double xmax, double alpha) const;

	// called by: crossection
	// - input: photon energy [eV], threshold [eV], unknown [eV]
	// - output: unknown [no unit]
	double Ef(double x, double th, double w) const;

	// called by: crossection
	// - input: cross section [µbarn], width [GeV], mass [GeV/c^2], rest frame photon energy [GeV]
	// - output: Breit-Wigner crossection of a resonance of width Gamma
	double breitwigner(double sigma_0, double Gamma, double DMM, double epsPrime, bool onProton) const;
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_PHOTONBACKGROUND_H
