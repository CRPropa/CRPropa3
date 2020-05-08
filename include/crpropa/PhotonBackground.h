#ifndef CRPROPA_PHOTONBACKGROUND_H
#define CRPROPA_PHOTONBACKGROUND_H

#include "crpropa/Common.h"
#include "crpropa/Referenced.h"

#include <string>

namespace crpropa {
/**
 * \addtogroup PhotonFields
 * @{
 */

/**
 @class PhotonField
 @brief Abstract base class for photon fields.
 */
class PhotonField: public Referenced {
public:
	PhotonField() {
		this->fieldName = "AbstractPhotonField";
		this->isRedshiftDependent = false;
	}

	/**
	 returns comoving photon density [1/m^3].
	 multiply with (1+z^3) for physical number density.
	 @param ePhoton		photon energy [J]
	 @param z			redshift
	 */
	virtual double getPhotonDensity(double ePhoton, double z = 0.) const {
		return 0.;
	};

	/**
	 returns overall comoving scaling factor
	 (cf. CRPropa3-data/calc_scaling.py)
	 @param z		redshift
	 */
	virtual double getRedshiftScaling(double z) const {
		return 1.;
	};

	bool hasRedshiftDependence() const {
		return this->isRedshiftDependent;
	}

	std::string getFieldName() const {
		return this->fieldName;
	}

protected:
	std::string fieldName;
	bool isRedshiftDependent;
};

/**
 @class TabularPhotonField
 @brief Photon field decorator for tabulated photon fields.
 */
class TabularPhotonField: public PhotonField {
public:
	TabularPhotonField(std::string fieldName, bool hasRedshiftDependence = true);
	double getPhotonDensity(double ePhoton, double z = 0.) const;
	double getRedshiftScaling(double z) const;
	bool hasRedshiftDependence() const;

protected:
	void initPhotonEnergy(std::string fieldName);
	void initPhotonDensity(std::string fieldName);
	void initRedshift(std::string fieldName);
	void initRedshiftScaling();
	void checkInputData() const;

	std::vector<double> photonEnergies;
	std::vector<double> photonDensity;
	std::vector<double> redshifts;
	std::vector<double> redshiftScalings;
};

/**
 @class BlackbodyPhotonField
 @brief Photon field decorator for black body photon fields.
 */
class BlackbodyPhotonField: public PhotonField {
public:
	BlackbodyPhotonField(std::string fieldName, double blackbodyTemperature);
	double getPhotonDensity(double ePhoton) const;

protected:
	double blackbodyTemperature;
};

/**
 @class PhotonFieldSampling
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
