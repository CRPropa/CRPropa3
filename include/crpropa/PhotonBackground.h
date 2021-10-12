#ifndef CRPROPA_PHOTONBACKGROUND_H
#define CRPROPA_PHOTONBACKGROUND_H

#include "crpropa/Common.h"
#include "crpropa/Referenced.h"

#include <vector>
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
	 @param z			redshift (if redshift dependent, default = 0.)
	 */
	virtual double getPhotonDensity(double ePhoton, double z = 0.) const = 0;
	virtual double getMinimumPhotonEnergy(double z) const = 0;
	virtual double getMaximumPhotonEnergy(double z) const = 0;

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

 This class reads photon field data from files;
 The first file must be a list of photon energies [J], named fieldName_photonEnergy.txt
 The second file must be a list of comoving photon field densities [1/m^3], named fieldName_photonDensity.txt
 Optionally, a third file contains redshifts, named fieldName_redshift.txt
 */
class TabularPhotonField: public PhotonField {
public:
	TabularPhotonField(const std::string fieldName, const bool isRedshiftDependent = true);
	
	double getPhotonDensity(double ePhoton, double z = 0.) const;
	double getRedshiftScaling(double z) const;
	double getMinimumPhotonEnergy(double z) const;
	double getMaximumPhotonEnergy(double z) const;

protected:
	void readPhotonEnergy(std::string filePath);
	void readPhotonDensity(std::string filePath);
	void readRedshift(std::string filePath);
	void initRedshiftScaling();
	void checkInputData() const;

	std::vector<double> photonEnergies;
	std::vector<double> photonDensity;
	std::vector<double> redshifts;
	std::vector<double> redshiftScalings;
};

/**
 @class IRB_Kneiske04
 @brief Extragalactic background light model from Kneiske et al. 2004

 Source info:
 DOI:10.1051/0004-6361:20031542,
 https://www.aanda.org/articles/aa/pdf/2004/03/aa3848.pdf, figure 1 ("Best-fit" model)
 */
class IRB_Kneiske04: public TabularPhotonField {
public:
	IRB_Kneiske04() : TabularPhotonField("IRB_Kneiske04", true) {}
};

/**
 @class IRB_Stecker05
 @brief Extragalactic background light model by Stecker at al. 2005

 Source info:
 DOI:10.1086/506188, astro-ph/0510449
 https://iopscience.iop.org/article/10.1086/506188/pdf
 */
class IRB_Stecker05: public TabularPhotonField {
public:
	IRB_Stecker05() : TabularPhotonField("IRB_Stecker05", true) {}
};

/**
 @class IRB_Franceschini08
 @brief Extragalactic background light model from Franceschini et al. 2008

 Source info:
 DOI:10.1051/0004-6361:200809691
 https://arxiv.org/pdf/0805.1841.pdf, tables 1 and 2
 */
class IRB_Franceschini08: public TabularPhotonField {
public:
	IRB_Franceschini08() : TabularPhotonField("IRB_Franceschini08", true) {}
};

/**
 @class IRB_Finke10
 @brief Extragalactic background light model from Finke et al. 2010

 Source info:
 DOI:10.1088/0004-637X/712/1/238
 https://iopscience.iop.org/article/10.1088/0004-637X/712/1/238/pdf
 */
class IRB_Finke10: public TabularPhotonField {
public:
	IRB_Finke10() : TabularPhotonField("IRB_Finke10", true) {}
};

/**
 @class IRB_Dominguez11
 @brief Extragalactic background light model from Dominguez et al. 2011

 Source info:
 DOI:10.1111/j.1365-2966.2010.17631.x
 https://academic.oup.com/mnras/article/410/4/2556/1008012
 */
class IRB_Dominguez11: public TabularPhotonField {
public:
	IRB_Dominguez11() : TabularPhotonField("IRB_Dominguez11", true) {}
};

/**
 @class IRB_Gilmore12
 @brief Extragalactic background light model from Gilmore et al. 2012

 Source info:
 DOI:10.1111/j.1365-2966.2012.20841.x
 https://academic.oup.com/mnras/article/422/4/3189/1050758
 */
class IRB_Gilmore12: public TabularPhotonField {
public:
	IRB_Gilmore12() : TabularPhotonField("IRB_Gilmore12", true) {}
};

/**
 @class IRB_Stecker16_upper
 @brief Extragalactic background light model from Stecker et al. 2016 (upper-bound model)

 Source info:
 DOI:10.3847/0004-637X/827/1/6
 https://iopscience.iop.org/article/10.3847/0004-637X/827/1/6
 */
class IRB_Stecker16_upper: public TabularPhotonField {
public:
	IRB_Stecker16_upper() : TabularPhotonField("IRB_Stecker16_upper", true) {}
};

/**
 @class IRB_Stecker16_lower
 @brief Extragalactic background light model from Stecker et al. 2016 (lower-bound model)

 Source info:
 DOI:10.3847/0004-637X/827/1/6
 https://iopscience.iop.org/article/10.3847/0004-637X/827/1/6
 */
class IRB_Stecker16_lower: public TabularPhotonField {
public:
	IRB_Stecker16_lower() : TabularPhotonField("IRB_Stecker16_lower", true) {}
};

/**
 @class URB
 @brief Extragalactic background light model from Protheroe & Biermann 1996

 Source info:
 DOI:10.1016/S0927-6505(96)00041-2
 https://www.sciencedirect.com/science/article/abs/pii/S0927650596000412
 */
class URB_Protheroe96: public TabularPhotonField {
public:
	URB_Protheroe96() : TabularPhotonField("URB_Protheroe96", false) {}
};

/**
 @class URB
 @brief Extragalactic background light model based on ARCADE2 observations, by Fixsen et al.
 Note that this model does not cover the same energy range as other URB models. Here, only ~10 MHz - 10 GHz is considered. 
 Therefore, it only makes sense to use this model in very specific studies.
  
 Source info:
 DOI:10.1088/0004-637X/734/1/5 
 https://iopscience.iop.org/article/10.1088/0004-637X/734/1/5
 */
class URB_Fixsen11: public TabularPhotonField {
public:
	URB_Fixsen11() : TabularPhotonField("URB_Fixsen11", false) {}
};

/**
 @class URB
 @brief Extragalactic background light model by Nitu et al. 
  
 Source info:
 DOI:10.1016/j.astropartphys.2020.102532
 https://www.sciencedirect.com/science/article/pii/S0927650520301043?
 */
class URB_Nitu21: public TabularPhotonField {
public:
	URB_Nitu21() : TabularPhotonField("URB_Nitu21", false) {}
};

/**
 @class BlackbodyPhotonField
 @brief Photon field decorator for black body photon fields.
 */
class BlackbodyPhotonField: public PhotonField {
public:
	BlackbodyPhotonField(const std::string fieldName, const double blackbodyTemperature);
	double getPhotonDensity(double ePhoton, double z = 0.) const;
	double getMinimumPhotonEnergy(double z) const;
	double getMaximumPhotonEnergy(double z) const;
	void setQuantile(double q);

protected:
	double blackbodyTemperature;
	double quantile;
};

/**
 @class CMB
 @brief Cosmic mircowave background photon field

 Source info:
 This field is an isotropic blackbody photon field with temperature T = 2.73 K
 */
class CMB: public BlackbodyPhotonField {
public:
	CMB() : BlackbodyPhotonField("CMB", 2.73) {}
};


/**
 @class PhotonFieldSampling
 @brief Reimplementation of SOPHIA photon sampling. Naming and unit conventions are taken from SOPHIA to ease comparisions.
 */
class PhotonFieldSampling :public Referenced {
public:
	PhotonFieldSampling();
	explicit PhotonFieldSampling(ref_ptr<PhotonField> photonField);

	/**
	 SOPHIA's photon sampling method. Returns energy [J] of a photon of the photon field.
	 @param onProton	particle type: proton or neutron
	 @param E		energy of incoming nucleon [J]
	 @param z		redshift of incoming nucleon
	 */
	double sampleEps(bool onProton, double E, double z) const;
	
	// using log or lin spacing of photons in the range between epsMin and
	// epsMax for computing the maximum probability of photons in field
	void setSampleLog(bool log);

	// given the discrete steps to compute the maximum interaction probability pEpsMax 
	// of photons in field, the real pEpsMax may lie between the descrete tested photon energies.
	// A correction factor can be set to increase pEpsMax by that factor
	void setCorrectionFactor(double factor);

	bool getSampleLog() const;
	double getCorrectionFactor() const;

protected:
	// called by: sampleEps
	// - input: s [GeV^2]
	// - output: (s-p^2) * sigma_(nucleon/gamma) [GeV^2 * mubarn]
	double functs(double s, bool onProton) const;

	// called by: sampleEps, gaussInt
	// - input: photon energy eps [eV], Ein [GeV]
	// - output: probability to encounter photon of energy eps
	double probEps(double eps, bool onProton, double Ein, double z) const;

	/** called by: sampleEps
	@param onProton	particle type: proton or neutron
	@param Ein		energy of incoming nucleon
	@param z		redshift of incoming nucleon
	@param epsMin   minimum photon energy of field
	@param epsMax   maximum photon energy of field
	- output: maximum probability of all photons in field
	 */
	double probEpsMax(bool onProton, double Ein, double z, double epsMin, double epsMax) const;

	/** called by: sampleEps
	@param onProton	particle type: proton or neutron
	@param Ein		energy of incoming nucleon
	- output: labframe energy [eV] of least energetic photon where PPP can occur
	 */
	double epsMinInteraction(bool onProton, double Ein) const;

	/** called by: probEps, epsMinInteraction
	@param onProton	particle type: proton or neutron
	@param Ein		energy of incoming nucleon
	- output: hadron momentum [GeV/c]
	 */
	double momentum(bool onProton, double Ein) const;
	
	// called by: functs
	// - input: photon energy [eV]
	// - output: crossection of nucleon-photon-interaction [mubarn]
	double crossection(double eps, bool onProton) const;

	// called by: crossection
	// - input: photon energy [eV], threshold [eV], max [eV], unknown [no unit]
	// - output: unknown [no unit]
	double Pl(double eps, double xth, double xMax, double alpha) const;

	// called by: crossection
	// - input: photon energy [eV], threshold [eV], unknown [eV]
	// - output: unknown [no unit]
	double Ef(double eps, double epsTh, double w) const;

	// called by: crossection
	// - input: cross section [µbarn], width [GeV], mass [GeV/c^2], rest frame photon energy [GeV]
	// - output: Breit-Wigner crossection of a resonance of width Gamma
	double breitwigner(double sigma0, double gamma, double DMM, double epsPrime, bool onProton) const;

	// called by: probEps, crossection, breitwigner, functs
	// - input: is proton [bool]
	// - output: mass [Gev/c^2]
	double mass(bool onProton) const;

	// - output: [GeV^2] head-on collision 
	double sMin() const;

	ref_ptr<PhotonField> photonField;
	bool sampleLog = true;
	double correctionFactor = 1.6; // increeses the maximum of the propability function
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_PHOTONBACKGROUND_H
