#ifndef CRPROPA_PHOTONBACKGROUND_H
#define CRPROPA_PHOTONBACKGROUND_H

#include "crpropa/Common.h"
#include "crpropa/Referenced.h"
#include "crpropa/Vector3.h"
#include "crpropa/Geometry.h"

#include <vector>
#include <string>
#include <unordered_map>

namespace crpropa {
/**
 * \addtogroup PhotonFields
 * @{
 */

/**
 @class PhotonField
 @brief Abstract base class for photon fields.
 */
class PhotonField {
public:
	
	/** Constructor
	 * Sets the following values per default:
	 * @var fieldName = "AbstractPhotonField";
	 * @var isRedshiftDependent = false;
	 * @var isPositionDependent = false;
	 * @var surface = nullptr;
	 */
	PhotonField() {
		this->fieldName = "AbstractPhotonField";
		this->isRedshiftDependent = false;
		this->isPositionDependent = false;
		this->surface = nullptr;
	}
	virtual ~PhotonField() = default; 

	/**
	 returns comoving photon density [1/m^3].
	 multiply with (1+z^3) for physical number density.
	 @param ePhoton		photon energy [J]
	 @param z			redshift (if redshift dependent, default = 0.)
	 @param pos			position (if position dependent, default = Vector3d(0,0,0))
	 */
	virtual double getPhotonDensity(double ePhoton, double z = 0., const Vector3d &pos = Vector3d(0.,0.,0.)) const = 0;
	/** Returns the minimum photon energy
	 * @param z  Redshift
	 * @param pos  Position 
	 */
	virtual double getMinimumPhotonEnergy(double z, const Vector3d &pos = Vector3d(0.,0.,0.)) const = 0;
	/** Returns the maximum photon energy
	 * @param z  Redshift
	 * @param pos  Position 
	 */
	virtual double getMaximumPhotonEnergy(double z, const Vector3d &pos = Vector3d(0.,0.,0.)) const = 0;
	/// Returns name of the current used photon field
	virtual inline std::string getFieldName() const {
		return this->fieldName;
	}
	
	/**
	 returns overall comoving scaling factor
	 (cf. CRPropa3-data/calc_scaling.py)
	 @param z		redshift
	 */
	virtual inline double getRedshiftScaling(double z) const {
		return 1.;
	};
	
	/// Returns if the derived field has a redshift dependency or not
	inline bool hasRedshiftDependence() const {
		return this->isRedshiftDependent;
	}
	/// Returns if the derived field has a position dependency or not
	inline bool hasPositionDependence() const {
		return this->isPositionDependent;
	}
	/// Returns the set surface
	inline ref_ptr<Surface> getSurface() const {
		return this->surface;
	}
	/** Sets the name of the currently used photon field
	 * @param fieldName  Name of the currently used photon field
	 */
	void setFieldName(std::string fieldName) {
		this->fieldName = fieldName;
	}
	
protected:
	std::string fieldName;  /**< Name of the currently used field */
	bool isRedshiftDependent;  /**< If photon field is redshift dependent */
	bool isPositionDependent;  /**< If photon field is position dependent */
	ref_ptr<Surface> surface;  /**< Currently used Surface */
	
};

/**
 @class TabularPhotonField
 @brief Photon field decorator for tabulated photon fields.

 This class reads photon field data from files;
 The first file must be a list of photon energies [J], named fieldName_photonEnergy.txt
 The second file must be a list of comoving photon field densities [1/m^3], named fieldName_photonDensity.txt
 Optionally, a third file contains redshifts, named fieldName_redshift.txt.
 */
class TabularPhotonField: public PhotonField {
	public:
		/** Constructor
		 * @param fieldName  String field identifier
		 * @param isRedshiftDependent  Whether or not the given field is redshift dependent
		 */
		TabularPhotonField(const std::string fieldName, const bool isRedshiftDependent = true);

		/** Returns photon density dependend on energy, redshift and position
		 * Returns the photon density for a specific photon energy, redshift and position
		 * @param ePhoton  Photon energy
		 * @param z  Redshift
		 * @param pos  Position
		 */
		double getPhotonDensity(double ePhoton, double z = 0., const Vector3d &pos = Vector3d(0.,0.,0.)) const;
		double getRedshiftScaling(double z) const;
		/** Returns the minimum possible photon energy
		 * Returns the minimum possible photon energy for a given redshift and position.
		 * @param z  Redshift
		 * @param pos Position
		 */
		double getMinimumPhotonEnergy(double z, const Vector3d &pos = Vector3d(0.,0.,0.)) const;
		/** Returns the maximum possible photon energy
		 * Returns the maximum possible photon energy for a given redshift and position.
		 * @param z  Redshift
		 * @param pos Position
		 */
		double getMaximumPhotonEnergy(double z, const Vector3d &pos = Vector3d(0.,0.,0.)) const;

	protected:
		/** Reads the photon energies from the given folder
		 * @param filePath  Path containing photon energy files
		 */
		void readPhotonEnergy(std::string filePath);
		/** Reads the photon densities from the given folder
		 * @param filePath  Path containing photon densitiy files
		 */
		void readPhotonDensity(std::string filePath);
		/** Reads the photon redshift from the given folder
		 * @param filePath  Path containing photon redshift files
		 */
		void readRedshift(std::string filePath);
		/// Initializes the redshift scaling from the read photon redshift
		void initRedshiftScaling();
		/// Checks if read data is valid
		void checkInputData() const;

		std::vector<double> photonEnergies;
		std::vector<double> photonDensity;
		std::vector<double> redshifts;
		std::vector<double> redshiftScalings;
};

/**
 @class TabularSpatialPhotonField
 @brief Position dependent photon field decorator for tabulated photon fields.

 This class reads photon field data from files in the appropriate directory;
 The first files must be lists of photon energies [J], named fieldName_photonEnergy.txt and contained in the subdirectory /photonEnegy/;
 The second files must be lists of comoving photon field densities [1/m^3], named fieldName_photonDensity.txt and contained in the subdirectory /photonDensity/;
 The generated files through the CRPropa procedure (https://crpropa.github.io/CRPropa3/pages/example_notebooks/custom_photonfield/custom-photon-field.html) have a different ordering: the energy bins from the larger to the lower.
 No redshift dependence is available.
 The surface is defined to include the nodes of the grid contained within.
 */
class TabularSpatialPhotonField: public PhotonField {
	public:
		TabularSpatialPhotonField(const std::string fieldName, ref_ptr<Surface> surface = nullptr);
		
		/** Returns photon density dependend on energy, redshift and position
		 * Returns the photon density for a specific photon energy, redshift and position
		 * @param ePhoton  Photon energy
		 * @param z  Redshift
		 * @param pos  Position
		 */
		double getPhotonDensity(double ePhoton = 0., double z = 0., const Vector3d &pos = Vector3d(0.,0.,0.)) const;
		/** Returns the minimum possible photon energy
		 * Returns the minimum possible photon energy for a given redshift and position.
		 * @param z  Redshift
		 * @param pos Position
		 */
		double getMinimumPhotonEnergy(double z, const Vector3d &pos = Vector3d(0.,0.,0.)) const;
		/** Returns the maximum possible photon energy
		 * Returns the maximum possible photon energy for a given redshift and position.
		 * @param z  Redshift
		 * @param pos Position
		 */
		double getMaximumPhotonEnergy(double z, const Vector3d &pos = Vector3d(0.,0.,0.)) const;

	protected:
		/** Reads the photon energies from the given folder
		 * @param filePath  Path containing photon energy files
		 */
		std::vector<double> readPhotonEnergy(std::string filePath);
		/** Reads the photon densities from the given folder
		 * @param filePath  Path containing photon densitiy files
		 */
		std::vector<double> readPhotonDensity(std::string filePath);
		/// Checks if read data is valid
		void checkInputData() const;
		
		/** Apply a surface that confine the position dependent photon field
		 * @param surface closed surface to confine the nodes of grid to be uploaded */
		void setSurface(ref_ptr<Surface> surface);
		
		// assuming all the nodes in the grid have the same energy binning
		std::vector<double> photonEnergies;
		std::vector<std::vector<double>> photonDensity;
		std::unordered_map<int, Vector3d> photonDict;
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
 @class IRB_Saldana21
 @brief Extragalactic background light model from Saldana-Lopez et al. 2021

 Source info:
 DOI:10.1093/mnras/stab2393
 https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.5144S/abstract
 */
class IRB_Saldana21: public TabularPhotonField {
public:
	IRB_Saldana21() : TabularPhotonField("IRB_Saldana21", true) {}
};

/**
 @class IRB_Saldana21_upper
 @brief Extragalactic background light model from Saldana-Lopez et al. 2021 (upper-bound model)

 Source info:
 DOI:10.1093/mnras/stab2393
 https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.5144S/abstract
 */
class IRB_Saldana21_upper: public TabularPhotonField {
public:
	IRB_Saldana21_upper() : TabularPhotonField("IRB_Saldana21_upper", true) {}
};

/**
 @class IRB_Saldana21_lower
 @brief Extragalactic background light model from Saldana-Lopez et al. 2021 (lower-bound model)

 Source info:
 DOI:10.1093/mnras/stab2393
 https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.5144S/abstract
 */
class IRB_Saldana21_lower: public TabularPhotonField {
public:
	IRB_Saldana21_lower() : TabularPhotonField("IRB_Saldana21_lower", true) {}
};

/**
 @class IRB_Finke22
 @brief Extragalactic background light model from Finke et al. 2022

 Source info:
 DOI:10.3847/1538-4357/ac9843
 https://iopscience.iop.org/article/10.3847/1538-4357/ac9843/pdf
 */
class IRB_Finke22: public TabularPhotonField {
public:
	IRB_Finke22() : TabularPhotonField("IRB_Finke22", true) {}
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
 @class ISRF
 @brief Interstellar radiation field model by Freudenreich et al. (1998) implemented in Porter et al. (2017)
 
 Source info:
 DOI:
 https://iopscience.iop.org/article/10.3847/1538-4357/aa844d
 */
class ISRF_Freudenreich98: public TabularSpatialPhotonField {
public:
	ISRF_Freudenreich98(ref_ptr<Surface> surface) : TabularSpatialPhotonField("ISRF_Freudenreich98", surface) {}
};

/**
 @class ISRF
 @brief Interstellar radiation field model by Robitaille et al. (2012) implemented in Porter et al. (2017)
 
 Source info:
 DOI:
 https://iopscience.iop.org/article/10.3847/1538-4357/aa844d
 */
class ISRF_Robitaille12: public TabularSpatialPhotonField {
public:
	ISRF_Robitaille12(ref_ptr<Surface> surface) : TabularSpatialPhotonField("ISRF_Robitaille12", surface) {}
};

/**
 @class BlackbodyPhotonField
 @brief Photon field decorator for black body photon fields.
 */
class BlackbodyPhotonField: public PhotonField {
public:
	/** Constructor
	 * @param fieldName  String identifier of the desired field
	 * @param blackbodyTemperatur  Blackbody temperature
	 */
	BlackbodyPhotonField(const std::string fieldName, const double blackbodyTemperature);

	/** Returns photon density dependend on energy, redshift and position
	 * Returns the photon density for a specific photon energy, redshift and position
	 * @param ePhoton  Photon energy
	 * @param z  Redshift
	 * @param pos  Position
	 */
	double getPhotonDensity(double ePhoton, double z = 0., const Vector3d &pos = Vector3d(0.,0.,0.)) const;
	/** Returns the minimum possible photon energy
	 * Returns the minimum possible photon energy for a given redshift and position.
	 * @param z  Redshift
	 * @param pos Position
	 */
	double getMinimumPhotonEnergy(double z, const Vector3d &pos = Vector3d(0.,0.,0.)) const;
	/** Returns the maximum possible photon energy
	 * Returns the maximum possible photon energy for a given redshift and position.
	 * @param z  Redshift
	 * @param pos Position
	 */
	double getMaximumPhotonEnergy(double z, const Vector3d &pos = Vector3d(0.,0.,0.)) const;
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


} // namespace crpropa

#endif // CRPROPA_PHOTONBACKGROUND_H
