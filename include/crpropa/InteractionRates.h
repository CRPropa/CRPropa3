#ifndef CRPROPA_INTERACTIONRATES_H
#define CRPROPA_INTERACTIONRATES_H

#include "crpropa/Common.h"
#include "crpropa/Referenced.h"
#include "crpropa/Vector3.h"
#include "crpropa/Geometry.h"

#include <nanoflann.hpp>

#include <vector>
#include <string>
#include <unordered_map>

namespace crpropa {

struct PointCloud {
	std::vector<Vector3d> points;
	std::vector<int> ids;

	inline size_t kdtree_get_point_count() const { return points.size(); }

	inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
		if (dim == 0) return points[idx].x;
		if (dim == 1) return points[idx].y;
		return points[idx].z;
	}

	// optional bounding-box computation (required by nanoflann)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX& /*bb*/) const {
		return false;  // no bounding box optimization
	}
	
};

using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
	nanoflann::L2_Simple_Adaptor<double, PointCloud>,
	PointCloud,
	3
>;

/**
 * \addtogroup InteractionRates
 * @{
 */

/**
 @class Interaction Rates
 @brief Abstract base class for photon fields interaction rates.
 */
class InteractionRates: public Referenced {
public:
	virtual double getProcessRate(const double E, const Vector3d &position) const = 0;
	virtual void loadPerformInteractionTabs(const Vector3d &position, std::vector<double> &tabE, std::vector<double> &tabs, std::vector<std::vector<double>> &tabCDF) const = 0;
	
	std::string getRatesName() const {
		return this->ratesName;
	}
	
	bool hasPositionDependence() const {
		return this->isPositionDependent;
	}
	
	void setRatesName(std::string ratesName) {
		this->ratesName = ratesName;
	}

	virtual void initRate(std::string path) = 0;
	virtual void initCumulativeRate(std::string path) = 0;

protected: 

  std::string ratesName = "AbstractInteractionRates";
  bool isPositionDependent = false; 

};

/**
 @class InteractionRateHomogeneous
 @brief Interaction rates decorator for tabulated homogeneous photon fields.
 */
class InteractionRatesHomogeneous: public InteractionRates {
public:
	/** Constructor of InteractionRatesHomogeneous
	 * @param RateFile Path to the file containing the interaction rate data
	 * @param CumulativeRateFile Path to the file containing the cumulative interaction rate data
	 */
	InteractionRatesHomogeneous(std::string RateFile = "", std::string CumulativeRateFile = "");
	
	std::vector<double> getTabulatedEnergy() const;
	std::vector<double> getTabulatedRate() const;
	std::vector<double> getTabulatedE() const;
	std::vector<double> getTabulateds() const;
	std::vector<std::vector<double>> getTabulatedCDF() const;
	
	double getProcessRate(const double E, const Vector3d &position) const;
	void loadPerformInteractionTabs(const Vector3d &position, std::vector<double> &tabE, std::vector<double> &tabs, std::vector<std::vector<double>> &tabCDF) const;
	
	void setTabulatedEnergy (std::vector<double>& tabEnergy);
	void setTabulatedRate (std::vector<double>& tabRate);
	void setTabulatedE (std::vector<double>& tabE);
	void setTabulateds (std::vector<double>& tabs);
	void setTabulatedCDF (std::vector<std::vector<double>>& tabCDF);

	/** Loads the interaction rate
	 * This function loads the interaction rate
	 * @param filename The name of the file containing the interaction rates
	 */
	void initRate(std::string filename);
	/** Loads the cumulative interaction rate
	 * This function loads the interaction rate
	 * @param filename The name of the file containing the interaction rates
	 */
	void initCumulativeRate(std::string filename);
	
protected:
	
	// tabulated interaction rates 1/lambda(E)
	std::vector<double> tabEnergy; //!< electron energy in [J]
	std::vector<double> tabRate; //!< interaction rate in [1/m]
	
	// tabulated CDF(s_kin, E) = cumulative differential interaction rate
	std::vector<double> tabE; //!< electron energy in [J]
	std::vector<double> tabs; //!< s_kin = s - m^2 in [J**2]
	std::vector<std::vector<double>> tabCDF; //!< cumulative interaction rate
	
};

/**
 @class InteractionRatePositionDependent
 @brief Interaction rates decorator for tabulated position dependent photon fields.
 */
class InteractionRatesPositionDependent: public InteractionRates {

public:
	/** Constructor of InteractioNRatesPositionDependent
	 * @param RateFilePath Path containing the interaction rates files (* /Rate)
	 * @param CumulativeRateFilePath Path containing the cumulative interaction rates files (* /CumulativeRate)
	 * @param surface Closed surface to confine the grid nodes to be uploaded (optional)
	 */
	InteractionRatesPositionDependent(std::string RateFilePath = "", std::string CumulativeRateFilePath = "", ref_ptr<Surface> surface = NULL);
	
	int findClosestGridPoint(const Vector3d &position) const;
	
	std::vector<double> getTabulatedEnergy() const;
	std::vector<std::vector<double>> getTabulatedRate() const;
	std::vector<double> getTabulatedE() const;
	std::vector<std::vector<double>> getTabulateds() const;
	std::vector<std::vector<std::vector<double>>> getTabulatedCDF() const;
	std::unordered_map<int, Vector3d> getPhotonDict() const;
	std::vector<double> getClosestRate(const Vector3d &position) const;
	std::vector<double> getClosests(const Vector3d &position) const;
	std::vector<std::vector<double>> getClosestCDF(const Vector3d &position) const;
	
	double getProcessRate(const double E, const Vector3d &position) const;
	void loadPerformInteractionTabs(const Vector3d &position, std::vector<double> &tabE, std::vector<double> &tabs, std::vector<std::vector<double>> &tabCDF) const;
	
	void setTabulatedEnergy (std::vector<double>& tabEnergy);
	void setTabulatedRate (std::vector<std::vector<double>>& tabRate);
	void setTabulatedE (std::vector<double>& tabE);
	void setTabulateds (std::vector<std::vector<double>>& tabs);
	void setTabulatedCDF (std::vector<std::vector<std::vector<double>>>& tabCDF);
	void setPhotonDict (std::unordered_map<int, Vector3d>& photonDict);

	/** Apply a surface that confine the position dependent photon field
	 * @param surface Closed surface to confine the grid nodes to be uploaded
	 */
	void setSurface(ref_ptr<Surface> surface);
	ref_ptr<Surface> getSurface() const;

	/** Loads the interaction rate
	 * This function loads the position dependent interaction rate
	 * @param filepath The name of the folder containing the interaction rates (* /Rate)
	 */
	void initRate(std::string filepath);
	/** Loads the interaction rate
	 * This function loads the cumulative position dependent interaction rate
	 * @param filepath The name of the folder containing the cumulative interaction rates (* /CumulativeRate)
	 */
	void initCumulativeRate(std::string filepath);

protected:
	
	// tabulated interaction rates 1/lambda(E)
	std::vector<double> tabEnergy; //!< electron energy in [J], assuming the same energy binning in each node
	std::vector<std::vector<double>> tabRate; //!< interaction rate in [1/m]
	
	// tabulated CDF(s_kin, E) = cumulative differential interaction rate
	std::vector<double> tabE; //!< electron energy in [J], assuming the same energy binning in each node
	std::vector<std::vector<double>> tabs; //!< s_kin = s - m^2 in [J**2]
	std::vector<std::vector<std::vector<double>>> tabCDF; //!< cumulative interaction rate
	std::unordered_map<int, Vector3d> photonDict; //!< dictionary to link tables to spatial coordinates
	
	PointCloud cloud; //!< point cloud for nanoflann KD-tree
	KDTree* tree = nullptr; //!< pointer to the KD Tree
	ref_ptr<Surface> surface;
	
};

} // namespace crpropa

#endif // CRPROPA_INTERACTIONRATES_H
