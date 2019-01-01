#ifndef CRPROPA_PHOTOPIONPRODUCTION_H
#define CRPROPA_PHOTOPIONPRODUCTION_H

// #include <unordered_map>  // ! c++11 and above only
#include <vector>
#include <string>

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"


namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 @class PhotoPionProduction
 @brief Photo-pion interactions of nuclei with background photons.
 */
class PhotoPionProduction: public Module {
protected:
	PhotonField photonField;
	ScalarGrid4d geometryGrid;
	Photon_Field phtnfld;
	std::vector<std::string> hashMap;  // contains histogram hashtags (workaround until c++11)
	std::vector< std::vector<double> > histData;  // contains histogram data (workaround until c++11)
	// std::unordered_map<std::string, std::vector<double> > particleMap;  // if c++11
	std::vector<double> tabLorentz; ///< Lorentz factor of nucleus
	std::vector<double> tabRedshifts;  ///< redshifts
	std::vector<double> tabProtonRate; ///< interaction rate in [1/m] for protons
	std::vector<double> tabNeutronRate; ///< interaction rate in [1/m] for neutrons
	double limit; ///< fraction of mean free path to limit the next step
	bool havePhotons;
	bool haveNeutrinos;
	bool haveElectrons;
	bool haveAntiNucleons;
	bool useTabulatedData;

public:
	PhotoPionProduction(
		PhotonField photonField = CMB,
		ScalarGrid4d geometryGrid = ScalarGrid4d(Vector3d(0.),0., 1,1,1,1, Vector3d(1.),1.),
		bool photons = false,
		bool neutrinos = false,
		bool electrons = false,
		bool antiNucleons = false,
		bool useTabulatedData = false,
		double limit = 0.1);
	void setPhotonField(PhotonField photonField);
	void setHavePhotons(bool b);
	void setHaveNeutrinos(bool b);
	void setHaveElectrons(bool b);
	void setHaveAntiNucleons(bool b);
	void setUseTabulatedData(bool b);
	void setLimit(double limit);
	void initRate(std::string filename);
	double nucleonMFP(double gamma, double z, bool onProton, Vector3d pos, double time) const;
	double nucleiModification(int A, int X) const;
	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate, bool onProton) const;
	// related to histogram version of SOPHIA
	void initHistogram(std::string filename);
	std::string hashTag(int nature, double Ein, double eps, int ID, int multiplicity) const;
	int produce(const std::vector<double> &particle) const;
	double drawEnergy(const std::vector<double> &data) const;
	double snapToHalfLog(double x) const;
	std::vector<double> sophiaEvent(bool onProton, double E, double e) const;
	/**
	 Calculates the loss length E dx/dE in [m].
	 This is not used in the simulation.
	 @param	id		PDG particle id
	 @param gamma	Lorentz factor of particle
	 @param z		redshift
	 */
	// double lossLength(int id, double gamma, double z = 0);
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_PHOTOPIONPRODUCTION_H
