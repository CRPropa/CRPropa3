#ifndef CRPROPA_EMCASCADE_H
#define CRPROPA_EMCASCADE_H

#include "crpropa/Module.h"

namespace crpropa {

/**
 @class EMCascade
 @brief Collects and deactivates photons, electrons and positrons. Uses DINT to calculate the EM cascade.
 */
class EMCascade: public Module {
private:
	// energy and distance binning
	int nE, nD;
	double logEmin, logEmax, dlogE, Dmax, dD;

	// histograms (distance,energy) of photons, electrons and positrons
	mutable std::vector<uint64_t> photonHist;
	mutable std::vector<uint64_t> electronHist;
	mutable std::vector<uint64_t> positronHist;
	void init();

public:
	EMCascade();

	/** Change the distance binning */
	void setDistanceBinning(
		double Dmax,  //!< maximum distance [m]
		int nD        //!< number of distance bins
		);

	/** Collect and deactivate photons, electrons and positrons */
	void process(Candidate *candidate) const;

	/** Save the unpropagated histogram of EM particles */
	void save(const std::string &filename);

	/** Load the unpropagated histogram of EM particles */
	void load(const std::string &filename);

	/** Calculates the EM cascade with DINT */
	void runCascade(
		const std::string &filename,  //!< output filename
		int IRBFlag = 4,        //!< EBL background 0: high, 1: low, 2: Primack, 4: Stecker'06
		int RadioFlag = 4,      //!< radio background 0: high, 1: medium, 2: obs, 3: none, 4: Protheroe'96
		double Bfield = 1E-13,  //!< magnetic field strength [T], default = 1 nG
		double cutCascade = 0   //!< a-parameter, see CRPropa 2 paper
		);

	std::string getDescription() const;
};

} // namespace crpropa

#endif // CRPROPA_EMCASCADE_H
