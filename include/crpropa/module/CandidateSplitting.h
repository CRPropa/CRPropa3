#ifndef CRPROPA_CANDIDATEPLITTING_H
#define CRPROPA_CANDIDATEPLITTING_H

#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>

#include "crpropa/Vector3.h"
#include "crpropa/Module.h"
#include "crpropa/Units.h"
#include "kiss/logger.h"


namespace crpropa {

/**
@class CandidateSplitting
@brief Candidates are split into n copies when they cross specified energy bins. Weights are set accordingly.
		In case of Diffusice Shock Acceleration, splitting number can be adapted to expected spectral index to 
		compensate for the loss of particles per magnitude in energy
*/

class CandidateSplitting: public Module {
private:
	double n_split;
	double minWeight;
	std::vector<double> Ebins;

public:

	CandidateSplitting();
	
	CandidateSplitting(int n_split, double Emin, double Emax, double n_bins, double minWeight);
	/** Constructor
	 @param n_split 	Number of copies candidates are split 
	 @param Emin 		Minimal energy for splitting
	 @param Emax		Maximal energy for splitting
	 @param minWeight   Mimimal Weight
	 @param n_bins		Number of energy bins 
	 */
	CandidateSplitting(int n_split, double Emin, double Emax, double n_bins, double minWeight, bool log);
	/** Constructor
	 @param n_split 	Number of copies candidates are split 
	 @param Emin 		Minimal energy for splitting
	 @param Emax		Maximal energy for splitting
	 @param n_bins		Number of energy bins 
	 @param minWeight   Mimimal Weight
	 @param log 		Energy bins in log
	 */

	CandidateSplitting(double SpectralIndex, double Emin, int nBins);
	/** Constructor
	 @param SpectralIndex    Absolute value of expected spectral index determines splitting number 
	 @param Emin 			 Minimal energy for splitting
	 @param nBins            Number of bins in energy, with dE(spectralIndex) it determines Emax 
	 */

	void process(Candidate *c) const;

	void setEnergyBins(double Emin, double Emax, double n_bins, bool log);

	void setEnergyBinsDSA(double Emin, double dE, int n);

	void setNsplit(int n);

	void setMinimalWeight(double w);

	int getNsplit() const;

	double getMinimalWeight() const;

	const std::vector<double>& getEnergyBins() const;

};
/** @}*/

} // namespace crpropa
#endif // CRPROPA_PARTICLESPLITTING_H
