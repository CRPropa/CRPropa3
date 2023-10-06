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
@brief Candidates are split into n copies when they gain energy and cross specified energy bins. Weights are set accordingly.
		In case of Diffusice Shock Acceleration, splitting number can be adapted to expected spectral index to 
		compensate for the loss of particles per magnitude in energy
*/
class CandidateSplitting: public Module {
private:
	double nSplit;
	double minWeight;
	std::vector<double> Ebins;

public:

	CandidateSplitting();
	
	/** Constructor
	 @param nSplit 	Number of copies candidates are split 
	 @param Emin 		Minimal energy for splitting
	 @param Emax		Maximal energy for splitting
	 @param minWeight   Mimimal Weight
	 @param nBins		Number of energy bins 
	 */
	CandidateSplitting(int nSplit, double Emin, double Emax, double nBins, double minWeight);
	
	/** Constructor
	 @param nSplit 	Number of copies candidates are split 
	 @param Emin 		Minimal energy for splitting
	 @param Emax		Maximal energy for splitting
	 @param nBins		Number of energy bins 
	 @param minWeight   Mimimal Weight
	 @param log 		Energy bins in log
	 */
	CandidateSplitting(int nSplit, double Emin, double Emax, double nBins, double minWeight, bool log);
	
	/** Constructor
	 @param spectralIndex    Absolute value of expected spectral index determines splitting number 
	 @param Emin 			 Minimal energy for splitting
	 @param nBins            Number of bins in energy, with dE(spectralIndex) it determines Emax 
	 */
	CandidateSplitting(double spectralIndex, double Emin, int nBins);

	void process(Candidate *c) const;

	void setEnergyBins(double Emin, double Emax, double nBins, bool log);

	void setEnergyBinsDSA(double Emin, double dE, int n);

	void setNsplit(int n);

	void setMinimalWeight(double w);

	int getNsplit() const;

	double getMinimalWeight() const;

	const std::vector<double>& getEnergyBins() const;

};
/** @}*/

} // namespace crpropa
#endif // CRPROPA_CANDIDATESPLITTING_H
