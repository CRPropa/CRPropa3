#include "crpropa/module/CandidateSplitting.h"

namespace crpropa {

CandidateSplitting::CandidateSplitting() {
	// no particle splitting if EnergyBins and NSplit are not specified
	setNsplit(0);
	setMinimalWeight(1.);
}

CandidateSplitting::CandidateSplitting(int nSplit, double Emin, double Emax,  double nBins, double minWeight, bool log) {
	setNsplit(nSplit);
	setEnergyBins(Emin, Emax, nBins, log);
	setMinimalWeight(minWeight);
}

CandidateSplitting::CandidateSplitting(double spectralIndex, double Emin, int nBins)  {
	// to use with Diffusive Shock Acceleration
	if (spectralIndex > 0){
		throw std::runtime_error(
				"CandidateSplitting: spectralIndex > 0 !"); 
	}

	setNsplit(2); // always split in 2, calculate bins in energy for given spectrum:
	double dE = pow(1. / 2, 1. / (spectralIndex + 1)); 
	setEnergyBinsDSA(Emin, dE, nBins);
	setMinimalWeight(1. / pow(2, nBins));
}

void CandidateSplitting::process(Candidate *c) const {
	double currE = c->current.getEnergy(); 
	double prevE = c->previous.getEnergy();

	if (c->getWeight() <= minWeight){
		// minimal weight reached, no splitting
		return;
	}
	if (currE < Ebins[0] || nSplit == 0 ){
		// current energy is smaller than first bin -> no splitting
		// or, number of splits = 0
		return;
	}
	for (size_t i = 0; i < Ebins.size(); ++i){
		
		if( prevE < Ebins[i] ){
			// previous energy is in energy bin [i-1, i]
			if(currE < Ebins[i]){
				//assuming that dE greater than 0, prevE and E in same energy bin -> no splitting
				return;
			}
			// current energy is in energy bin [i,i+1] or higher -> particle splitting for each crossing
			for (size_t j = i; j < Ebins.size(); ++j ){

				// adapted from Acceleration Module:
				c->updateWeight(1. / nSplit); // * 1/n_split

				for (int i = 1; i < nSplit; i++) {
				
					ref_ptr<Candidate> new_candidate = c->clone(false);
					//InteractionTag is PRIM, physically no new particles are created
					new_candidate->parent = c;
					new_candidate->previous.setEnergy(currE); // so that new candidate is not split again in next step!
					c->addSecondary(new_candidate);
				}
				if (j < Ebins.size()-1 && currE < Ebins[j+1]){
					// candidate is in energy bin [j, j+1] -> no further splitting
					return;
				}
			}
			return;
		}
	}
}

void CandidateSplitting::setEnergyBins(double Emin, double Emax, double nBins, bool log) {
	Ebins.resize(0);
	if (Emin > Emax){
		throw std::runtime_error(
				"CandidateSplitting: Emin > Emax!");
	}
	double dE = (Emax-Emin)/nBins;
	for (size_t i = 0; i < nBins; ++i) {
		if (log == true) {
			Ebins.push_back(Emin * pow(Emax / Emin, i / (nBins - 1.0)));
		} else {
			Ebins.push_back(Emin + i * dE);
		}
	}
}

void CandidateSplitting::setEnergyBinsDSA(double Emin, double dE, int n) {
	Ebins.resize(0);
	for (size_t i = 1; i < n + 1; ++i) {
		Ebins.push_back(Emin * pow(dE, i));
	}
}

const std::vector<double>& CandidateSplitting::getEnergyBins() const {
	return Ebins;
}

void CandidateSplitting::setNsplit(int n) {
	nSplit = n;
}

void CandidateSplitting::setMinimalWeight(double w) {
	minWeight = w;
}

int CandidateSplitting::getNsplit() const {
	return nSplit;
}

double CandidateSplitting::getMinimalWeight() const {
	return minWeight;
}

} // end namespace crpropa

