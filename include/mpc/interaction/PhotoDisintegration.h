#ifndef PHOTODISINTEGRATION_H_
#define PHOTODISINTEGRATION_H_

#include "mpc/Module.h"
#include "mpc/ParticleState.h"

#include <fstream>

namespace mpc {

struct IndexStruct {
	int id;
	size_t start;
	size_t end;
};

class PhotoDisintegration: public Module {
private:
	int n;
	double dx;
	double xMin;
	std::vector<IndexStruct> nucleusIndex;
	std::vector<IndexStruct> channelIndex;
	std::vector<double> tabulatedMeanFreePath;

public:
	PhotoDisintegration() {
		IndexStruct index;
		char str[256];
		size_t dev0;

		// load nucleus index
		std::ifstream infile("data/PDInitNucleusId.cmt");
		infile.getline(str, 255); // skip header
		while (!infile.eof()) {
			infile >> index.id >> dev0 >> dev0 >> index.start >> index.end;
			nucleusIndex.push_back(index);
		}

		// load channel index
		 infile.open("data/PDExclTabMnFrPthCrossId.cmt");
		infile.getline(str, 255); // skip header
		while (!infile.eof()) {
			infile >> index.id >> index.start;
			channelIndex.push_back(index);
		}

		// load mean free path values
		double value;
		infile.open("data/PDExclTabMnFrPthCross.cmt");
		infile.getline(str, 255); // skip header
		while (!infile.eof()) {
			infile >> value;
			tabulatedMeanFreePath.push_back(value);
		}
	}

	~PhotoDisintegration();

	void process(Candidate *candidate, std::vector<Candidate *> &secondaries) {
	}

	double getLambdaSum(Candidate candidate) {
		int id = candidate.current.getId();
		double logGamma = log(candidate.current.getLorentzFactor());

		size_t start = nucleusIndex[id].start;
		size_t end = nucleusIndex[id].end;

		double mfp = 0;
		for (size_t i = start; i < end; i++)
			mfp += interpolate(channelIndex[i].start, logGamma);

		return 1. / mfp;
	}

	double interpolate(size_t start, double x) {
		// bin below
		size_t i = size_t((x - xMin) / dx) + start;
		//energy below / above
		double x0 = xMin + i * dx;
		double x1 = x0 + dx;
		//value below / above
		double y0 = tabulatedMeanFreePath[i];
		double y1 = tabulatedMeanFreePath[i];
		return (y0 + (y0 - y1) * (x - x0) / dx);
	}

	void disintegrate(Candidate *candidate,
			std::vector<Candidate *> &secondaries, size_t channel) {

		int id = candidate->current.getId();

		// disintegration
		int nNeutron = channel / 100000;
		int nProton = (channel % 100000) / 10000;
		int nDeuterium = (channel % 10000) / 1000;
		int nTritium = (channel % 1000) / 100;
		int nHelium3 = (channel % 100) / 10;
		int nHelium4 = (channel % 10);
		int dA = nNeutron + nProton + 2 * nDeuterium + 3 * nTritium
				+ 3 * nHelium3 + 4 * nHelium4;
		int dZ = nProton + nDeuterium + nTritium + 2 * nHelium3 + 2 * nHelium4;

		// particle
		int A = getMassNumberFromNucleusId(id);
		int Z = getChargeNumberFromNucleusId(id);
		double energyPerNucleon = candidate->current.getEnergy() / double(A);

		// update particle
		A -= dA;
		Z -= dA;
		candidate->current.setId(getNucleusId(A, Z));
		candidate->current.setEnergy(energyPerNucleon * A);

		// create secondaries
		for (int i = 0; i < nNeutron; i++) {
			createSecondary(candidate, secondaries, 2112, energyPerNucleon);
		}
		for (int i = 0; i < nProton; i++) {
			createSecondary(candidate, secondaries, 2212, energyPerNucleon);
		}
		for (int i = 0; i < nDeuterium; i++) {
			createSecondary(candidate, secondaries, 1000010020,
					energyPerNucleon * 2);
		}
		for (int i = 0; i < nTritium; i++) {
			createSecondary(candidate, secondaries, 1000010030,
					energyPerNucleon * 3);
		}
		for (int i = 0; i < nHelium3; i++) {
			createSecondary(candidate, secondaries, 1000020030,
					energyPerNucleon * 3);
		}
		for (int i = 0; i < nHelium4; i++) {
			createSecondary(candidate, secondaries, 1000020040,
					energyPerNucleon * 4);
		}
	}

	void createSecondary(Candidate *candidate,
			std::vector<Candidate *> &secondaries, int id, double energy) {
		ParticleState initial = candidate->current;
		initial.setEnergy(energy);
		initial.setId(id);
		Candidate secondary;
		secondary.current = initial;
		secondary.initial = initial;
		secondary.setNextStep(candidate->getCurrentStep());
		secondaries.push_back(&secondary);
	}

};

} // namespace mpc

#endif /* PHOTODISINTEGRATION_H_ */
