#include "crpropa/PhotonPropagation.h"
#include "crpropa/Common.h"
#include "crpropa/Units.h"
#include "crpropa/Cosmology.h"
#include "crpropa/ProgressBar.h"

#include "EleCa/Propagation.h"
#include "EleCa/Particle.h"
#include "EleCa/Common.h"

#include "dint/prop_second.h"
#include "dint/DintEMCascade.h"

#include <fstream>
#include <cstdio>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <cmath>

namespace crpropa {

void ElecaPropagation(
		const std::string &inputfile,
		const std::string &outputfile,
		bool showProgress,
		double lowerEnergyThreshold,
		double magneticFieldStrength,
		const std::string &background) {

	std::ifstream infile(inputfile.c_str());
	std::streampos startPosition = infile.tellg();

	infile.seekg(0, std::ios::end);
	std::streampos endPosition = infile.tellg();
	infile.seekg(startPosition);


	ProgressBar progressbar(endPosition);
	if (showProgress) {
		progressbar.start("Run ElecaPropagation");
	}

	if (!infile.good())
		throw std::runtime_error(
				"ElecaPropagation: could not open file " + inputfile);

	eleca::setSeed();
	eleca::Propagation propagation;
	propagation.SetEthr(lowerEnergyThreshold / eV );
	propagation.ReadTables(getDataPath("EleCa/eleca.dat"));
	propagation.InitBkgArray(background);

	propagation.SetB(magneticFieldStrength / gauss);

	std::ofstream output(outputfile.c_str());
	output << "# ID\tE\tiID\tiE\tgeneration\n";
	output << "# ID          Id of particle (photon, electron, positron)\n";
	output << "# E           Energy [EeV]\n";
	output << "# iID         Id of source particle\n";
	output << "# iE          Energy [EeV] of source particle\n";
	output << "# Generation  number of interactions during propagation before particle is created\n";
	while (infile.good()) {
		if (infile.peek() != '#') {
			double E, D, pE, iE;
			int Id, pId, iId;
			infile >> Id >> E >> D >> pId >> pE >> iId >> iE;
			if (showProgress) {
				progressbar.setPosition(infile.tellg());
			}

			if (infile) {

				double z = eleca::Mpc2z(D);
				eleca::Particle p0(Id, E * 1e18, z);

				std::vector<eleca::Particle> ParticleAtMatrix;
				std::vector<eleca::Particle> ParticleAtGround;
				ParticleAtMatrix.push_back(p0);

				while (ParticleAtMatrix.size() > 0) {

					eleca::Particle p1 = ParticleAtMatrix.back();
					ParticleAtMatrix.pop_back();

					if (p1.IsGood()) {
						propagation.Propagate(p1, ParticleAtMatrix,
								ParticleAtGround);
					}
				}

				for (int i = 0; i < ParticleAtGround.size(); ++i) {
					eleca::Particle &p = ParticleAtGround[i];
					if (p.GetType() != 22)
						continue;
					char buffer[256];
					size_t bufferPos = 0;
					bufferPos += std::sprintf(buffer + bufferPos, "%i\t", p.GetType());
					bufferPos += std::sprintf(buffer + bufferPos, "%.4E\t", p.GetEnergy() / 1E18 );
					bufferPos += std::sprintf(buffer + bufferPos, "%i\t", iId);
					bufferPos += std::sprintf(buffer + bufferPos, "%.4E\t", iE );
					bufferPos += std::sprintf(buffer + bufferPos, "%i", p.Generation());
					bufferPos += std::sprintf(buffer + bufferPos, "\n");

					output.write(buffer, bufferPos);
				}
			}
		}

		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	infile.close();
	output.close();
}

typedef struct _Secondary {
	double E, D;
	int Id;
} _Secondary;

bool _SecondarySortPredicate(const _Secondary& s1, const _Secondary& s2) {
	return s1.D < s2.D;
}

void FillInSpectrum(Spectrum *a, const _Secondary &s) {
	double logE = log10(s.E) + 18;  // log10(E/eV)
	int iBin = floor((logE - MIN_ENERGY_EXP) / 0.1);  // bin number from 0 - NUM_MAIN_BINS-1
	if (iBin >= NUM_MAIN_BINS) {
		std::cout << "DintPropagation: Energy too high " << logE << std::endl;
		return;
	}
	if (iBin < 0) {
		std::cout << "DintPropagation: Energy too low " << logE << std::endl;
		return;
	}
	if (s.Id == 22)
		a->spectrum[PHOTON][iBin] += 1.;
	else if (s.Id == 11)
		a->spectrum[ELECTRON][iBin] += 1.;
	else if (s.Id == -11)
		a->spectrum[POSITRON][iBin] += 1.;
	else
		std::cout << "DintPropagation: Unhandled particle ID " << s.Id << std::endl;
}

void DintPropagation(
		const std::string &inputfile,
		const std::string &outputfile,
		int IRBFlag,
		int RadioFlag,
		double magneticFieldStrength,
		double aCutcascade_Magfield) {

	// initialize the energy grids for DINT
	dCVector energyGrid, energyWidth;
	New_dCVector(&energyGrid, NUM_MAIN_BINS);
	New_dCVector(&energyWidth, NUM_MAIN_BINS);
	SetEnergyBins(MIN_ENERGY_EXP, &energyGrid, &energyWidth);

	std::ofstream outfile(outputfile.c_str());
	if (!outfile.good())
		throw std::runtime_error(
				"DintPropagation: could not open file " + outputfile);

	std::ifstream infile(inputfile.c_str());
	if (!infile.good())
		throw std::runtime_error(
				"DintPropagation: could not open file " + inputfile);

	// initialize the spectrum
	Spectrum finalSpectrum;
	NewSpectrum(&finalSpectrum, NUM_MAIN_BINS);
	InitializeSpectrum(&finalSpectrum);

	std::string dataPath = getDataPath("dint");
	double B = magneticFieldStrength / gauss;
	double h = H0() * Mpc / 1000;
	DintEMCascade dint(IRBFlag, RadioFlag, dataPath, B, h, omegaM(), omegaL());

	const size_t nBuffer = 7.5E7;  // maximum number of simultaneously processed particles, keep memory requirement < 1GB
	const double dMargin = 0.1;  // distance bin width in [Mpc]

	while (infile.good()) {
		// read up to nBuffer secondaries from input file
		std::vector<_Secondary> secondaries;
		secondaries.reserve(nBuffer);
		while (infile.good() && (secondaries.size() < nBuffer)) {
			if (infile.peek() != '#') {
				_Secondary s;
				infile >> s.Id >> s.E >> s.D;
				s.D = comoving2LightTravelDistance(s.D * Mpc) / Mpc;  // DintEMCascade expects light travel distance
				if (infile)
					secondaries.push_back(s);
			}
			infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}

		if (secondaries.empty())
			break;  // all secondaries processed, nothing left to do

		// sort by distance
		std::sort(secondaries.begin(), secondaries.end(),
				_SecondarySortPredicate);

		Spectrum inputSpectrum, outputSpectrum;
		NewSpectrum(&inputSpectrum, NUM_MAIN_BINS);
		NewSpectrum(&outputSpectrum, NUM_MAIN_BINS);
		InitializeSpectrum(&inputSpectrum);

		// process secondaries
		while ((secondaries.back().D > 0) and (secondaries.size() > 0)) {
			double Dmax = secondaries.back().D;  // upper bound of distance bin
			double Dmin = max(Dmax - dMargin, 0.);  // lower bound of distance bin

			// add all secondaries within the current distance bin
			while ((secondaries.back().D > Dmin) and (secondaries.size() > 0)) {
				FillInSpectrum(&inputSpectrum, secondaries.back());
				secondaries.pop_back();
			}

			// propagate to next closest particle or to D=0
			double D = 0;
			if (secondaries.size() > 0)
				D = secondaries.back().D;

			// propagate distance step and make the output the new input spectrum
			InitializeSpectrum(&outputSpectrum);
			dint.propagate(Dmax, D, &inputSpectrum, &outputSpectrum, aCutcascade_Magfield);
			SetSpectrum(&inputSpectrum, &outputSpectrum);
		}

		// add remaining secondaries at D=0 to output spectrum
		while (secondaries.size() > 0) {
			FillInSpectrum(&inputSpectrum, secondaries.back());
			secondaries.pop_back();
		}

		AddSpectrum(&finalSpectrum, &inputSpectrum);
		DeleteSpectrum(&outputSpectrum);
		DeleteSpectrum(&inputSpectrum);
	}

	// output
	outfile << "# logE photons electrons positrons\n";
	outfile << "#   - logE: energy bin center <log10(E/eV)>\n";
	outfile << "#   - photons, electrons, positrons: total flux weights\n";
	for (int j = 0; j < finalSpectrum.numberOfMainBins; j++) {
		double logEc = MIN_ENERGY_EXP + 0.05 + j * 1. / BINS_PER_DECADE;
		outfile << std::setw(5) << logEc;
		for (int i = 0; i < 3; i++) {
			outfile << std::setw(13) << finalSpectrum.spectrum[i][j];
		}
		outfile << "\n";
	}
	outfile.close();

	DeleteSpectrum(&finalSpectrum);
	Delete_dCVector(&energyGrid);
	Delete_dCVector(&energyWidth);
}



bool _ParticlesAtGroundSortPredicate(const eleca::Particle& p1, const eleca::Particle& p2) {
	return p1.Getz() < p2.Getz();
}

void DintElecaPropagation(
		const std::string &inputfile,
		const std::string &outputfile,
		bool showProgress,
		double crossOverEnergy,
		double magneticFieldStrength,
		double aCutcascade_Magfield) {

	////////////////////////////////////////////////////////////////////////
	//Initialize EleCa
	std::ifstream infile(inputfile.c_str());
	std::streampos startPosition = infile.tellg();

	infile.seekg(0, std::ios::end);
	std::streampos endPosition = infile.tellg();
	infile.seekg(startPosition);

	ProgressBar progressbar(endPosition);
	if (showProgress) {
		progressbar.start("Run EleCa propagation");
	}

	if (!infile.good())
		throw std::runtime_error(
				"EleCaPropagation: could not open file " + inputfile);

	eleca::setSeed();
	eleca::Propagation propagation;
	propagation.SetEthr(crossOverEnergy / eV );
	propagation.ReadTables(getDataPath("EleCa/eleca.dat"));
	propagation.InitBkgArray("ALL");
	propagation.SetB(magneticFieldStrength / gauss);
	std::vector<eleca::Particle> ParticleAtGround;

	////////////////////////////////////////////////////////////////////////
	//Initialize DINT
	dCVector energyGrid, energyWidth;
	// Initialize the energy grids for dint
	New_dCVector(&energyGrid, NUM_MAIN_BINS);
	New_dCVector(&energyWidth, NUM_MAIN_BINS);
	SetEnergyBins(MIN_ENERGY_EXP, &energyGrid, &energyWidth);

	std::ofstream outfile(outputfile.c_str());
	if (!outfile.good())
		throw std::runtime_error(
				"DintPropagation: could not open file " + outputfile);

	Spectrum finalSpectrum;
	NewSpectrum(&finalSpectrum, NUM_MAIN_BINS);
	InitializeSpectrum(&finalSpectrum);

	std::string dataPath = getDataPath("dint");
	double h = H0() * Mpc / 1000;
	double ol = omegaL();
	double om = omegaM();
	DintEMCascade dint(4, 4, dataPath, magneticFieldStrength/gauss, h, om, ol);

	////////////////////////////////////////////////////////////////////////
	// Loop over infile

	while (infile.good()) {
		/// Eleca Propagation
		if (infile.peek() == '#') {
			infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			continue;
		}

		double E, D, pE, iE;
		int Id, pId, iId;
		infile >> Id >> E >> D >> pId >> pE >> iId >> iE;

		if (showProgress) {
			progressbar.setPosition(infile.tellg());
		}
		if (infile) { // stop at last line
			double z = eleca::Mpc2z(D);
			eleca::Particle p0(Id, E * 1e18, z);

			std::vector<eleca::Particle> ParticleAtMatrix;
			ParticleAtMatrix.push_back(p0);

			while (ParticleAtMatrix.size() > 0) {
				eleca::Particle p1 = ParticleAtMatrix.back();
				ParticleAtMatrix.pop_back();

				if (p1.IsGood()) {
					propagation.Propagate(p1, ParticleAtMatrix,
							ParticleAtGround, false);
				}
			}
		}

		// The vector is larger than ~1GB, or the infile is completely read - better call DINT.
		if (ParticleAtGround.size() > 1000000 || !infile) {
			const double dMargin = 0.1 * Mpc;

			std::sort(ParticleAtGround.begin(), ParticleAtGround.end(), _ParticlesAtGroundSortPredicate);

			Spectrum inputSpectrum, outputSpectrum;
			NewSpectrum(&inputSpectrum, NUM_MAIN_BINS);
			NewSpectrum(&outputSpectrum, NUM_MAIN_BINS);

			InitializeSpectrum(&inputSpectrum);
			// process secondaries
			while (ParticleAtGround.size() > 0) {
				double currentDistance =  redshift2LightTravelDistance(ParticleAtGround.back().Getz());  // dint expects light travel distance
				bool lastStep = (currentDistance == 0.);
				// add secondaries at the current distance to spectrum
				while ((ParticleAtGround.size() > 0) && (redshift2LightTravelDistance(ParticleAtGround.back().Getz()) >= (currentDistance - dMargin)))	{
					if (redshift2LightTravelDistance(ParticleAtGround.back().Getz()) > 0. || lastStep) {
						double criticalEnergy = ParticleAtGround.back().GetEnergy() / (ELECTRON_MASS); // units of dint
						int maxBin = (int) ((log10(criticalEnergy * ELECTRON_MASS) - MIN_ENERGY_EXP) * BINS_PER_DECADE + 0.5 + 1); // +1 line before to avoid conversion error to int for negative values (int(-0.7) = 0)
						maxBin -= 1; // remove the additional 1 from line before
						if (maxBin >= NUM_MAIN_BINS) {
							std::cout << "DintPropagation: Energy too high " <<
								ParticleAtGround.back().GetEnergy() << " eV"  <<
								std::endl;
							ParticleAtGround.pop_back();
							continue;
						}
						if (maxBin < 0) {
							std::cout << "DintPropagation: Energy too low " <<
								ParticleAtGround.back().GetEnergy() << " eV"  << std::endl;
							ParticleAtGround.pop_back();
							continue;
						}
						int Id = ParticleAtGround.back().GetType();
						if (Id == 22)
							inputSpectrum.spectrum[PHOTON][maxBin] += 1.;
						else if (Id == 11)
							inputSpectrum.spectrum[ELECTRON][maxBin] += 1.;
						else if (Id == -11)
							inputSpectrum.spectrum[POSITRON][maxBin] += 1.;
						else {
							std::cout << "DintPropagation: Unhandled particle ID " << Id
								<< std::endl;
						}
						ParticleAtGround.pop_back();
					} else
						break;
				}

				double D = 0;
				// only propagate to next particle
				if (ParticleAtGround.size() > 0)
					D = redshift2LightTravelDistance(ParticleAtGround.back().Getz());

				InitializeSpectrum(&outputSpectrum);
				dint.propagate(currentDistance / Mpc, D / Mpc, &inputSpectrum,
						&outputSpectrum, aCutcascade_Magfield);
				SetSpectrum(&inputSpectrum, &outputSpectrum);
			} // while (secondaries.size() > 0)

			AddSpectrum(&finalSpectrum, &inputSpectrum);

			DeleteSpectrum(&outputSpectrum);
			DeleteSpectrum(&inputSpectrum);
		} // dint call
	}

	infile.close();

	// output
	outfile << "# logE photons electrons positrons\n";
	outfile << "#   - logE: energy bin center <log10(E/eV)>\n";
	outfile << "#   - photons, electrons, positrons: total flux weights\n";
	for (int j = 0; j < finalSpectrum.numberOfMainBins; j++) {
		double logEc = MIN_ENERGY_EXP + 0.05 + j * 1. / BINS_PER_DECADE;
		outfile << std::setw(5) << logEc;
		for (int i = 0; i < 3; i++) {
			outfile << std::setw(13) << finalSpectrum.spectrum[i][j];
		}
		outfile << "\n";
	}
	outfile.close();

	DeleteSpectrum(&finalSpectrum);
	Delete_dCVector(&energyGrid);
	Delete_dCVector(&energyWidth);
}

} // namespace crpropa
