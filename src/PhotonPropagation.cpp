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
#include <stdio.h>
#include <stdexcept>
#include <limits>
#include <iostream>
#include <algorithm>

namespace crpropa {




void EleCaPropagation(const std::string &inputfile,
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
		progressbar.start("Run EleCa propagation");
	}

	if (!infile.good())
		throw std::runtime_error(
				"EleCaPropagation: could not open file " + inputfile);

	eleca::setSeed();
	eleca::Propagation propagation;
  propagation.SetEthr(lowerEnergyThreshold / eV );
	propagation.ReadTables(getDataPath("EleCa/eleca.dat"));
	propagation.InitBkgArray(background);

	propagation.SetB(magneticFieldStrength / gauss);

	std::ofstream output(outputfile.c_str());
	output << "# ID\tE\tiID\tiE\n";
	output << "# ID          Id of particle (photon, electron, positron)\n";
	output << "# E           Energy [EeV]\n";
	output << "# iID         Id of source particle\n";
	output << "# iE          Energy [EeV] of source particle\n";
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
					bufferPos += sprintf(buffer + bufferPos, "%i\t", p.GetType());
					bufferPos += sprintf(buffer + bufferPos, "%.4E\t", p.GetEnergy() / 1E18 );
					bufferPos += sprintf(buffer + bufferPos, "%i\t", iId);
					bufferPos += sprintf(buffer + bufferPos, "%.4E", iE );
					bufferPos += sprintf(buffer + bufferPos, "\n");

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
void AddSpectrum(Spectrum *a, const Spectrum *b) {
	for (int i = 0; i < NUM_SPECIES; i++) {
		for (int j = 0; j < a->numberOfMainBins; j++)
			a->spectrum[i][j] += b->spectrum[i][j];
	}
}


void DintPropagation(const std::string &inputfile,
		const std::string &outputfile, double magneticFieldStrength,  int IRFlag, int RadioFlag, double Zmax) {
	// Initialize the spectrum
	dCVector energyGrid, energyWidth;
	// Initialize the energy grids for dint
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

	// Initialize the bField
	dCVector bField;
	New_dCVector(&bField, 5);
	for (size_t i = 0; i < 5; i++)	bField.vector[i] = magneticFieldStrength / gauss;  

	Spectrum finalSpectrum;
	NewSpectrum(&finalSpectrum, NUM_MAIN_BINS);
	InitializeSpectrum(&finalSpectrum);

	std::string dataPath = getDataPath("dint");
	double h = H0() * Mpc / 1000;
	double ol = omegaL();
	double om = omegaM();
	DintEMCascade dint(IRFlag, RadioFlag, dataPath, magneticFieldStrength/gauss, h, om, ol);

	const size_t nBuffer = 1 << 20;
	const double dMargin = 0.1; // Mpc;

	size_t cnt = 0;
	while (infile.good()) {
		// buffer for secondaries
		std::vector<_Secondary> secondaries;
		secondaries.reserve(nBuffer);

		// read secondaries from file
		size_t n = 0;
		while (infile.good() && (n < nBuffer)) {
			if (infile.peek() != '#') {
				double pE, iE;
				int pId, iId;
				_Secondary s;
				infile >> s.Id >> s.E >> s.D >> pId >> pE >> iId >> iE;
				if (infile) {
					n++;
					secondaries.push_back(s);
				}
			}

			infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}

		if (secondaries.size() == 0)
			continue;

		// sort by D
		std::sort(secondaries.begin(), secondaries.end(),
				_SecondarySortPredicate);


		Spectrum inputSpectrum, outputSpectrum;
		NewSpectrum(&inputSpectrum, NUM_MAIN_BINS);
		NewSpectrum(&outputSpectrum, NUM_MAIN_BINS);
		InitializeSpectrum(&inputSpectrum);
		// process secondaries
		while (secondaries.size() > 0) {
			double currentDistance = secondaries.back().D;
			// add secondaries at the current distance to spectrum
			while ((secondaries.size() > 0) && (secondaries.back().D >= (currentDistance - dMargin))) {
				double criticalEnergy = secondaries.back().E * EeV / (eV * ELECTRON_MASS); // units of dint
				int maxBin = (int) ((log10(criticalEnergy * ELECTRON_MASS)
						- MAX_ENERGY_EXP) * BINS_PER_DECADE + NUM_MAIN_BINS);
				if (maxBin >= NUM_MAIN_BINS) {
					std::cout << "DintPropagation: Energy too high " << secondaries.back().E
							<< std::endl;
					secondaries.pop_back();
					continue;
				}
				if (maxBin < 0) {
					std::cout << "DintPropagation: Energy too low " << secondaries.back().E
							<< std::endl;
					secondaries.pop_back();
					continue;
				}
				int Id = secondaries.back().Id;
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
				secondaries.pop_back();
			}

			double D = 0;
			// only propagate to next particle
			if (secondaries.size() > 0)
				D = secondaries.back().D;

			InitializeSpectrum(&outputSpectrum);
			//prop_second(currentDistance, &bField, &energyGrid, &energyWidth, &inputSpectrum,
			//		&outputSpectrum, dataPath, IRFlag, Zmax, RadioFlag, h, om,
			//		ol, 0);

			dint.propagate(currentDistance, D, &inputSpectrum, &outputSpectrum);
			SetSpectrum(&inputSpectrum, &outputSpectrum);
		}

		AddSpectrum(&finalSpectrum, &inputSpectrum);

		DeleteSpectrum(&outputSpectrum);
		DeleteSpectrum(&inputSpectrum);
	}

	outfile << "# BinCenter [EeV] BinWidth [EeV] Flux-Weights for photons electrons positrons ... \n";
	for (int j = 0; j < finalSpectrum.numberOfMainBins; j++) {
		outfile << (energyGrid.vector[j] / EeV * (eV * ELECTRON_MASS)) << " ";
		outfile << (energyWidth.vector[j] / EeV * (eV * ELECTRON_MASS)) << " ";
		for (int i = 0; i < NUM_SPECIES; i++) {
			outfile << finalSpectrum.spectrum[i][j] << " ";
		}
		outfile << "\n";
	}


	DeleteSpectrum(&finalSpectrum);
	Delete_dCVector(&bField);

	Delete_dCVector(&energyGrid);
	Delete_dCVector(&energyWidth);

}



bool _ParticlesAtGroundSortPredicate(const eleca::Particle& p1, const eleca::Particle& p2) {
	return p1.Getz() < p2.Getz();
}


void DintElcaPropagation(const std::string &inputfile,
	const std::string &outputfile, 
	bool showProgress,
	double crossOverEnergy,
	double magneticFieldStrength)
{
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
		if (infile.peek() == '#')
		{
			infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			continue;
		}
	
		double E, D, pE, iE;
		int Id, pId, iId;
		infile >> Id >> E >> D >> pId >> pE >> iId >> iE;

		if (showProgress) {
			progressbar.setPosition(infile.tellg());
		}
		if (infile)
		{ // stop at last line
			double z = eleca::Mpc2z(D);
			eleca::Particle p0(Id, E * 1e18, z);

			std::vector<eleca::Particle> ParticleAtMatrix;
			ParticleAtMatrix.push_back(p0);

			while (ParticleAtMatrix.size() > 0) 
			{
				eleca::Particle p1 = ParticleAtMatrix.back();
				ParticleAtMatrix.pop_back();

				if (p1.IsGood()) {
					propagation.Propagate(p1, ParticleAtMatrix,
							ParticleAtGround, false);
				}
			}
		}

		if (ParticleAtGround.size() > 1000000 || !infile) // The vector is larger
			//than ~1GB, or the infile is completley read - better call DINT.
		{
				const double dMargin = 0.1 * Mpc; 
				size_t cnt = 0;
				
				std::sort(ParticleAtGround.begin(), ParticleAtGround.end(), _ParticlesAtGroundSortPredicate);

				Spectrum inputSpectrum, outputSpectrum;
				NewSpectrum(&inputSpectrum, NUM_MAIN_BINS);
				NewSpectrum(&outputSpectrum, NUM_MAIN_BINS);

				InitializeSpectrum(&inputSpectrum);
				// process secondaries
				while (ParticleAtGround.size() > 0) 
				{
					double currentDistance =  redshift2ComovingDistance(ParticleAtGround.back().Getz()) ;
					// add secondaries at the current distance to spectrum
					while ((ParticleAtGround.size() > 0) && (redshift2ComovingDistance(ParticleAtGround.back().Getz()) >= (currentDistance - dMargin))) 
					{
						double criticalEnergy = ParticleAtGround.back().GetEnergy() / (ELECTRON_MASS); // units of dint
						int maxBin = (int) ((log10(criticalEnergy * ELECTRON_MASS)
								- MAX_ENERGY_EXP) * BINS_PER_DECADE + NUM_MAIN_BINS);
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
					}

					double D = 0;
					// only propagate to next particle
					if (ParticleAtGround.size() > 0)
						D = redshift2ComovingDistance(ParticleAtGround.back().Getz());

					InitializeSpectrum(&outputSpectrum);
					dint.propagate(currentDistance / Mpc, D / Mpc, &inputSpectrum,
							&outputSpectrum);
					SetSpectrum(&inputSpectrum, &outputSpectrum);
				} // while (secondaries.size() > 0) 
	
				AddSpectrum(&finalSpectrum, &inputSpectrum);
	
				DeleteSpectrum(&outputSpectrum);
				DeleteSpectrum(&inputSpectrum);
		} // dint call
	}

	infile.close();
	// output spectrum
	outfile << "# BinCenter [EeV] BinWidth [EeV] Flux-Weights for photons electrons positrons ... \n";
	for (int j = 0; j < finalSpectrum.numberOfMainBins; j++) 
	{
		outfile << (energyGrid.vector[j] / EeV * (eV * ELECTRON_MASS)) << " ";
		outfile << (energyWidth.vector[j] / EeV * (eV * ELECTRON_MASS)) << " ";
		for (int i = 0; i < NUM_SPECIES; i++) {
			outfile << finalSpectrum.spectrum[i][j] << " ";
		}
		outfile << "\n";
	}

	DeleteSpectrum(&finalSpectrum);

	Delete_dCVector(&energyGrid);
	Delete_dCVector(&energyWidth);

}





} // namespace crpropa
