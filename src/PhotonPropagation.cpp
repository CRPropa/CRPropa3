#include "crpropa/PhotonPropagation.h"
#include "crpropa/Common.h"
#include "crpropa/Units.h"
#include "crpropa/Cosmology.h"

#include "EleCa/Propagation.h"
#include "EleCa/Particle.h"
#include "EleCa/Common.h"

#include "dint/prop_second.h"

#include <fstream>
#include <stdexcept>
#include <limits>
#include <iostream>

namespace crpropa {

void EleCaPropagation(const std::string &inputfile,
		const std::string &background, std::vector<double> &energy,
		std::vector<double> &spectrum) {
	std::ifstream infile(inputfile.c_str());

	const double emin = 16, emax = 22, step = 0.2;
	const size_t steps = (emax - emin) / step;
	energy.clear();
	energy.resize(steps);
	spectrum.resize(steps);
	for (size_t i = 0; i < steps; i++) {
		energy[i] = 16 + i * step;
		spectrum[i] = 0;
	}

	if (!infile.good())
		throw std::runtime_error(
				"EleCaPropagation: could not open file " + inputfile);

	eleca::Propagation propagation;
	propagation.ReadTables(getDataPath("EleCa/eleca.dat"));
	propagation.InitBkgArray(background);

	while (infile.good()) {
		if (infile.peek() != '#') {
			double E, D, pE, iE;
			int Id, pId, iId;
			infile >> Id >> E >> D >> pId >> pE >> iId >> iE;
			if (infile) {
				double z = eleca::Mpc2z(D);
				eleca::Particle p0(Id, E * 1e18, z);

				// TODO: find a motivated value!
				p0.SetB(1e-9);

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

				//propagation.WriteOutput(output, p0, ParticleAtGround);
				for (int i = 0; i < ParticleAtGround.size(); ++i) {
					eleca::Particle &p = ParticleAtGround[i];
					if (p.GetType() != 22)
						continue;
					size_t idx = (::log10(p.GetEnergy()) - emin) / step;
					spectrum.at(idx) += 1;
				}
			}
		}

		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	infile.close();
}

void EleCaPropagation(const std::string &inputfile,
		const std::string &outputfile, const std::string &background) {
	std::vector<double> energy, spectrum;
	EleCaPropagation(inputfile, background, energy, spectrum);
	std::ofstream output(outputfile.c_str());
	output << "# E N\n";
	for (size_t i = 0; i < energy.size(); i++) {
		output << energy[i] << " " << spectrum[i] << "\n";
	}
}

void DintPropagation(const std::string &inputfile,
		const std::string &outputfile, int IRFlag, int RadioFlag, double Zmax) {
	// Initialize the spectrum
	Spectrum inputSpectrum;
	NewSpectrum(&inputSpectrum, NUM_MAIN_BINS);

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
	New_dCVector(&bField, 1);

	// Initialize output spectrum
	Spectrum outputSpectrum;
	NewSpectrum(&outputSpectrum, NUM_MAIN_BINS);

	Spectrum finalSpectrum;
	NewSpectrum(&finalSpectrum, NUM_MAIN_BINS);

	std::string dataPath = getDataPath("dint");

	size_t cnt = 0;
	while (infile.good()) {
		if (infile.peek() != '#') {
			double E, D, pE, iE;
			int Id, pId, iId;
			infile >> Id >> E >> D >> pId >> pE >> iId >> iE;
			cnt++;
			std::cerr << cnt << std::endl;
			if (infile) {

				double criticalEnergy = E * EeV / (eV * ELECTRON_MASS); // units of dint
				int maxBin = (int) ((log10(criticalEnergy * ELECTRON_MASS)
						- MAX_ENERGY_EXP) * BINS_PER_DECADE + NUM_MAIN_BINS);
				if (Id == 22)
					inputSpectrum.spectrum[PHOTON][maxBin] = 1.;
				else if (Id == 11)
					inputSpectrum.spectrum[ELECTRON][maxBin] = 1.;
				else if (Id == -11)
					inputSpectrum.spectrum[POSITRON][maxBin] = 1.;
				else {
					std::cerr << "DintPropagation: Unhandled particle ID " << Id
							<< std::endl;
					continue;
				}

				double h = H0() * Mpc / 1000;
				double ol = omegaL();
				double om = omegaM();

				InitializeSpectrum(&outputSpectrum);
				prop_second(D, &bField, &energyGrid, &energyWidth,
						&inputSpectrum, &outputSpectrum, dataPath, IRFlag, Zmax,
						RadioFlag, h, om, ol, 0);
				// add spectrum
				for (int i = 0; i < NUM_SPECIES; i++) {
					for (int j = 0; j < outputSpectrum.numberOfMainBins; j++)
						finalSpectrum.spectrum[i][j] +=
								outputSpectrum.spectrum[i][j];
				}

			}
		}

		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}

	for (int j = 0; j < finalSpectrum.numberOfMainBins; j++) {
		outfile << (energyGrid.vector[j] / EeV * (eV * ELECTRON_MASS)) << " ";
		for (int i = 0; i < NUM_SPECIES; i++) {
			outfile << finalSpectrum.spectrum[i][j] << " ";
		}
		outfile << "\n";
	}

	DeleteSpectrum(&finalSpectrum);
	DeleteSpectrum(&outputSpectrum);
	DeleteSpectrum(&inputSpectrum);
	Delete_dCVector(&bField);

	Delete_dCVector(&energyGrid);
	Delete_dCVector(&energyWidth);

}

} // namespace crpropa
