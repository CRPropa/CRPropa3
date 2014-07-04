#include "crpropa/PhotonPropagation.h"
#include "crpropa/Common.h"

#include "EleCa/Propagation.h"
#include "EleCa/Particle.h"
#include "EleCa/Common.h"

#include <fstream>
#include <stdexcept>
#include <limits>
#include <iostream>

namespace crpropa {

void EleCaPropagation(const std::string &background,
		const std::string &inputfile, std::vector<double> &energy,
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

void EleCaPropagation(const std::string &background,
		const std::string &inputfile, const std::string &outputfile) {
	std::vector<double> energy, spectrum;
	EleCaPropagation(background, inputfile, energy, spectrum);
	std::ofstream output(outputfile.c_str());
	output << "# E N\n";
	for (size_t i = 0; i < energy.size(); i++) {
		output << energy[i] << " " << spectrum[i] << "\n";
	}
}

} // namespace crpropa

