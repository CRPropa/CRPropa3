#include "crpropa/module/PhotonEleCa.h"

#include "EleCa/Propagation.h"
#include "EleCa/Particle.h"
#include "EleCa/Common.h"

#include <vector>

namespace crpropa {

PhotonEleCa::PhotonEleCa() :
		propagation(new eleca::Propagation) {
	propagation->ReadTables(getDataPath("eleca_lee.txt"));
	propagation->InitBkgArray("CMB");
}

PhotonEleCa::~PhotonEleCa() {
}

void PhotonEleCa::process(Candidate *candidate) const {
	if (candidate->current.getId() != 22)
		return; // do nothing if not a photon

	eleca::Particle p0(candidate->current.getId(),
			candidate->current.getEnergy() / eV, candidate->getRedshift());

	std::vector<eleca::Particle> ParticleAtMatrix;
	std::vector<eleca::Particle> ParticleAtGround;
	ParticleAtMatrix.push_back(p0);

	while (ParticleAtMatrix.size() > 0) {

		eleca::Particle p1 = ParticleAtMatrix.back();
		ParticleAtMatrix.pop_back();

		if (p1.IsGood()) {
			propagation->Propagate(p1, ParticleAtMatrix, ParticleAtGround);
		}
	}

#pragma omp critical
	{
		propagation->WriteOutput(p0, ParticleAtGround, 0);
	}

	candidate->setActive(false);
	return;
}

std::string PhotonEleCa::getDescription() const {
	std::stringstream s;
	s << "PhotonEleCa";
	return s.str();
}

} // namespace crpropa
