#include "crpropa/module/PhotonEleCa.h"

namespace crpropa {

PhotonEleCa::PhotonEleCa(ref_ptr<MagneticField> field) {
}

PhotonEleCa::~PhotonEleCa() {
}

void PhotonEleCa::process(Candidate *candidate) const {
	if (candidate->current.getId() != 22)
		return; // do nothing if not a photon

	return;
}

std::string PhotonEleCa::getDescription() const {
	std::stringstream s;
	s << "PhotonEleCa";
	return s.str();
}

} // namespace crpropa
