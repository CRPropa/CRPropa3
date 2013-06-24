#ifndef CRPROPA_PHOTONELECA_H
#define CRPROPA_PHOTONELECA_H

#include "crpropa/Module.h"
#include "crpropa/magneticField/MagneticField.h"

namespace crpropa {

class PhotonEleCa: public Module {
private:
	ref_ptr<MagneticField> field;

public:
	PhotonEleCa(ref_ptr<MagneticField> field);
	~PhotonEleCa();
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

} // namespace crpropa

#endif // CRPROPA_PHOTONELECA_H
