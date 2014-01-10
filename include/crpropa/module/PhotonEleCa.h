#ifndef CRPROPA_PHOTONELECA_H
#define CRPROPA_PHOTONELECA_H

#include "crpropa/Module.h"
#include "crpropa/magneticField/MagneticField.h"

#include <memory>
#include <fstream>

// forward declaration
namespace eleca {
class Propagation;
}

namespace crpropa {

class PhotonEleCa: public Module {
private:
	std::auto_ptr<eleca::Propagation> propagation;
	mutable std::ofstream output;
public:
	PhotonEleCa(const std::string background, const std::string &filename);
	~PhotonEleCa();
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

} // namespace crpropa

#endif // CRPROPA_PHOTONELECA_H
