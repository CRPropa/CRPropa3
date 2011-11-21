#ifndef SPHMAGNETICFIELD_H_
#define SPHMAGNETICFIELD_H_

#include "mpc/Vector3.h"
#include "mpc/Units.h"
#include "mpc/MagneticField.h"

#include "gadget/MagneticField.hpp"
#include "gadget/SmoothParticle.hpp"

#include <memory>

namespace mpc {

class SphMagneticField: public MagneticField {
public:
	SphMagneticField(Vector3 &originKpc, double &sizeKpc,
			double &gridStepSizeKpc, std::string &sphFileName) {
		this->originKpc = originKpc;
		this->sizeKpc = sizeKpc;
		this->gridStepSizeKpc = gridStepSizeKpc;
		this->sphFileName;
		initialize();
	}
	~SphMagneticField() {
	}
	Vector3 getField(Vector3 position) const {
		Vector3f r = Vector3f(position.x(), position.y(), position.z());
		Vector3f b = magneticField->getField(r / kpc);
		Vector3 bField = Vector3(b.x, b.y, b.z) * gauss;
		return bField;
	}
private:
	Vector3 originKpc;
	double sizeKpc;
	double gridStepSizeKpc;
	std::string sphFileName;
	std::auto_ptr< ::MagneticField > magneticField;

	void initialize() {
		std::cout << "[SampledSphMagneticField] initializing field from "
				<< sphFileName << std::endl;
		Vector3f orig = Vector3f(originKpc.x(), originKpc.y(), originKpc.z());

		SampledMagneticField * sampled = new SampledMagneticField(orig,
				sizeKpc);
		sampled->init(gridStepSizeKpc);
		sampled->load(sphFileName);
		magneticField.reset(sampled);
	}
};

} // namespace mpc

#endif /* SPHMAGNETICFIELD_H_ */
