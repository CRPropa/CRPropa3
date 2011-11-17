#ifndef MAGNETICFIELD_H_
#define MAGNETICFIELD_H_

#include "ThreeVector.h"
#include "gadget/MagneticField.hpp"
#include "gadget/SmoothParticle.hpp"
//#include "SystemOfUnits.h"
//#include <math.h>
//#include <vector>
#include <memory>
//#include <stdint.h>
//#include <sstream>
//#include <unistd.h>

namespace mpc{

class MagneticField {
public:
	virtual Hep3Vector getField(const Hep3Vector &position) const = 0;
};

class HomogeneousMagneticField: public MagneticField {
public:
	HomogeneousMagneticField(const Hep3Vector &value) :
			value(value) {
	}
	Hep3Vector getField(const Hep3Vector &position) const {
		return value;
	}
private:
	Hep3Vector value;
};

class SampledSphMagneticField: public MagneticField {
public:
	SampledSphMagneticField(Hep3Vector &originKpc, double &sizeKpc,
			double &gridStepSizeKpc, std::string &sphFileName) {
		this->originKpc = originKpc;
		this->sizeKpc = sizeKpc;
		this->gridStepSizeKpc = gridStepSizeKpc;
		this->sphFileName;
		initialize();
	}
	~SampledSphMagneticField() {
	}
	Hep3Vector getField(Hep3Vector position) const {
		Vector3f r = Vector3f(position.x(), position.y(), position.z());
		Vector3f b = sampledMagneticField->getField(r / kpc);
		Hep3Vector bField = Hep3Vector(b.x, b.y, b.z) * gauss;
		return bField;
	}
private:
	Hep3Vector originKpc;
	double sizeKpc;
	double gridStepSizeKpc;
	std::string sphFileName;
	std::auto_ptr<SampledMagneticField> sampledMagneticField;

	void initialize() {
		std::cout << "[SampledSphMagneticField] initializing field from "
				<< sphFileName << std::endl;
		Vector3f orig = Vector3f(originKpc.x(), originKpc.y(), originKpc.z());
		sampledMagneticField.reset(new SampledMagneticField(orig, sizeKpc));
		sampledMagneticField.init(gridStepSizeKpc);
		sampledMagneticField.load(sphFileName);
	}
};

}

#endif /* MAGNETICFIELD_H_ */
