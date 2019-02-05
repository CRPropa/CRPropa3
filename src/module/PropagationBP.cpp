#include "crpropa/module/PropagationBP.h"

#include <sstream>
#include <stdexcept>
#include <vector>

namespace crpropa {

	PropagationBP::PropagationBP(ref_ptr<MagneticField> field,
								 double propStep) : propStep(1 * Mpc)
	{
		setField(field);
		setStep(propStep);
	}


	void PropagationBP::process(Candidate *candidate) const {
		// save the new previous particle state
		ParticleState &current = candidate->current;
		candidate->previous = current;

		double step = getStep();
		candidate->setCurrentStep(step);
		double h = step / c_light;

		//get particle properties
		double q = current.getCharge();
		double m = current.getEnergy()/(c_light * c_light);
		double z = candidate->getRedshift();
		Vector3d pos = current.getPosition();
		Vector3d dir = current.getDirection();

		// half leap frog step in the position
		pos += c_light * dir * h / 2.;

		// get B field at particle position
		Vector3d B(0, 0, 0);
		try {
			B = field->getField(pos, z);
		} catch (std::exception &e) {
			std::cerr << "PropagationBP: Exception in getField." << std::endl;
			std::cerr << e.what() << std::endl;
		}
		// Boris help vectors
		Vector3d t = B * q / 2 / m * h;
		Vector3d s = t * 2 / (1 + t.dot(t));
		Vector3d v_help;

		// Boris push
		v_help = dir + dir.cross(t);
		dir = dir + v_help.cross(s);

		// full leap frog step in the velocity
		candidate->current.setDirection(dir);

		// the other half leap frog step in the position
		pos += c_light * dir * h / 2. ;
		candidate->current.setPosition(pos);
		candidate->setCurrentStep(h);
		candidate->setNextStep(h);
	}


	void PropagationBP::setField(ref_ptr<MagneticField> f) {
		field = f;
	}

	void PropagationBP::setStep(double step) {
		if (step < 0.)
			throw std::runtime_error("PropagationBP: step < 0");
		propStep = step;
	}

	double PropagationBP::getStep() const {
		return propStep;
	}

	std::string PropagationBP::getDescription() const {
		std::stringstream s;
		s << "Propagation in magnetic fields using the Boris push method.";
		s << ", step: " << propStep / kpc << " kpc";
		return s.str();
	}

} // namespace crpropa
