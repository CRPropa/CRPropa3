#include "crpropa/module/PropagationBP.h"

#include <sstream>
#include <stdexcept>
#include <vector>

namespace crpropa {
	void PropagationBP::tryStep(const Y &y, Y &out, Y &error, double h,
			ParticleState &particle, double z) const {
		Y outHelp = dY(y.x, y.u, h/2, z, particle);  // 2 steps with h/2
		Y outCompare = dY(outHelp.x, outHelp.u, h/2, z, particle);

		out = dY(y.x, y.u, h, z, particle);  // 1 step with h

		error = errorEstimation(out.x , outCompare.x , h*particle.getVelocity().getR());
	}

	PropagationBP::Y PropagationBP::dY(Vector3d pos, Vector3d dir, double dt, 
		double z, ParticleState &current) const {

		Vector3d vel = dir*current.getVelocity().getR();
			
		// do nothing if velocity is zero to avoid dividing by zero in getUnitVector
		if (vel==Vector3(0))
			return Y(pos, dir);
			
		Vector3d B = getFieldAtPosition(pos, z);
		double q = current.getCharge();
		// lorentz factor is between 1 and infinity (but never actually infinity)
		double gamma = current.getLorentzFactor();
		double m = gamma*current.getMass();

		// do a half leapfrog step
		pos += vel * dt / 2.;

		// Boris help vectors:
		Vector3d t = B * q / 2. / m * dt;
		Vector3d s = t * 2. / (1. + t.dot(t));

		Vector3d v_prime = vel + vel.cross(t);
		vel = vel + v_prime.cross(s);  // final velocity

		// the other half leapfrog step (only happens if vel!=0)
		pos += vel * dt / 2.;

		return Y(pos, vel.getUnitVector());
	}

	PropagationBP::PropagationBP(ref_ptr<MagneticField> BField, double fixedStep) :
		minStep(0) {
		setField(BField);
		setTolerance(0.42);
		setMaximumTimeStep(fixedStep/c_light);
		setMinimumTimeStep(fixedStep/c_light);
	}

	PropagationBP::PropagationBP(ref_ptr<MagneticField> BField, double tolerance, double minStep, double maxStep) :
		minStep(0) {
		setField(BField);
		setTolerance(tolerance);
		setMaximumTimeStep(maxStep/c_light);
		setMinimumTimeStep(minStep/c_light);
	}

	PropagationBP::PropagationBP(double fixedTimeStep, ref_ptr<MagneticField> BField) :
		minStep(0) {
		setField(BField);
		setTolerance(0.42);
		setMaximumTimeStep(fixedTimeStep);
		setMinimumTimeStep(fixedTimeStep);
	}

	PropagationBP::PropagationBP(double tolerance, double minTimeStep, double maxTimeStep, ref_ptr<MagneticField> BField) :
		minStep(0) {
		setField(BField);
		setTolerance(tolerance);
		setMaximumTimeStep(maxTimeStep);
		setMinimumTimeStep(minTimeStep);
	}


	void PropagationBP::process(Candidate *candidate) const {
		// save the new previous particle state
		ParticleState &current = candidate->current;
		candidate->previous = current;

		Y yIn(current.getPosition(), current.getDirection());

		// calculate charge of particle
		double step = maxStep;

		// rectilinear propagation for neutral particles
		if (current.getCharge() == 0) {
			step = clip(candidate->getNextStep(), minStep, maxStep);
			current.setPosition(yIn.x + yIn.u * candidate->getVelocity() * step);
			candidate->setCurrentStep(step);
			candidate->setNextStep(maxStep);
			return;
		}

		Y yOut, yErr;
		double newStep = step;
		double z = candidate->getRedshift();

		// if minStep is the same as maxStep the adaptive algorithm with its error
		// estimation is not needed and the computation time can be saved:
		if (minStep == maxStep){
			yOut = dY(yIn.x, yIn.u, step, z, current);
		} else {
			step = clip(candidate->getNextStep(), minStep, maxStep);
			newStep = step;
			double r = 42;  // arbitrary value

			// try performing step until the target error (tolerance) or the minimum/maximum step size has been reached
			while (true) {
				tryStep(yIn, yOut, yErr, step, current, z);
				r = yErr.u.getR() / tolerance;  // ratio of absolute direction error and tolerance
				if (r > 1) {  // large direction error relative to tolerance, try to decrease step size
					if (step == minStep)  // already minimum step size
						break;
					else {
						newStep = step * 0.95 * pow(r, -0.2);
						newStep = std::max(newStep, 0.1 * step); // limit step size decrease
						newStep = std::max(newStep, minStep); // limit step size to minStep
						step = newStep;
					}
				} else {  // small direction error relative to tolerance, try to increase step size
					if (step != maxStep) {  // only update once if maximum step size yet not reached
						newStep = step * 0.95 * pow(r, -0.2);
						newStep = std::min(newStep, 5 * step); // limit step size increase
						newStep = std::min(newStep, maxStep); // limit step size to maxStep
					}
					break;
				}
			}
		}

		current.setPosition(yOut.x);
		current.setDirection(yOut.u);
		candidate->setCurrentStep(step);
		candidate->setNextStep(newStep);
	}


	void PropagationBP::setField(ref_ptr<MagneticField> f) {
		field = f;
	}


	ref_ptr<MagneticField> PropagationBP::getField() const {
		return field;
	}


	Vector3d PropagationBP::getFieldAtPosition(Vector3d pos, double z) const {
		Vector3d B(0, 0, 0);
		try {
			// check if field is valid and use the field vector at the
			// position pos with the redshift z
			if (field.valid())
				B = field->getField(pos, z);
		} catch (std::exception &e) {
			KISS_LOG_ERROR 	<< "PropagationBP: Exception in PropagationBP::getFieldAtPosition.\n"
					<< e.what();
		}	
		return B;
	}


	double PropagationBP::errorEstimation(const Vector3d x1, const Vector3d x2, double step) const {
		// compare the position after one step with the position after two steps with step/2.
		// 1/4 = (1/2)²  number of steps for x1 divided by number of steps for x2 to the power of p (order)
		return (x1 - x2).getR() / (step * (1 - 1/4.) );
	}

	void PropagationBP::setTolerance(double tol) {
		if ((tol > 1) or (tol < 0))
			throw std::runtime_error(
					"PropagationBP: target error not in range 0-1");
		tolerance = tol;
	}

	void PropagationBP::setMinimumStep(double min) {
		if (min < 0)
			throw std::runtime_error("PropagationBP: minStep < 0 ");
		if (min/c_light > maxStep)
			throw std::runtime_error("PropagationBP: minStep > maxStep");
		minStep = min/c_light;
	}

	void PropagationBP::setMaximumStep(double max) {
		if (max/c_light < minStep)
			throw std::runtime_error("PropagationBP: maxStep < minStep");
		maxStep = max/c_light;
	}

	void PropagationBP::setMinimumTimeStep(double min) {
		if (min < 0)
			throw std::runtime_error("PropagationBP: minStep < 0 ");
		if (min > maxStep)
			throw std::runtime_error("PropagationBP: minStep > maxStep");
		minStep = min;
	}

	void PropagationBP::setMaximumTimeStep(double max) {
		if (max < minStep)
			throw std::runtime_error("PropagationBP: maxStep < minStep");
		maxStep = max;
	}

	std::string PropagationBP::getDescription() const {
		std::stringstream s;
		s << "Propagation in magnetic fields using the adaptive Boris push method.";
		s << " Target error: " << tolerance;
		s << ", Minimum Step: " << minStep / kiloyear << " kiloyear";
		s << ", Maximum Step: " << maxStep / kiloyear << " kiloyear";
		return s.str();
	}
} // namespace crpropa
