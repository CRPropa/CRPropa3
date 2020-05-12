#include "crpropa/module/PropagationCK.h"

#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace crpropa {

// Cash-Karp coefficients
const double cash_karp_a[] = {
	0., 0., 0., 0., 0., 0.,
	1. / 5., 0., 0., 0., 0., 0.,
	3. / 40., 9. / 40., 0., 0., 0., 0.,
	3. / 10., -9. / 10., 6. / 5., 0., 0., 0.,
	-11. / 54., 5. / 2., -70. / 27., 35. / 27., 0., 0.,
	1631. / 55296., 175. / 512., 575. / 13824., 44275. / 110592., 253. / 4096., 0.
};

const double cash_karp_b[] = {
	37. / 378., 0, 250. / 621., 125. / 594., 0., 512. / 1771.
};

const double cash_karp_bs[] = {
	2825. / 27648., 0., 18575. / 48384., 13525. / 55296., 277. / 14336., 1. / 4.
};

void PropagationCK::tryStep(const Y &y, Y &out, Y &error, double h,
		ParticleState &particle, double z) const {
	std::vector<Y> k;
	k.reserve(6);

	out = y;
	error = Y(0);

	// calculate the sum of b_i * k_i
	for (size_t i = 0; i < 6; i++) {

		Y y_n = y;
		for (size_t j = 0; j < i; j++)
			y_n += k[j] * a[i * 6 + j] * h;

		// update k_i
		k[i] = dYdt(y_n, particle, z);

		out += k[i] * b[i] * h;
		error += k[i] * (b[i] - bs[i]) * h;
	}
}

PropagationCK::Y PropagationCK::dYdt(const Y &y, ParticleState &p, double z) const {
	// normalize direction vector to prevent numerical losses
	Vector3d velocity = y.u.getUnitVector() * c_light;
	Vector3d B(0, 0, 0);
	try {
		B = field->getField(y.x, z);
	} catch (std::exception &e) {
		std::cerr << "PropagationCK: Exception in getField." << std::endl;
		std::cerr << e.what() << std::endl;
	}
	// Lorentz force: du/dt = q*c/E * (v x B)
	Vector3d dudt = p.getCharge() * c_light / p.getEnergy() * velocity.cross(B);
	return Y(velocity, dudt);
}

PropagationCK::PropagationCK(ref_ptr<MagneticField> field, double tolerance,
		double minStep, double maxStep) :
		minStep(0) {
	setField(field);
	setTolerance(tolerance);
	setMaximumStep(maxStep);
	setMinimumStep(minStep);

	// load Cash-Karp coefficients
	a.assign(cash_karp_a, cash_karp_a + 36);
	b.assign(cash_karp_b, cash_karp_b + 6);
	bs.assign(cash_karp_bs, cash_karp_bs + 6);
}

void PropagationCK::process(Candidate *candidate) const {
	// save the new previous particle state
	ParticleState &current = candidate->current;
	candidate->previous = current;

	double step = clip(candidate->getNextStep(), minStep, maxStep);

	// rectilinear propagation for neutral particles
	if (current.getCharge() == 0) {
		Vector3d pos = current.getPosition();
		Vector3d dir = current.getDirection();
		current.setPosition(pos + dir * step);
		candidate->setCurrentStep(step);
		candidate->setNextStep(maxStep);
		return;
	}

	Y yIn(current.getPosition(), current.getDirection());
	Y yOut, yErr;
	double newStep = step;
	double r = 42;  // arbitrary value > 1
	double z = candidate->getRedshift();

	// try performing step until the target error (tolerance) or the minimum step size has been reached
	while (r > 1) {
		step = newStep;
		tryStep(yIn, yOut, yErr, step / c_light, current, z);

		r = yErr.u.getR() / tolerance;  // ratio of absolute direction error and tolerance
		newStep = step * 0.95 * pow(r, -0.2);  // update step size to keep error close to tolerance
		newStep = clip(newStep, 0.1 * step, 5 * step);  // limit the step size change
		newStep = clip(newStep, minStep, maxStep);

		if (step == minStep)
			break;  // performed step already at the minimum
	}

	current.setPosition(yOut.x);
	current.setDirection(yOut.u.getUnitVector());
	candidate->setCurrentStep(step);
	candidate->setNextStep(newStep);
}

void PropagationCK::setField(ref_ptr<MagneticField> f) {
	field = f;
}

void PropagationCK::setTolerance(double tol) {
	if ((tol > 1) or (tol < 0))
		throw std::runtime_error(
				"PropagationCK: target error not in range 0-1");
	tolerance = tol;
}

void PropagationCK::setMinimumStep(double min) {
	if (min < 0)
		throw std::runtime_error("PropagationCK: minStep < 0 ");
	if (min > maxStep)
		throw std::runtime_error("PropagationCK: minStep > maxStep");
	minStep = min;
}

void PropagationCK::setMaximumStep(double max) {
	if (max < minStep)
		throw std::runtime_error("PropagationCK: maxStep < minStep");
	maxStep = max;
}

double PropagationCK::getTolerance() const {
	return tolerance;
}

double PropagationCK::getMinimumStep() const {
	return minStep;
}

double PropagationCK::getMaximumStep() const {
	return maxStep;
}

std::string PropagationCK::getDescription() const {
	std::stringstream s;
	s << "Propagation in magnetic fields using the Cash-Karp method.";
	s << " Target error: " << tolerance;
	s << ", Minimum Step: " << minStep / kpc << " kpc";
	s << ", Maximum Step: " << maxStep / kpc << " kpc";
	return s.str();
}

} // namespace crpropa
