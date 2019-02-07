#include "crpropa/module/PropagationBP_step.h"

#include <sstream>
#include <stdexcept>
#include <vector>

namespace crpropa {


void PropagationBP_step::tryStep(const Y &y, Y &out, Y &error, double h,
                            ParticleState &particle, double z, double m, double q) const {
    out = dY(q, m, y.x, y.u, z, h);  // 1 step with h

    Y outHelp = dY(q, m, y.x, y.u, z, h/2);  // 2 steps with h/2
    Y outCompare = dY(q, m, outHelp.x, outHelp.u, z, h/2);

    error = errorEstimation(out.x , outCompare.x , h);
}


PropagationBP_step::Y PropagationBP_step::dY(double q, double m, Vector3d x, Vector3d v, double z, double step) const {

    double h = step;

    // half leap frog step in the position
    x += c_light * v * h /2. ;

    // get B field at particle position
    Vector3d B(0, 0, 0);
    try {
        B = field->getField(x, z);
    } catch (std::exception &e) {
        std::cerr << "PropagationBP_step: Exception in getField." << std::endl;
        std::cerr << e.what() << std::endl;
    }
    // Boris help vectors
    Vector3d t = B * q/2/m * h;
    Vector3d s = t *2. /(1+t.dot(t));
    Vector3d v_help;

    // Boris Push
    v_help = v + v.cross(t);
    v = v + v_help.cross(s);

    // the other half leap frog step in the position
    x += c_light * v * h /2. ;

    return Y(x, v);

}

PropagationBP_step::PropagationBP_step(ref_ptr<MagneticField> field, double tolerance,
                             double minStep, double maxStep) :
        minStep(0) {
    setField(field);
    setTolerance(tolerance);
    setMaximumStep(maxStep);
    setMinimumStep(minStep);
}

void PropagationBP_step::process(Candidate *candidate) const {
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
    // further particle parameters
    double m = current.getEnergy()/(c_light * c_light);
    double q = current.getCharge();


    // try performing step until the target error (tolerance) or the minimum step size has been reached
    while (true) {
        tryStep(yIn, yOut, yErr, step / c_light, current, z, m, q);
        r = yErr.u.getR() / tolerance;  // ratio of absolute direction error and tolerance
        if (r > 1) {  // large direction error relative to tolerance, try to decrease step size
            if (step == minStep)  // already minimum step size
                break;
            else {
                newStep = step * 0.95 * pow(r, -0.2);
                newStep = std::max(newStep, 0.1 * step); // limit step size decrease
                newStep = std::max(newStep, minStep); // limit step size to minStep
                step = newStep;
                std::cout << "error " << r << "  oldStep " << step/Mpc << "  newStep " << newStep/Mpc << std::endl;
            }
        } else {  // small direction error relative to tolerance, try to increase step size
            if (step != maxStep) {  // already maximum step size
                newStep = step * 0.95 * pow(r, -0.2);
                newStep = std::min(newStep, 5 * step); // limit step size increase
                newStep = std::min(newStep, maxStep); // limit step size to max Step
                std::cout << "error " << r << "  oldStep " << step/Mpc << "  newStep " << newStep/Mpc << std::endl;
            }
            break;
        }
    }

    current.setPosition(yOut.x);
    current.setDirection(yOut.u.getUnitVector());
    candidate->setCurrentStep(step);
    candidate->setNextStep(newStep);
}

void PropagationBP_step::setField(ref_ptr<MagneticField> f) {
    field = f;
}

double PropagationBP_step::errorEstimation(const Vector3d mu, const Vector3d muh, double h) const {

    Vector3d diff = (mu - muh);
    //~ std::cout << " diff = " << diff.getR() << std::endl;
    double S = diff.getR() / (h * c_light * (1 - 1/4 ) );    // 1/4 = (1/2)² = mu hoch p

    return S;
}

void PropagationBP_step::setTolerance(double tol) {
    if ((tol > 1) or (tol < 0))
        throw std::runtime_error(
                "PropagationBP: target error not in range 0-1");
    tolerance = tol;
}

void PropagationBP_step::setMinimumStep(double min) {
    if (min < 0)
        throw std::runtime_error("PropagationBP: minStep < 0 ");
    if (min > maxStep)
        throw std::runtime_error("PropagationBP: minStep > maxStep");
    minStep = min;
}

void PropagationBP_step::setMaximumStep(double max) {
    if (max < minStep)
        throw std::runtime_error("PropagationBP: maxStep < minStep");
    maxStep = max;
}

double PropagationBP_step::getTolerance() const {
    return tolerance;
}

double PropagationBP_step::getMinimumStep() const {
    return minStep;
}

double PropagationBP_step::getMaximumStep() const {
    return maxStep;
}

std::string PropagationBP_step::getDescription() const {
    std::stringstream s;
    s << "Propagation in magnetic fields using the adaptive Boris push method.";
    s << " Target error: " << tolerance;
    s << ", Minimum Step: " << minStep / kpc << " kpc";
    s << ", Maximum Step: " << maxStep / kpc << " kpc";
    return s.str();
}

} // namespace crpropa