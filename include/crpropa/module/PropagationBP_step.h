#ifndef CRPROPA_PROPAGATIONBP_STEP_H
#define CRPROPA_PROPAGATIONBP_STEP_H

#include "crpropa/Module.h"
#include "crpropa/Units.h"
#include "crpropa/magneticField/MagneticField.h"

namespace crpropa {
/**
 * \addtogroup Propagation
 * @{
 */

/**
 @class PropagationBP
 @brief Rectilinear propagation through magnetic fields using the Boris method.

 This module solves the equations of motion of a relativistic charged particle when propagating through a magnetic field.\n
 It uses the Boris push integration method.\n
 The step size control tries to keep the relative error close to, but smaller than the designated tolerance.
 Additionally a minimum and maximum size for the steps can be set.
 For neutral particles a rectilinear propagation is applied and a next step of the maximum step size proposed.
 */
    class PropagationBP_step: public Module {
    public:
        class Y {
        public:
            Vector3d x, u; /*< phase-point: position and direction */

            Y() {
            }

            Y(const Vector3d &x, const Vector3d &u) :
                    x(x), u(u) {
            }

            Y(double f) :
                    x(Vector3d(f, f, f)), u(Vector3d(f, f, f)) {
            }

            Y operator *(double f) const {
                return Y(x * f, u * f);
            }

            Y &operator +=(const Y &y) {
                x += y.x;
                u += y.u;
                return *this;
            }
        };

    private:
        std::vector<double> a, b, bs; /*< Cash-Karp coefficients */
        ref_ptr<MagneticField> field;
        double tolerance; /*< target relative error of the numerical integration */
        double minStep; /*< minimum step size of the propagation */
        double maxStep; /*< maximum step size of the propagation */

    public:
        PropagationBP_step(ref_ptr<MagneticField> field = NULL, double tolerance = 1e-4,
                      double minStep = (0.1 * kpc), double maxStep = (1 * Gpc));
        void process(Candidate *candidate) const;

        Y dY(Vector3d  pos, Vector3d  dir, double step, double z, double q, double m) const;
        double errorEstimation(const Vector3d mu, const Vector3d muh, double h) const;

        void tryStep(const Y &y, Y &out, Y &error, double t,
                     ParticleState &p, double z, double m, double q) const;

        void setField(ref_ptr<MagneticField> field);
        void setTolerance(double tolerance);
        void setMinimumStep(double minStep);
        void setMaximumStep(double maxStep);

        ref_ptr<MagneticField> getField() const;
        double getTolerance() const;
        double getMinimumStep() const;
        double getMaximumStep() const;
        std::string getDescription() const;
    };
/** @}*/

} // namespace crpropa

#endif // CRPROPA_PROPAGATIONBP_STEP_H
