#ifndef CRPROPA_PROPAGATIONBP_STEP_H
#define CRPROPA_PROPAGATIONBP_STEP_H

#include "crpropa/Module.h"
#include "crpropa/Units.h"
#include "crpropa/magneticField/MagneticField.h"

namespace crpropa {
	
	class PropagationBP_step: public Module {

public:
	struct phasepoint{
		
		Vector3d x;
		Vector3d v;	
	};
	
private:

	ref_ptr<MagneticField> field;
	double minStep;
	double maxStep;
	double epsilon;
	
	
public:
	PropagationBP_step(ref_ptr<MagneticField> field = NULL,
			double minStep = (0.1 * kpc), double maxStep = (1. * Gpc), double epsilon = 0);
			
	void process(Candidate *candidate) const;
	
	phasepoint testStep(double q, double m, Vector3d x, Vector3d v, double z, double testStep) const;
	double errorEstimation(const Vector3d mu, const Vector3d muh, double h) const;
	
	void setField(ref_ptr<MagneticField> field);
	void setMinStep(double minStep);
	void setMaxStep(double maxStep);
	void setEpsilon(double epsilon);
	double getMinStep() const;
	double getMaxStep() const;
	double getEpsilon() const;
	
};
	
}// namespace crpropa

#endif // CRPROPA_PROPAGATIONBP_STEP_H
