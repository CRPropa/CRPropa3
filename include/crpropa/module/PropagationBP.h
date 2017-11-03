#ifndef CRPROPA_PROPAGATIONBP_H
#define CRPROPA_PROPAGATIONBP_H

#include "crpropa/Module.h"
#include "crpropa/Units.h"
#include "crpropa/magneticField/MagneticField.h"

namespace crpropa {
	
	class PropagationBP: public Module {

	ref_ptr<MagneticField> field;
	double Step;
	
	
public:
	PropagationBP(ref_ptr<MagneticField> field = NULL,
			double Step = (0.1 * Mpc));
			
	void process(Candidate *candidate) const;
	
	void setField(ref_ptr<MagneticField> field);
	void setStep(double Step);
	double getStep() const;
	
};
	
}// namespace crpropa

#endif // CRPROPA_PROPAGATIONBP_H
