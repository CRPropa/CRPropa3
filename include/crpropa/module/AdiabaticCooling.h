#pragma once
#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>

#include "crpropa/Vector3.h"
#include "crpropa/Module.h"
#include "crpropa/Units.h"
#include "crpropa/advectionField/AdvectionField.h"

namespace crpropa {

/**
@class AdiabaticCooling
@brief Implements adiabatic cooling/heating due to advection.
*/

class AdiabaticCooling: public Module {
	private:
		ref_ptr<AdvectionField> advectionField;
		double limit;

	public:
		AdiabaticCooling(ref_ptr<AdvectionField> advectionField);
		AdiabaticCooling(ref_ptr<AdvectionField> advectionField, double limit);		
		void process(Candidate *c) const;
	
		void setLimit(double l);

		double getLimit() const;

};	



}; // end namesspace crpropa
