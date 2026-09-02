#ifndef SIMPLEPROPAGATION_H
#define SIMPLEPROPAGATION_H

#include "crpropa/Module.h"
#include "crpropa/Units.h"

namespace crpropa {
/**
 * \addtogroup Propagation 
 * @{
 */

/**
 @class SimplePropagation
 @brief Simple rectilinear propagation in absence of magnetic fields.

 This module implements rectilinear propagation.
 The step size is guaranteed to be larger than minStep and smaller than maxStep.
 It always proposes a next step size of maxStep.
 */
class SimplePropagation: public Module {
private:
	double minStep, maxStep;

public:
	/** Constructor
	 * @param minStep	minimum stepsize in [m], is converted to [s] over 1/c_light
	 * @param maxStep	maximum stepsize in [m], is converted to [s] over 1/c_light
	 */
	SimplePropagation(double minStep = (1 * kiloparsec), double maxStep = (1 * gigaparsec));
	/** Constructor
	 * @param useTimeStep   bool to differenciate from other constructor
	 * @param minTimeStep	minimum timestep in [s]
	 * @param maxTimeStep	maximum timestep in [s]
	 */
	SimplePropagation(bool useTimeSteps, double minStep, double maxStep);
	void process(Candidate *candidate) const;
	void setMinimumStep(double minStep);
	void setMaximumStep(double maxStep);
	void setMinimumTimeStep(double minStep);
	void setMaximumTimeStep(double maxStep);
	inline double getMinimumStep() const {return minStep*c_light;}
	inline double getMaximumStep() const {return maxStep*c_light;}
	inline double getMinimumTimeStep() const {return minStep;}
	inline double getMaximumTimeStep() const {return maxStep;}
	std::string getDescription() const;
};
/** @}*/

} // namespace crpropa

#endif // SIMPLEPROPAGATION_H

