#ifndef MPC_BREAKCONDITION_H_
#define MPC_BREAKCONDITION_H_

#include "mpc/Module.h"

namespace mpc {

/**
 @class MaximumTrajectoryLength
 @brief Stops propagation after a maximum trajectory length.

 This modules deactivates the candidate at a given maximum trajectory length.
 In that case the property ("Deactivated", module::description) is set.
 It also limits the candidates next step size to ensure the maximum trajectory length is no exceeded.
 */
class MaximumTrajectoryLength: public Module {
private:
	double maxLength;
	void updateDescription();

public:
	MaximumTrajectoryLength(double length = 0);
	void setMaximumTrajectoryLength(double length);
	double getMaximumTrajectoryLength() const;
	void process(Candidate *candidate) const;
};

/**
 @class MinimumEnergy
 @brief Stops propagation if energy drops under the set minimum energy.

 This modules deactivates the candidate below a give minimum energy.
 In that case the property ("Deactivated", module::description) is set.
 */
class MinimumEnergy: public Module {
private:
	double minEnergy;
	void updateDescription();

public:
	MinimumEnergy(double minEnergy = 0);
	void setMinimumEnergy(double energy);
	double getMinimumEnergy() const;
	void process(Candidate *candidate) const;
};

} // namespace mpc

#endif /* MPC_BREAKCONDITION_H_ */
