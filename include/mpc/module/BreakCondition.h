#ifndef MPC_BREAKCONDITION_H_
#define MPC_BREAKCONDITION_H_

#include "mpc/Module.h"

namespace mpc {

/**
 @class MaximumTrajectoryLength
 @brief Deactivates the candidate beyond a maximum trajectory length

 This modules deactivates the candidate at a given maximum trajectory length.
 In that case the property ("Deactivated", module::description) is set.
 It also limits the candidates next step size to ensure the maximum trajectory length is no exceeded.
 */
class MaximumTrajectoryLength: public Module {
private:
	double maxLength;

public:
	MaximumTrajectoryLength(double length = 0);
	void setMaximumTrajectoryLength(double length);
	double getMaximumTrajectoryLength() const;
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

/**
 @class MinimumEnergy
 @brief Deactivates the candidate below a minimum energy

 This modules deactivates the candidate below a given minimum energy.
 In that case the property ("Deactivated", module::description) is set.
 */
class MinimumEnergy: public Module {
private:
	double minEnergy;

public:
	MinimumEnergy(double minEnergy = 0);
	void setMinimumEnergy(double energy);
	double getMinimumEnergy() const;
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

/**
 @class MinimumRedshift
 @brief Deactivates the candidate below a minimum redshift

 This modules deactivates the candidate below a given minimum redshift.
 In that case the property ("Deactivated", module::description) is set.
 */
class MinimumRedshift: public Module {
private:
	double zmin;

public:
	MinimumRedshift(double zmin = 0);
	void setMinimumRedshift(double z);
	double getMinimumRedshift();
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

} // namespace mpc

#endif /* MPC_BREAKCONDITION_H_ */
