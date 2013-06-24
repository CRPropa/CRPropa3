#ifndef CRPROPA_BREAKCONDITION_H
#define CRPROPA_BREAKCONDITION_H

#include "crpropa/Module.h"

namespace crpropa {

/**
 @class MaximumTrajectoryLength
 @brief Deactivates the candidate beyond a maximum trajectory length

 This modules deactivates the candidate at a given maximum trajectory length.
 In that case the property ("Deactivated", module::description) is set.
 It also limits the candidates next step size to ensure the maximum trajectory length is no exceeded.
 */
class MaximumTrajectoryLength: public Module {
	double maxLength;
	std::string flag;
public:
	MaximumTrajectoryLength(double length = 0,
			std::string flag = "Deactivated");
	void setMaximumTrajectoryLength(double length);
	double getMaximumTrajectoryLength() const;
	void setFlag(std::string flag);
	std::string getFlag() const;
	std::string getDescription() const;
	void process(Candidate *candidate) const;
};

/**
 @class MinimumEnergy
 @brief Deactivates the candidate below a minimum energy

 This modules deactivates the candidate below a given minimum energy.
 In that case the property ("Deactivated", module::description) is set.
 */
class MinimumEnergy: public Module {
	double minEnergy;
	std::string flag;
public:
	MinimumEnergy(double minEnergy = 0, std::string flag = "Deactivated");
	void setMinimumEnergy(double energy);
	double getMinimumEnergy() const;
	void setFlag(std::string flag);
	std::string getFlag() const;
	std::string getDescription() const;
	void process(Candidate *candidate) const;
};

/**
 @class MinimumRedshift
 @brief Deactivates the candidate below a minimum redshift

 This modules deactivates the candidate below a given minimum redshift.
 In that case the property ("Deactivated", module::description) is set.
 */
class MinimumRedshift: public Module {
	double zmin;
	std::string flag;
public:
	MinimumRedshift(double zmin = 0, std::string flag = "Deactivated");
	void setMinimumRedshift(double z);
	double getMinimumRedshift();
	void setFlag(std::string flag);
	std::string getFlag() const;
	std::string getDescription() const;
	void process(Candidate *candidate) const;
};

} // namespace crpropa

#endif // CRPROPA_BREAKCONDITION_H
