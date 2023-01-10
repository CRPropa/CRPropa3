#ifndef CRPROPA_BREAKCONDITION_H
#define CRPROPA_BREAKCONDITION_H

#include "crpropa/Module.h"

namespace crpropa {
/**
 * \addtogroup Condition 
 * @{
 */

/**
 @class MaximumTrajectoryLength
 @brief Deactivates the candidate beyond a maximum trajectory length

 This module deactivates the candidate at a given maximum trajectory length.
 In that case the property ("Deactivated", module::description) is set.
 It also limits the candidates next step size to ensure the maximum trajectory length is not exceeded.
 */
class MaximumTrajectoryLength: public AbstractCondition {
	double maxLength;
	std::vector<Vector3d> observerPositions;
public:
	MaximumTrajectoryLength(double length = 0);
	void setMaximumTrajectoryLength(double length);
	double getMaximumTrajectoryLength() const;
	void addObserverPosition(const Vector3d &position);
	const std::vector<Vector3d>& getObserverPositions() const;
	std::string getDescription() const;
	void process(Candidate *candidate) const;
};

/**
 @class MinimumEnergy
 @brief Deactivates the candidate below a minimum energy

 This module deactivates the candidate below a given minimum energy.
 In that case the property ("Deactivated", module::description) is set.
 */
class MinimumEnergy: public AbstractCondition {
	double minEnergy;
public:
	MinimumEnergy(double minEnergy = 0);
	void setMinimumEnergy(double energy);
	double getMinimumEnergy() const;
	std::string getDescription() const;
	void process(Candidate *candidate) const;
};


/**
 @class MinimumRigidity
 @brief Deactivates the candidate below a minimum rigidity

 This module deactivates the candidate below a given minimum rigidity (E/Z in EeV).
 In that case the property ("Deactivated", module::description) is set.
 */
class MinimumRigidity: public AbstractCondition {
	double minRigidity;
public:
	MinimumRigidity(double minRigidity = 0);
	void setMinimumRigidity(double minRigidity);
	double getMinimumRigidity() const;
	std::string getDescription() const;
	void process(Candidate *candidate) const;
};

/**
 @class MinimumRedshift
 @brief Deactivates the candidate below a minimum redshift

 This module deactivates the candidate below a given minimum redshift.
 In that case the property ("Deactivated", module::description) is set.
 */
class MinimumRedshift: public AbstractCondition {
	double zmin;
public:
	MinimumRedshift(double zmin = 0);
	void setMinimumRedshift(double z);
	double getMinimumRedshift();
	std::string getDescription() const;
	void process(Candidate *candidate) const;
};

/**
 @class MinimumChargeNumber
 @brief Deactivates the candidate below a minimum number

 This module deactivates the candidate below a given minimum charge number.
 A minimum charge number of 26 deactivates all (anti-) isotopes which 
 are ranked in the periodic table before iron (Fe). 
 In that case the property ("Deactivated", module::description) is set.
 */
class MinimumChargeNumber: public AbstractCondition {
	int minChargeNumber;
public:
	MinimumChargeNumber(int minChargeNumber = 0);
	void setMinimumChargeNumber(int chargeNumber);
	int getMinimumChargeNumber() const;
	std::string getDescription() const;
	void process(Candidate *candidate) const;
};

/**
 @class MinimumEnergyPerParticleId
 @brief Deactivates the candidate below a minimum energy for specific particle Ids.

 This module deactivates the candidate below a given minimum energy for specific particle types.
 In that case the property ("Deactivated", module::description) is set.
 All particles whose minimum energy is not specified follow the more general minEnergyOthers condition.
 */
class MinimumEnergyPerParticleId: public AbstractCondition {
	std::vector<double> minEnergies;
	std::vector<int> particleIds;
	double minEnergyOthers;
public:
	MinimumEnergyPerParticleId(double minEnergyOthers = 0);
	void setMinimumEnergyOthers(double energy);
	double getMinimumEnergyOthers() const;
	void add(int id, double energy);
	std::string getDescription() const;
	void process(Candidate *candidate) const;
};


/**
 @class DetectionLength
 @brief Detects the candidate at a given trajectoryLength
 
 This break condition can be used for non-regular time observation of the particle density. See also ObserverTimeEvolution.
 */
class DetectionLength: public AbstractCondition {
	double detLength;
public:
	DetectionLength(double length = 0);
	void setDetectionLength(double length);
	double getDetectionLength() const;
	std::string getDescription() const;
	void process(Candidate *candidate) const;
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_BREAKCONDITION_H
