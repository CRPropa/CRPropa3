#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "mpc/Module.h"
#include <fstream>

namespace mpc {

/**
 @class TrajectoryOutput
 @brief Saves trajectories to CSV file.
 */
class TrajectoryOutput: public Module {
private:
	mutable std::ofstream outfile;

public:
	TrajectoryOutput(std::string name);
	~TrajectoryOutput();
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

/**
 @class ConditionalOutput
 @brief Saves particles with a given property to a CSV file.
 */
class ConditionalOutput: public Module {
private:
	mutable std::ofstream outfile;
	std::string propertyName;
	bool removeProperty;
public:
	ConditionalOutput(std::string filename, std::string propName);
	~ConditionalOutput();
	void setRemoveProperty(bool removeProperty);
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

/**
 @class ShellOutput
 @brief Write properties of the candidate to the shell.
 */
class ShellOutput: public Module {
public:
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

} // namespace mpc

#endif /* OUTPUT_H_ */
