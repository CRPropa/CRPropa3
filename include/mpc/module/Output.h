#ifndef MPC_OUTPUT_H_
#define MPC_OUTPUT_H_

#include "mpc/Module.h"
#include <fstream>

namespace mpc {

/**
 @class ShellOutput
 @brief Write properties of the candidate to the shell.
 */
class ShellOutput: public Module {
public:
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

/**
 @class TrajectoryOutput
 @brief Saves trajectories to plain text file.
 */
class TrajectoryOutput: public Module {
	mutable std::ofstream outfile;
public:
	TrajectoryOutput(std::string filename);
	~TrajectoryOutput();
	void process(Candidate *candidate) const;
};

/**
 @class ConditionalOutput
 @brief Saves particles with a given property to a plain text file.
 */
class ConditionalOutput: public Module {
	mutable std::ofstream outfile;
	std::string condition;
	bool removeProperty;
public:
	ConditionalOutput(std::string filename, std::string condition = "Detected",
			bool removeProperty = true);
	~ConditionalOutput();
	void process(Candidate *candidate) const;
	void setRemoveProperty(bool removeProperty);
};

/**
 @class TrajectoryOutput1D
 @brief Saves 1D trajectories to plain text file.
 */
class TrajectoryOutput1D: public Module {
	mutable std::ofstream outfile;
public:
	TrajectoryOutput1D(std::string filename);
	~TrajectoryOutput1D();
	void process(Candidate *candidate) const;
};

/**
 @class EventOutput1D
 @brief Records particles that are inactive and have the property 'Detected' to a plain text file.
 */
class EventOutput1D: public Module {
	mutable std::ofstream outfile;
public:
	EventOutput1D(std::string filename);
	~EventOutput1D();
	void process(Candidate *candidate) const;
};

} // namespace mpc

#endif // MPC_OUTPUT_H_
