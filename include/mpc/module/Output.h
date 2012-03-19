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
	std::ofstream outfile;

public:
	TrajectoryOutput(std::string name);
	~TrajectoryOutput();
	void process(Candidate *candidate);
	std::string getDescription() const;
};

/**
 @class FlaggedOutput
 @brief Saves particles with a given flag to a CSV file.
 */
class FlaggedOutput: public Module {
private:
	mutable std::ofstream outfile;
	Candidate::Status flag;

public:
	FlaggedOutput(std::string name, Candidate::Status flag);
	~FlaggedOutput();
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
