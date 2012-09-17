#ifndef OUTPUT_H_
#define OUTPUT_H_

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
 @brief Saves trajectories to CSV file.
 */
class TrajectoryOutput: public Module {
	mutable std::ofstream outfile;
public:
	TrajectoryOutput(std::string name);
	~TrajectoryOutput();
	void process(Candidate *candidate) const;
};

/**
 @class ConditionalOutput
 @brief Saves particles with a given property to a CSV file.
 */
class ConditionalOutput: public Module {
	mutable std::ofstream outfile;
	std::string propertyName;
	bool removeProperty;
public:
	ConditionalOutput(std::string filename, std::string propName = "Detected", bool removeProperty = false);
	~ConditionalOutput();
	void setRemoveProperty(bool removeProperty);
	void process(Candidate *candidate) const;
};

/**
 @class CRPropa2EventOutput
 @brief Saves events in CRPropa2 format.
 */
class CRPropa2EventOutput: public Module {
	mutable std::ofstream outfile;
public:
	CRPropa2EventOutput(std::string name);
	~CRPropa2EventOutput();
	void process(Candidate *candidate) const;
};

/**
 @class CRPropa2TrajectoryOutput
 @brief Saves trajectory in CRPropa2 format.
 */
class CRPropa2TrajectoryOutput: public Module {
	mutable std::ofstream outfile;
public:
	CRPropa2TrajectoryOutput(std::string name);
	~CRPropa2TrajectoryOutput();
	void process(Candidate *candidate) const;
};

} // namespace mpc

#endif /* OUTPUT_H_ */
