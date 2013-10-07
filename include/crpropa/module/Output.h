#ifndef CRPROPA_OUTPUT_H
#define CRPROPA_OUTPUT_H

#include "crpropa/Module.h"
#include "crpropa/AssocVector.h"
#include <fstream>

namespace crpropa {

/**
 @class TrajectoryOutput
 @brief Saves trajectories to plain text file.
 */
class TrajectoryOutput: public Module {
	mutable std::ofstream fout;
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
	mutable std::ofstream fout;
	std::string condition;
public:
	ConditionalOutput(std::string filename, std::string condition = "Detected");
	~ConditionalOutput();
	void process(Candidate *candidate) const;
};

/**
 @class TrajectoryOutput1D
 @brief Saves 1D trajectories to plain text file.
 */
class TrajectoryOutput1D: public Module {
	mutable std::ofstream fout;
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
	mutable std::ofstream fout;
public:
	EventOutput1D(std::string filename);
	~EventOutput1D();
	void process(Candidate *candidate) const;
};

} // namespace crpropa

#endif // CRPROPA_OUTPUT_H
