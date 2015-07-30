#ifndef CRPROPA_OUTPUTTXT_H
#define CRPROPA_OUTPUTTXT_H

#include "crpropa/Module.h"
#include "crpropa/AssocVector.h"

#include <fstream>
#include <bitset>

namespace crpropa {

enum OutputColumn {
	TrajectoryLengthColumn,
	RedshiftColumn,
	CurrentIdColumn,
	CurrentEnergyColumn,
	CurrentPositionColumn,
	CurrentDirectionColumn,
	SourceIdColumn,
	SourceEnergyColumn,
	SourcePositionColumn,
	SourceDirectionColumn,
	CreatedIdColumn,
	CreatedEnergyColumn,
	CreatedPositionColumn,
	CreatedDirectionColumn
};

class TextOutput: public Module {
protected:
	double lengthScale, energyScale;
	std::ostream *out;
	std::ofstream outfile;
	std::bitset<64> fields;
	std::string filename;
	bool oneDimensional;
public:
	TextOutput();
	~TextOutput();
	TextOutput(std::ostream &out);
	TextOutput(const std::string &filename);
	void setEnergyScale(double scale);
	void setLengthScale(double scale);

	void set(OutputColumn field, bool value);
	void enable(OutputColumn field);
	void disable(OutputColumn field);
	void enableAll();
	void disableAll();
	void resetFile();
	void printHeader();
	void set1D(bool value);
	void process(Candidate *candidate) const;
	void endRun();

	void gzip();
};

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
	void endRun();
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
	void endRun();
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
	void endRun();
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
	void endRun();
};

} // namespace crpropa

#endif // CRPROPA_OUTPUTTXT_H
