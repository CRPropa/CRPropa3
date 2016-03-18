#ifndef CRPROPA_OUTPUTCRPROPA2_H
#define CRPROPA_OUTPUTCRPROPA2_H

#include "crpropa/Module.h"
#include <fstream>

namespace crpropa {

/**
 @class CRPropa2EventOutput3D
 @brief Saves events to a plain text file in CRPropa2 format.
 */
class CRPropa2EventOutput3D: public Module {
	mutable std::ofstream outfile;
public:
	CRPropa2EventOutput3D(std::string filename);
	~CRPropa2EventOutput3D();
	void process(Candidate *candidate) const;
	void close();
};

/**
 @class CRPropa2TrajectoryOutput3D
 @brief Saves trajectories to a plain text file in CRPropa2 format.
 */
class CRPropa2TrajectoryOutput3D: public Module {
	mutable std::ofstream outfile;
public:
	CRPropa2TrajectoryOutput3D(std::string filename);
	~CRPropa2TrajectoryOutput3D();
	void process(Candidate *candidate) const;
	void close();
};

/**
 @class CRPropa2EventOutput1D
 @brief Saves 1D events to a plain text file in CRPropa2 format.
 */
class CRPropa2EventOutput1D: public Module {
	mutable std::ofstream outfile;
public:
	CRPropa2EventOutput1D(std::string filename);
	~CRPropa2EventOutput1D();
	void process(Candidate *candidate) const;
	void close();
};

/**
 @class CRPropa2TrajectoryOutput1D
 @brief Saves 1D trajectories to a plain text file in CRPropa2 format.
 */
class CRPropa2TrajectoryOutput1D: public Module {
	mutable std::ofstream outfile;
public:
	CRPropa2TrajectoryOutput1D(std::string filename);
	~CRPropa2TrajectoryOutput1D();
	void process(Candidate *candidate) const;
	void close();
};

} // namespace crpropa

#endif // CRPROPA_OUTPUTCRPROPA2_H
