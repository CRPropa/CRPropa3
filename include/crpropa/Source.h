#ifndef CRPROPA_SOURCE_H
#define CRPROPA_SOURCE_H

#include "crpropa/Referenced.h"
#include "crpropa/Candidate.h"
#include "crpropa/Grid.h"

#include <vector>

namespace crpropa {

/**
 @class SourceProperty
 @brief Abstract class for properties of cosmic ray sources
 */
class SourceProperty: public Referenced {
public:
	virtual void prepare(ParticleState& particle) const;
	virtual void prepare(Candidate& candidate) const;
};

/**
 @class Source
 @brief General cosmic ray source

 This class is a container for source properties.
 The source prepares a new candidate by passing it to all its source properties to be modified accordingly.
 */
class Source: public Referenced {
	std::vector<ref_ptr<SourceProperty> > properties;
public:
	void addProperty(SourceProperty* property);
	ref_ptr<Candidate> getCandidate() const;
};

/**
 @class SourceList
 @brief List of cosmic ray sources of individual total lumosities.

 The SourceList is a source itself. It can be used if several UHECR sources are needed in one simulation.
 */
class SourceList: public Source {
	std::vector<ref_ptr<Source> > sources;
	std::vector<double> luminosities;
public:
	void addSource(Source* source, double luminosity = 1);
	ref_ptr<Candidate> getCandidate() const;
};

/**
 @class SourceParticleType
 @brief Particle type at the source
 */
class SourceParticleType: public SourceProperty {
	double id;
public:
	SourceParticleType(int id);
	void prepare(ParticleState &particle) const;
};

/**
 @class SourceMultipleParticleTypes
 @brief Multiple particle types with individual (relative) total abundances
 */
class SourceMultipleParticleTypes: public SourceProperty {
	std::vector<int> ids; /**< particle ids */
	std::vector<double> abundances; /**< (relative) total abundances */
public:
	void add(int id, double abundance = 1);
	void prepare(ParticleState &particle) const;
};

/**
 @class SourceEnergy
 @brief Sets the initial energy to a given value
 */
class SourceEnergy: public SourceProperty {
	double E;
public:
	SourceEnergy(double energy);
	void prepare(ParticleState &particle) const;
};

/**
 @class SourcePowerLawSpectrum
 @brief Particle energy following a power law spectrum
 */
class SourcePowerLawSpectrum: public SourceProperty {
	double Emin;
	double Emax;
	double index;
public:
	/** Constructor
	 @param Emin	minimum energy
	 @param Emax	maximum energy
	 @param index	differential spectral index
	 */
	SourcePowerLawSpectrum(double Emin, double Emax, double differentialIndex);
	/** Set particle with a random energy from a power law distribution */
	void prepare(ParticleState &particle) const;
};

/**
 @class SourceComposition
 @brief Nuclei with given abundances and a uniform power law spectrum between Emin and Z * Rmax
 */
class SourceComposition: public SourceProperty {
	double Emin;
	double Rmax;
	double index;
	std::vector<int> isotope; /**< isotope id */
	std::vector<double> abundance; /**< relative abundance of source isotopes at equal energies */
	std::vector<double> probability; /**< cumulative probability of source isotopes */
	void normalize();
	double getSpectrumIntegral(int Z) const;

public:
	/** Constructor
	 @param Emin	minimum energy for cosmic rays
	 @param Rmax	maximum rigidity for cosmic rays
	 @param index	differential spectral index
	 */
	SourceComposition(double Emin, double Rmax, double index);
	/** Add a species to the composition
	 @param id 			particle id
	 @param abundance	absolute or relative abundance at a fixed value of energy/nucleons
	 */
	void add(int id, double abundance);
	/** Add a species to the composition
	 @param A 			mass number
	 @param Z			charge number
	 @param abundance	absolute or relative abundance at a fixed value of energy/nucleons
	 */
	void add(int A, int Z, double abundance);
	/** Randomly select a species and energy */
	void prepare(ParticleState &particle) const;
};

/**
 @class SourcePosition
 @brief Position of a point source
 */
class SourcePosition: public SourceProperty {
	Vector3d position; /**< Source position */
public:
	SourcePosition(Vector3d position);
	void prepare(ParticleState &state) const;
};

/**
 @class SourceMultiplePositions
 @brief Multiple point source positions with individual luminosities
 */
class SourceMultiplePositions: public SourceProperty {
	std::vector<Vector3d> positions;
	std::vector<double> luminosities;
public:
	void add(Vector3d position, double luminosity = 1);
	void prepare(ParticleState &particle) const;
};

/**
 @class SourceUniformSphere
 @brief Uniform random source positions inside a sphere
 */
class SourceUniformSphere: public SourceProperty {
	Vector3d center;
	double radius;
public:
	SourceUniformSphere(Vector3d center, double radius);
	void prepare(ParticleState &particle) const;
};

/**
 @class SourceUniformShell
 @brief Uniform random source positions on a sphere
 */
class SourceUniformShell: public SourceProperty {
	Vector3d center;
	double radius;
public:
	SourceUniformShell(Vector3d center, double radius);
	void prepare(ParticleState &particle) const;
};

/**
 @class SourceUniformBox
 @brief Uniform random source positions inside a box
 */
class SourceUniformBox: public SourceProperty {
	Vector3d origin;
	Vector3d size;
public:
	/** Constructor
	 @param origin	lower box corner
	 @param size	upper box corner
	 */
	SourceUniformBox(Vector3d origin, Vector3d size);
	void prepare(ParticleState &particle) const;
};

/**
 @class SourceUniform1D
 @brief 1D-Positions from a uniform source distribution in an expanding universe

 This source property sets random x-coordinates according to a uniform source
 distribution in a given comoving distance interval.
 This is done by drawing a light travel distance from a flat distribution and
 converting to a comoving distance.
 */
class SourceUniform1D: public SourceProperty {
	double minD; // minimum light travel distance
	double maxD; // maximum light travel distance
	bool withCosmology;
public:
	/** Constructor
	 @param minD	minimum comoving distance
	 @param maxD 	maximum comoving distance
	 @param withCosmology	specify if universe expanding
	 */
	SourceUniform1D(double minD, double maxD, bool withCosmology=true);
	void prepare(ParticleState& particle) const;
};

/**
 @class SourceDensityGrid
 @brief Random source positions from a density grid

 This source property takes a density grid to compute random initial positions.
 To dial a source position, first a bin is drawn following the density distribution.
 Then a random position is drawn from a uniform distribution in the bin.
 */
class SourceDensityGrid: public SourceProperty {
	ref_ptr<ScalarGrid> grid;
	float sumDensity;
public:
	SourceDensityGrid(ref_ptr<ScalarGrid> densityGrid);
	void prepare(ParticleState &particle) const;
};

/**
 @class SourceDensityGrid1D
 @brief Random source positions from a 1D density grid

 This source property takes a N*1*1 grid to compute random initial positions.
 To dial a source position, first a bin is drawn following the density distribution.
 Then a random position is drawn from a uniform distribution in the bin.
 */
class SourceDensityGrid1D: public SourceProperty {
	ref_ptr<ScalarGrid> grid;
	float sumDensity;
public:
	SourceDensityGrid1D(ref_ptr<ScalarGrid> densityGrid);
	void prepare(ParticleState &particle) const;
};

/**
 @class SourceIsotropicEmission
 @brief Isotropic emission from a source
 */
class SourceIsotropicEmission: public SourceProperty {
public:
	void prepare(ParticleState &particle) const;
};

/**
 @class SourceDirection
 @brief Emission in a discrete direction
 */
class SourceDirection: public SourceProperty {
	Vector3d direction;
public:
	SourceDirection(Vector3d direction = Vector3d(-1, 0, 0));
	void prepare(ParticleState &particle) const;
};

/**
 @class SourceEmissionCone
 @brief Uniform random emission inside a cone
 */
class SourceEmissionCone: public SourceProperty {
	Vector3d direction;
	double aperture;
public:
	SourceEmissionCone(Vector3d direction, double aperture);
	void prepare(ParticleState &particle) const;
};

/**
 @class SourceRedshift
 @brief Discrete redshift (time of emission)
 */
class SourceRedshift: public SourceProperty {
	double z;
public:
	SourceRedshift(double z);
	void prepare(Candidate &candidate) const;
};

/**
 @class SourceUniformRedshift
 @brief Uniform redshift distribution (time of emission)
 */
class SourceUniformRedshift: public SourceProperty {
	double zmin, zmax;
public:
	SourceUniformRedshift(double zmin, double zmax);
	void prepare(Candidate &candidate) const;
};

/**
 @class SourceRedshift1D
 @brief Redshift according to the distance to 0

 This source property sets the redshift according to the distance to 0.
 It must be added after a position setting source property.
 */
class SourceRedshift1D: public SourceProperty {
public:
	void prepare(Candidate &candidate) const;
};

}// namespace crpropa

#endif // CRPROPA_SOURCE_H
