#ifndef MPC_SOURCE_H
#define MPC_SOURCE_H

#include "mpc/Referenced.h"
#include "mpc/Candidate.h"
#include "mpc/Grid.h"
#include "mpc/module/Redshift.h"

#include <vector>

namespace mpc {

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
	SourcePowerLawSpectrum(double Emin, double Emax, double index);
	void prepare(ParticleState &particle) const;
};

/**
 @class SourceNuclei
 @brief Nuclei with given total abundances
 */
class SourceNuclei: public SourceProperty {
	std::vector<int> ids; /**< nucleus id */
	std::vector<double> abundances; /**< relative abundance of source isotopes at equal energies */
public:
	void add(int id, double abundance = 1);
	void prepare(ParticleState &particle) const;
};

/**
 @class SourceComposition
 @brief Nuclei with given abundances and a uniform power law spectrum from Emin to Z * Rmax
 */
class SourceComposition: public SourceProperty {
	double Emin;
	double Rmax;
	double index;
	std::vector<int> isotope; /**< isotope id */
	std::vector<double> abundance; /**< relative abundance of source isotopes at equal energies */
	std::vector<double> probability; /**< cumulative probability of source isotopes */
public:
	SourceComposition(double Emin, double Rmax, double index);
	double getSpectrumIntegral(int Z = 1) const;
	void add(int id, double abundance);
	void add(int A, int Z, double abundance);
	void normalize();
	void prepare(ParticleState &particle) const;
};

/**
 @class SourcePosition
 @brief Position of a point source
 */
class SourcePosition: public SourceProperty {
	Vector3d position;
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
 @class SourceUniformDistributionSphere
 @brief Uniform random source positions inside a sphere
 */
class SourceUniformDistributionSphere: public SourceProperty {
	Vector3d center;
	double radius;
public:
	SourceUniformDistributionSphere(Vector3d center, double radius);
	void prepare(ParticleState &particle) const;
};

/**
 @class SourceUniformDistributionBox
 @brief Uniform random source positions inside a box
 */
class SourceUniformDistributionBox: public SourceProperty {
	Vector3d origin;
	Vector3d size;
public:
	SourceUniformDistributionBox(Vector3d origin, Vector3d size);
	void prepare(ParticleState &particle) const;
};

/**
 @class SourceUniformDistribution1D
 @brief Uniform random source positions in for 1D simulations
 */
class SourceUniformDistribution1D: public SourceProperty {
	double minDistance, maxDistance;
public:
	SourceUniformDistribution1D(double minDistance, double maxDistance);
	void prepare(ParticleState& particle) const;
};

/**
 @class SourceDensityGrid
 @brief Provides source positions from a density grid

 This module takes a density grid to compute random initial positions.
 First a discrete bin is drawn following the density distribution.
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
	ref_ptr<Redshift> redshift;
public:
	SourceRedshift1D(Redshift* redshift);
	void prepare(Candidate &candidate) const;
};

}// namespace mpc

#endif /* MPC_SOURCE_H */
