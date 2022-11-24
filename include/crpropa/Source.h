#ifndef CRPROPA_SOURCE_H
#define CRPROPA_SOURCE_H

#include "crpropa/Candidate.h"
#include "crpropa/Grid.h"
#include "crpropa/EmissionMap.h"


#include <vector>

namespace crpropa {
/** @addtogroup SourceFeatures
 *  @{
 */


/**
 @class SourceFeature
 @brief Abstract base class for specific source features
 */
class SourceFeature: public Referenced {
protected:
	std::string description;
public:
	virtual void prepareParticle(ParticleState& particle) const {};
	virtual void prepareCandidate(Candidate& candidate) const;
	std::string getDescription() const;
};


/**
 @class SourceInterface
 @brief Abstract base class for sources
 */
class SourceInterface : public Referenced {
public:
	virtual ref_ptr<Candidate> getCandidate() const = 0;
	virtual std::string getDescription() const = 0;
};


/**
 @class Source
 @brief General source of particles

 This class is a container for source features.
 The source prepares a new candidate by passing it to all its source features
 to be modified accordingly.
 */
class Source: public SourceInterface {
	std::vector<ref_ptr<SourceFeature> > features;
public:
	void add(SourceFeature* feature);
	ref_ptr<Candidate> getCandidate() const;
	std::string getDescription() const;
};


/**
 @class SourceList
 @brief List of particle sources of individual luminosities.

 The SourceList is a source itself. It can be used if several sources are
 needed in one simulation.
 */
class SourceList: public SourceInterface {
	std::vector<ref_ptr<Source> > sources;
	std::vector<double> cdf;
public:
	/** Add an individual source to the list.
	 @param source		source to be added
	 @param weight		weight of the source; defaults to 1.
	 */
	void add(Source* source, double weight = 1);
	ref_ptr<Candidate> getCandidate() const;
	std::string getDescription() const;
};


/**
 @class SourceParticleType
 @brief Particle type at the source

 This feature assigns a single particle type to the source.
 For multiple types, use e.g. SourceMultipleParticleTypes.
 Particles are identified following the PDG numbering scheme:
   https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
 */
class SourceParticleType: public SourceFeature {
	int id;
public:
	/** Constructor for a source with a sign
	 @param id		id of the particle following the PDG numbering scheme
	*/
	SourceParticleType(int id);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};


/**
 @class SourceMultipleParticleTypes
 @brief Multiple particle types with individual relative abundances

 This feature assigns particle types to the events emitted by the sources.
 It is possible to control the relative abundance of each particle species.
 Particles are identified following the PDG numbering scheme:
   https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
 */
class SourceMultipleParticleTypes: public SourceFeature {
	std::vector<int> particleTypes;
	std::vector<double> cdf;
public:
	/** Constructor
	 */
	SourceMultipleParticleTypes();
	/** Add an individual particle type.
	 @param id			id of the particle following the PDG numbering scheme
	 @param weight		relative abundance of individual particle species
	 */
	void add(int id, double weight = 1);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};


/**
 @class SourceEnergy
 @brief Sets the initial energy of the emitted particles to a specific value

 This feature assigns a monochromatic spectrum, i.e., a single energy to all particles.
 */
class SourceEnergy: public SourceFeature {
	double E;
public:
	/** Constructor
	 @param energy		energy of the particle (in Joules)
	 */
	SourceEnergy(double energy);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};


/**
 @class SourcePowerLawSpectrum
 @brief Particle energy following a power-law spectrum

 The power law is of the form: dN/dE ~ E^index, for energies in the interval [Emin, Emax].
 */
class SourcePowerLawSpectrum: public SourceFeature {
	double Emin;
	double Emax;
	double index;
public:
	/** Constructor
	 @param Emin		minimum energy (in Joules)
	 @param Emax		maximum energy (in Joules)
	 @param index		spectral index of the power law
	 */
	SourcePowerLawSpectrum(double Emin, double Emax, double index);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};


/**
 @class SourceComposition
 @brief Multiple nuclei species with a rigidity-dependent power-law spectrum

 The power law is of the form: E^index, for energies in the interval [Emin, Z * Rmax].
 */
class SourceComposition: public SourceFeature {
	double Emin;
	double Rmax;
	double index;
	std::vector<int> nuclei;
	std::vector<double> cdf;
public:
	/** Constructor
	 @param Emin		minimum energy (in Joules)
	 @param Rmax		maximum rigidity (in Volts)
	 @param index		spectral index of the power law
	 */
	SourceComposition(double Emin, double Rmax, double index);
	/** Add individual particle species with a given abundance
	 @param id			id of the particle following the PDG numbering scheme
	 @param abundance	relative abundance of the particle species
	 */
	void add(int id, double abundance);
	/** Add individual particle species with a given abundance
	 @param A			atomic mass of the cosmic-ray nucleus
	 @param Z			atomic number of the cosmic-ray nucleus
	 @param abundance	relative abundance of the particle species
	 */
	void add(int A, int Z, double abundance);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};


/**
 @class SourcePosition
 @brief Position of a point source
 */
class SourcePosition: public SourceFeature {
	Vector3d position; 
public:
	/** Constructor for a source in 3D
	 @param position	vector containing the coordinates of the point source [in meters]
	 */
	SourcePosition(Vector3d position);
	/** Constructor for a source in 1D
	 @param d	distance of the point source to the observer at x = 0 [in meters]; 
	 			internally this will be converted to a vector with x-coordinate equal to d
	 */
	SourcePosition(double d);
	void prepareParticle(ParticleState &state) const;
	void setDescription();
};


/**
 @class SourceMultiplePositions
 @brief Multiple point-source positions with individual luminosities
 */
class SourceMultiplePositions: public SourceFeature {
	std::vector<Vector3d> positions;
	std::vector<double> cdf;
public:
	/** Constructor.
	 The sources must be added individually to the object.
	 */
	SourceMultiplePositions();
	/** Add an individual source with a given luminosity/contribution.
	 @param position	vector containing the coordinates of the point source [in meters]
	 @param weight		luminosity/contribution of the individual source
	 */
	void add(Vector3d position, double weight = 1);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};


/**
 @class SourceUniformSphere
 @brief Uniform distribution of sources in a spherical volume
 */
class SourceUniformSphere: public SourceFeature {
	Vector3d center;
	double radius;
public:
	/** Constructor
	 @param center		vector containing the coordinates of the center of the sphere
	 @param radius		radius of the sphere
	 */
	SourceUniformSphere(Vector3d center, double radius);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};


/**
 @class SourceUniformHollowSphere
 @brief Uniform distribution of sources between two spheres
 */
class SourceUniformHollowSphere: public SourceFeature {
	Vector3d center;
	double radius_inner;
	double radius_outer;
public:
	/** Constructor
	 @param center			vector containing the coordinates of the center of the sphere
	 @param radius_inner	radius of the inner sphere
	 @param radius_outer	radius of the outer sphere
	 */
	SourceUniformHollowSphere(Vector3d center,
			double radius_inner, double double_outer);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};


/**
 @class SourceUniformShell
 @brief Uniform distribution of source positions on the surface of a sphere
 */
class SourceUniformShell: public SourceFeature {
	Vector3d center;
	double radius;
public:
	/** Constructor
	 @param center		vector containing the coordinates of the center of the sphere
	 @param radius		radius of the sphere
	 */
	SourceUniformShell(Vector3d center, double radius);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};


/**
 @class SourceUniformBox
 @brief Uniform random source positions inside a box. The box is aligned with the coordinate axes.
 */
class SourceUniformBox: public SourceFeature {
	Vector3d origin;	// lower box corner
	Vector3d size;		// sizes along each coordinate axes.
public:
	/** Constructor
	 @param origin	vector corresponding to the lower box corner
	 @param size	vector corresponding to the box sizes along each direction
	 */
	SourceUniformBox(Vector3d origin, Vector3d size);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};


/**
 @class SourceUniformCylinder
 @brief Uniform distribution of source positions inside the volume of a cylinder. 

 The circle of the cylinder lays in the xy-plane and the height is along the z-axis.
 */
class SourceUniformCylinder: public SourceFeature {
	Vector3d origin;	// central point of cylinder 
	double height;		// total height of the cylinder along z-axis. Half over/under the center.
	double radius;		// radius of the cylinder in the xy-plane
public:
	/** Constructor
	 @param origin	vector corresponding to the center of the cylinder axis
	 @param height	height of the cylinder, half lays over the origin, half is lower
	 @param radius	radius of the cylinder
	 */
	SourceUniformCylinder(Vector3d origin, double height, double radius);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};


/**
 @class SourceSNRDistribution
 @brief Source distribution that follows the Galactic SNR distribution in 2D

 The origin of the distribution is the Galactic center. The default maximum radius is set 
 to rMax=20 kpc and the default maximum height is zMax = 5 kpc.
 See G. Case and D. Bhattacharya (1996) for the details of the distribution.
 */
class SourceSNRDistribution: public SourceFeature {
	double rEarth; // parameter given by observation
	double alpha; // parameter to shift the maximum in R direction
	double beta; // parameter to shift the maximum in R direction
	double zg; // exponential cut parameter in z direction
	double frMax; // helper for efficient sampling
	double fzMax; // helper for efficient sampling
	double rMax; // maximum radial distance - default 20 kpc 
		      // (due to the extension of the JF12 field)
	double zMax; // maximum distance from galactic plane - default 5 kpc
	void setFrMax(); // calculate frMax with the current parameter. 

public:
	/** Default constructor. 
	 Default parameters are:
	 . rEarth = 8.5 kpc
	 . alpha = 2
	 . beta = 3.53
	 . zg = 300 pc
	 . rMax = 20 kpc
	 . zMax = 5 kpc
	*/ 
	SourceSNRDistribution();
	/** Generic constructor
	 @param rEarth	  distance from Earth to the Galactic centre [in meters]
	 @param alpha	  parameter that shifts radially the maximum of the distributions
	 @param beta	  parameter that shifts radially the maximum of the distributions 
	 @param zg		  exponential cut-off parameter in the z-direction [in meters]
	*/	
	SourceSNRDistribution(double rEarth,double alpha, double beta, double zg);

	void prepareParticle(ParticleState &particle) const;
	/**
	 radial distribution of the SNR density.
	 @param r	galactocentric radius in [meter]
	*/
	double fr(double r) const;
	/**
	 height distribution of the SNR density.
	 @param z	height over/under the galactic plane in [meter]
	*/
	double fz(double z) const;

	/**
	 Set the exponential cut-off parameter in the z-direction.
	 @param Zg	cut-off parameter
	*/
	void setFzMax(double Zg);

	/**
	 @param rMax maximal radius up to which sources are possible
	*/
	void setRMax(double rMax);

	/**
	 @param zMax maximal height up to which sources are possible
	*/
	void setZMax(double zMax);

	// parameter for the raidal distribution
	void setAlpha(double a);
	// parameter for the exponential cut-off in the radial distribution
	void setBeta(double b);
	double getFrMax() const;
	double getFzMax() const;
	double getRMax() const;
	double getZMax() const;
	double getAlpha() const;
	double getBeta() const;
	void setDescription();
};


/**
 @class SourcePulsarDistribution
 @brief Source distribution following the Galactic pulsar distribution

 A logarithmic spiral with four arms is used for the radial distribution.
 The z-distribution is a simple exponentially decaying distribution.
 The pulsar distribution is explained in detail in C.-A. Faucher-Giguere
 and V. M. Kaspi, ApJ 643 (May, 2006) 332. The radial distribution is 
 parametrized as in Blasi and Amato, JCAP 1 (Jan., 2012) 10.
 */
class SourcePulsarDistribution: public SourceFeature {
	double rEarth; // parameter given by observation
	double beta; // parameter to shift the maximum in R direction
	double zg; // exponential cut parameter in z direction
	double frMax; // helper for efficient sampling
	double fzMax; // helper for efficient sampling
	double rMax; // maximum radial distance - default 22 kpc 
	double zMax; // maximum distance from galactic plane - default 5 kpc
	double rBlur; // relative smearing factor for the radius
	double thetaBlur; // smearing factor for the angle. Unit = [1/length]
public:
	/** Default constructor. 
	 Default parameters are:
	 . rEarth = 8.5 kpc
	 . beta = 3.53
	 . zg = 300 pc
	 . Rmax = 22 kpc
	 . Zmax = 5 kpc
	 . rBlur = 0.07
	 . thetaBlur = 0.35 / kpc
	 */ 
	SourcePulsarDistribution();	
	/** Generic constructor
	 @param rEarth		distance from Earth to the Galactic centre [in meters]
	 @param beta		parameter that shifts radially the maximum of the distributions 
	 @param zg			exponential cut-off parameter in the z-direction [in meters]
	 @param rBlur		relative smearing factor for radius
	 @param thetaBlur	smearing factor for the angle [in 1 / meters]
	 */	
	SourcePulsarDistribution(double rEarth, double beta, double zg, double rBlur, double thetaBlur);
	void prepareParticle(ParticleState &particle) const;

	/** 
	 radial distribution of pulsars
	 @param r	galactocentric radius
	*/
	double fr(double r) const;
	/**
	 z distribution of pulsars
	 @param z	height over/under the galactic plane
	*/
	double fz(double z) const;
	double ftheta(int i, double r) const;
	double blurR(double r_tilde) const;
	double blurTheta(double theta_tilde, double r_tilde) const;
	void setFrMax(double R, double b);
	void setFzMax(double zg);
	void setRMax(double rMax);
	void setZMax(double zMax);
	void setRBlur(double rBlur);
	void setThetaBlur(double thetaBlur);
	double getFrMax();
	double getFzMax();
	double getRMax();
	double getZMax();
	double getRBlur();
	double getThetaBlur();
	void setDescription();
};


/**
 @class SourceUniform1D
 @brief Uniform source distribution in 1D

 This source property sets random x-coordinates according to a uniform source
 distribution in a given distance interval. If cosmological effects are included, 
 this is done by drawing a light-travel distance from a flat distribution and
 converting to a comoving distance. In the absence of cosmological effects, the
 positions are drawn uniformly in the light-travel distance interval (as opposed
 to a comoving interval).
 The source positions are assigned to the x-coordinate (Vector3d(distance, 0, 0))
 in this one-dimensional case.
 */
class SourceUniform1D: public SourceFeature {
	double minD; // minimum light-travel distance
	double maxD; // maximum light-travel distance
	bool withCosmology;	// whether to account for cosmological effects (expansion of the Universe)
public:
	/** Constructor
	 @param minD			minimum distance; comoving if withCosmology is True
	 @param maxD 			maximum distance; comoving if withCosmology is True
	 @param withCosmology	whether to account for cosmological effects (expansion of the Universe)
	 */
	SourceUniform1D(double minD, double maxD, bool withCosmology = true);
	void prepareParticle(ParticleState& particle) const;
	void setDescription();
};


/**
 @class SourceDensityGrid
 @brief Random source positions from a density grid
 */
class SourceDensityGrid: public SourceFeature {
	ref_ptr<Grid1f> grid;
public:
	/** Constructor
	 @param densityGrid 	3D grid containing the density of sources in each cell
	 */
	SourceDensityGrid(ref_ptr<Grid1f> densityGrid);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};


/**
 @class SourceDensityGrid1D
 @brief Random source positions from a 1D density grid
 */
class SourceDensityGrid1D: public SourceFeature {
	ref_ptr<Grid1f> grid;	// 1D grid with Ny = Nz = 1
public:
	/** Constructor
	 @param densityGrid 	1D grid containing the density of sources in each cell, Ny and Nz must be 1
	 */
	SourceDensityGrid1D(ref_ptr<Grid1f> densityGrid);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};


/**
 @class SourceIsotropicEmission
 @brief Isotropic emission from a source
 */
class SourceIsotropicEmission: public SourceFeature {
public:
	/** Constructor
	 */
	SourceIsotropicEmission();
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};


/**
 @class SourceDirectedEmission
 @brief Directed emission from a source from the von-Mises-Fisher distribution 
 
 The emission from the source is generated following the von-Mises-Fisher distribution
 with mean direction mu and concentration parameter kappa.
 The sampling from the vMF distribution follows this document by Julian Straub:
 http://people.csail.mit.edu/jstraub/download/straub2017vonMisesFisherInference.pdf
 The emitted particles are assigned a weight so that the detected particles can be
 reweighted to an isotropic emission distribution instead of a vMF distribution.
 For details, see PoS (ICRC2019) 447.
 */
class SourceDirectedEmission: public SourceFeature {
	Vector3d mu; // Mean emission direction in the vMF distribution
	double kappa; // Concentration parameter of the vMF distribution
	double ca; // helpers for the efficient calculation of frame rotation
	double sa;
	double cd;
	double sd;
public:
	/** Constructor
	 @param mu	mean direction of the emission, mu should be normelized
	 @param kappa	concentration parameter
	*/
	SourceDirectedEmission(Vector3d mu, double kappa);
	void prepareCandidate(Candidate &candidate) const;
	/**
	 set sampling parameter Ca
	 @param alpha	angle between x and y component of direction. alpha = arctan(mu.y / mu.x)
	*/
	void setCa(double alpha);
	/**
	 set sampling parameter Sa
	 @param alpha	angle between x and y component of direction. alpha = arctan(mu.y / mu.x)
	*/
	void setSa(double alpha);
	/**
	 set sampling parameter Cd
	 @param delta	angle between mu vector and z-axis. delta = arcsin(mu.z) 
	*/
	void setCd(double delta);
	/**
	 set sampling parameter Sd
	 @param delta	angle between mu vector and z-axis. delta = arcsin(mu.z) 
	*/
	void setSd(double delta);
	double getCa() const;
	double getSa() const;
	double getCd() const;
	double getSd() const;
	void setDescription();
};

/**
 @class SourceLambertDistributionOnSphere
 @brief Uniform random position on a sphere with isotropic Lamberts distributed directions.

 This function should be used for crosschecking the arrival distribution for a
 Galactic propagation with an isotropic arrival distribution at the Edge of our
 Galaxy. Note, that for simulation speed you should rather use the backtracking
 technique: see e.g. http://physik.rwth-aachen.de/parsec
 */
class SourceLambertDistributionOnSphere: public SourceFeature {
	Vector3d center;	// center of the sphere
	double radius;		// radius of the sphere
	bool inward;		// if true, direction point inwards
public:
	/** Constructor
	 @param center		vector containing the coordinates of the center of the sphere
	 @param radius		radius of the sphere
	 @param inward		if true, the directions point inwards
	 */
	SourceLambertDistributionOnSphere(const Vector3d &center, double radius, bool inward);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};


/**
 @class SourceDirection
 @brief Collimated emission along a specific direction
 */
class SourceDirection: public SourceFeature {
	Vector3d direction;
public:
	/** Constructor
	 @param direction	Vector3d corresponding to the direction of emission
	 */
	SourceDirection(Vector3d direction = Vector3d(-1, 0, 0));
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};


/**
 @class SourceEmissionMap
 @brief Deactivate Candidate if it has zero probability in provided EmissionMap. 

	This feature does not change the direction of the candidate. Therefore a usefull direction feature (isotropic or directed emission)
	must be added to the sources before. The propability of the emission map is not taken into account. 
 */
class SourceEmissionMap: public SourceFeature {
	ref_ptr<EmissionMap> emissionMap;
public:
	/** Constructor
	 @param emissionMap		emission map containing probabilities of emission in various directions
	 */
	SourceEmissionMap(EmissionMap *emissionMap);
	void prepareCandidate(Candidate &candidate) const;
	void setEmissionMap(EmissionMap *emissionMap);
	void setDescription();
};


/**
 @class SourceEmissionCone
 @brief Uniform emission within a cone
 */
class SourceEmissionCone: public SourceFeature {
	Vector3d direction;
	double aperture;
public:
	/** Constructor
	 @param direction		Vector3d corresponding to the cone axis 
	 @param aperture		opening angle of the cone
	 */
	SourceEmissionCone(Vector3d direction, double aperture);
	void prepareParticle(ParticleState &particle) const;

	/**
	 @param direction Vector3d corresponding to the cone axis
	*/
	void setDirection(Vector3d direction);
	void setDescription();
};


/**
 @class SourceRedshift
 @brief Emission of particles at a specific redshift (or time)

 The redshift coordinate is used to treat cosmological effects and as a time coordinate.
 Consider, for instance, a source located at a distance corresponding to a redshift z. 
 In the absence of processes that cause time delays (e.g., magnetic deflections), particles
 from this source could arrive after a time corresponding to the source redshift. Charged 
 particles, on the other hand, can arrive at a time later than the corresponding straight-
 line travel duration. 
 This treatment is also useful for time-dependent studies (e.g. transient sources).
 */
class SourceRedshift: public SourceFeature {
	double z;
public:
	/** Constructor
	 @param z		redshift of emission
	 */
	SourceRedshift(double z);
	void prepareCandidate(Candidate &candidate) const;
	void setDescription();
};


/**
 @class SourceUniformRedshift
 @brief Random redshift (time of emission) from a uniform distribution

 This function assigns random redshifts to the particles emitted by a given source.
 These values are drawn from a uniform redshift distribution in the interval [zmin, zmax].
 The redshift coordinate is used to treat cosmological effects and as a time coordinate.
 Consider, for instance, a source located at a distance corresponding to a redshift z. 
 In the absence of processes that cause time delays (e.g., magnetic deflections), particles
 from this source could arrive after a time corresponding to the source redshift. Charged 
 particles, on the other hand, can arrive at a time later than the corresponding straight-
 line travel duration. 
 This treatment is also useful for time-dependent studies (e.g. transient sources).
 */
class SourceUniformRedshift: public SourceFeature {
	double zmin, zmax;
public:
	/** Constructor
	 @param zmin	minimum redshift
	 @param zmax	maximum redshift
	 */
	SourceUniformRedshift(double zmin, double zmax);
	void prepareCandidate(Candidate &candidate) const;
	void setDescription();
};


/**
 @class SourceRedshiftEvolution
 @brief Random redshift (time of emission) from (1+z)^m distribution

 This assigns redshifts to a given source according to a typical power-law distribution.
 The redshift coordinate is used to treat cosmological effects and as a time coordinate.
 Consider, for instance, a source located at a distance corresponding to a redshift z. 
 In the absence of processes that cause time delays (e.g., magnetic deflections), particles
 from this source could arrive after a time corresponding to the source redshift. Charged 
 particles, on the other hand, can arrive at a time later than the corresponding straight-
 line travel duration. 
 This treatment is also useful for time-dependent studies (e.g. transient sources).
 */
class SourceRedshiftEvolution: public SourceFeature {
	double zmin, zmax;
	double m;
public:
	/** Constructor
	 @param m		index of the power law (1 + z)^m
	 @param zmin	minimum redshift
	 @param zmax	maximum redshift
	 */
	SourceRedshiftEvolution(double m, double zmin, double zmax);
	void prepareCandidate(Candidate &candidate) const;
};


/**
 @class SourceRedshift1D
 @brief Redshift according to the distance to 0

 This source property sets the redshift according to the distance from 
 the source to the origin (0, 0, 0). 
 It must be added after the position of the source is set because it
 computes the redshifts based on the source distance.
 */
class SourceRedshift1D: public SourceFeature {
public:
	/** Constructor
	 */
	SourceRedshift1D();
	void prepareCandidate(Candidate &candidate) const;
	void setDescription();
};


#ifdef CRPROPA_HAVE_MUPARSER
/**
 @class SourceGenericComposition
 @brief Add multiple cosmic rays with energies described by an expression string

 This is particularly useful if an arbitrary combination of nuclei types with 
 specific energy spectra. The strings parsed may contain 'A' (atomic mass), 
 'Z' (atomic number).  The following units are recognized as part of the strings:
 GeV, TeV, PeV, EeV.  The variable for energy is 'E', with limits 'Emin', 'Emax'.
 This property only works if muparser is available.
 For details about the library see:
 	https://beltoforion.de/en/muparser/
 */
class SourceGenericComposition: public SourceFeature {
public:
	struct Nucleus {
		int id;
		std::vector<double> cdf;
	};
	/** Constructor
	 @param Emin		minimum energy [in Joules]
	 @param Emax		maximum energy [in Joules]
	 @param expression	string containing the expression to generate the composition
	 @param bins		number of energy bins
	 */
	SourceGenericComposition(double Emin, double Emax, std::string expression, size_t bins = 1024);
	/** Add an individual particle id.
	 @param id			id of the particle following the PDG numbering scheme
	 @param abundance	relative abundance of individual particle species
	 */
	void add(int id, double abundance);
	/** Add an individual particle id.
	 @param A			atomic mass of the cosmic-ray nucleus
	 @param Z			atomic number of the cosmic-ray nucleus
	 @param weight		relative abundance of individual particle species
	 */
	void add(int A, int Z, double abundance);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();

	const std::vector<double> *getNucleusCDF(int id) const {
		for (size_t i = 0; i < nuclei.size(); i++) {
			if (nuclei[i].id == id)
				return &nuclei[i].cdf;
		}
		return 0;
	}

protected:
	double Emin, Emax;
	size_t bins;
	std::string expression;
	std::vector<double> energy;

	std::vector<Nucleus> nuclei;
	std::vector<double> cdf;

};
#endif

/**  @} */ // end of group SourceFeature

}// namespace crpropa

#endif // CRPROPA_SOURCE_H
