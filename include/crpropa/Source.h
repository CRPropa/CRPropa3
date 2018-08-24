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
 @brief Abstract base class cosmic ray source features
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
 @brief Abstract base class for cosmic ray sources
 */
class SourceInterface : public Referenced {
public:
	virtual ref_ptr<Candidate> getCandidate() const = 0;
	virtual std::string getDescription() const = 0;
};

/**
 @class Source
 @brief General cosmic ray source

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
 @brief List of cosmic ray sources of individual lumosities.

 The SourceList is a source itself. It can be used if several sources are
 needed in one simulation.
 */
class SourceList: public SourceInterface {
	std::vector<ref_ptr<Source> > sources;
	std::vector<double> cdf;
public:
	void add(Source* source, double weight = 1);
	ref_ptr<Candidate> getCandidate() const;
	std::string getDescription() const;
};




/**
 @class SourceParticleType
 @brief Particle type at the source
 */
class SourceParticleType: public SourceFeature {
	int id;
public:
	SourceParticleType(int id);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
 @class SourceMultipleParticleTypes
 @brief Multiple particle types with individual relative abundances
 */
class SourceMultipleParticleTypes: public SourceFeature {
	std::vector<int> particleTypes;
	std::vector<double> cdf;
public:
	SourceMultipleParticleTypes();
	void add(int id, double weight = 1);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
 @class SourceEnergy
 @brief Sets the initial energy to a given value
 */
class SourceEnergy: public SourceFeature {
	double E;
public:
	SourceEnergy(double energy);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
 @class SourcePowerLawSpectrum
 @brief Particle energy following a power law spectrum
 */
class SourcePowerLawSpectrum: public SourceFeature {
	double Emin;
	double Emax;
	double index;
public:
	SourcePowerLawSpectrum(double Emin, double Emax, double index);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
 @class SourceComposition
 @brief Multiple nuclei with power law spectrum between Emin and Z * Rmax

 See Allard et al. 2006, DOI 10.1088/1475-7516/2006/09/005
 */
class SourceComposition: public SourceFeature {
	double Emin;
	double Rmax;
	double index;
	std::vector<int> nuclei;
	std::vector<double> cdf;
public:
	SourceComposition(double Emin, double Rmax, double index);
	void add(int id, double abundance);
	void add(int A, int Z, double abundance);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
 @class SourcePosition
 @brief Position of a point source
 */
class SourcePosition: public SourceFeature {
	Vector3d position; /**< Source position */
public:
	SourcePosition(Vector3d position);
	SourcePosition(double d);
	void prepareParticle(ParticleState &state) const;
	void setDescription();
};

/**
 @class SourceMultiplePositions
 @brief Multiple point source positions with individual luminosities
 */
class SourceMultiplePositions: public SourceFeature {
	std::vector<Vector3d> positions;
	std::vector<double> cdf;
public:
	SourceMultiplePositions();
	void add(Vector3d position, double weight = 1);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
 @class SourceUniformSphere
 @brief Uniform random source positions inside a sphere
 */
class SourceUniformSphere: public SourceFeature {
	Vector3d center;
	double radius;
public:
	SourceUniformSphere(Vector3d center, double radius);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
 @class SourceUniformHollowSphere
 @brief Uniform random source positions inside of a hollow sphere wall
 */
class SourceUniformHollowSphere: public SourceFeature {
	Vector3d center;
	double radius_inner;
	double radius_outer;
public:
	SourceUniformHollowSphere(Vector3d center,
			double radius_inner, double double_outer);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
 @class SourceUniformShell
 @brief Uniform random source positions on a sphere
 */
class SourceUniformShell: public SourceFeature {
	Vector3d center;
	double radius;
public:
	SourceUniformShell(Vector3d center, double radius);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
 @class SourceUniformBox
 @brief Uniform random source positions inside a box
 */
class SourceUniformBox: public SourceFeature {
	Vector3d origin;
	Vector3d size;
public:
	/** Constructor
	 @param origin	lower box corner
	 @param size	upper box corner
	 */
	SourceUniformBox(Vector3d origin, Vector3d size);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
 @class SourceUniformCylinder
 @brief Uniform random source positions inside a Cylinder
 */

class SourceUniformCylinder: public SourceFeature {
	Vector3d origin;
	double height;
	double radius;
public:
	/** Constructor
	 @param origin	lower middle of cylinder
	 @param height	height of the cylinder
	 @param radius	radius of the cylinder
*/
	SourceUniformCylinder(Vector3d origin, double height, double radius);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
@class SourceSNRDistribution
@brief Source distribution that follows the Galactic SNR distribution

The origin of the distribution is the Galactic center. The default maximum radius is set 
to R_max=20 kpc and the default maximum height is Z_max = 5 kpc.
See G. Case and D. Bhattacharya (1996) for the details of the distribution.
*/

class SourceSNRDistribution: public SourceFeature {
	double R_earth; // parameter given by observation
	double beta; // parameter to shift the maximum in R direction
	double Zg; // exponential cut parameter in z direction
	double frMax; // helper for efficient sampling
	double fzMax; // helper for efficient sampling
	double R_max; // maximum radial distance - default 20 kpc 
		      // (due to the extension of the JF12 field)
	double Z_max; // maximum distance from galactic plane - default 5 kpc

public:
	SourceSNRDistribution();	
	SourceSNRDistribution(double R_earth, double beta, double Zg);
	void prepareParticle(ParticleState &particle) const;
	double f_r(double r) const;
	double f_z(double z) const;
	void set_frMax(double R, double b);
	void set_fzMax(double Zg);
	void set_RMax(double R_max);
	void set_ZMax(double Z_max);
	double get_frMax();
	double get_fzMax();
	double get_RMax();
	double get_ZMax();
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
	double R_earth; // parameter given by observation
	double beta; // parameter to shift the maximum in R direction
	double Zg; // exponential cut parameter in z direction
	double frMax; // helper for efficient sampling
	double fzMax; // helper for efficient sampling
	double R_max;// maximum radial distance - default 22 kpc 
	double Z_max; // maximum distance from galactic plane - default 5 kpc 
	double r_blur; // relative smearing factor for the radius
	double theta_blur; // smearing factor for the angle. Unit = [1/length]


	

public:
	SourcePulsarDistribution();	
	SourcePulsarDistribution(double R_earth, double beta, double Zg, double r_blur, double theta_blur);
	void prepareParticle(ParticleState &particle) const;
	double f_r(double r) const;
	double f_z(double z) const;
	double f_theta(int i, double r) const;
	double blur_r(double r_tilde) const;
	double blur_theta(double theta_tilde, double r_tilde) const;
	void set_frMax(double R, double b);
	void set_fzMax(double Zg);
	void set_RMax(double R_max);
	void set_ZMax(double Z_max);
	void set_rBlur(double r_blur);
	void set_thetaBlur(double theta_blur);
	double get_frMax();
	double get_fzMax();
	double get_RMax();
	double get_ZMax();
	double get_rBlur();
	double get_thetaBlur();
	void setDescription();
};

/**
 @class SourceUniform1D
 @brief 1D-Positions from a uniform source distribution in an expanding universe

 This source property sets random x-coordinates according to a uniform source
 distribution in a given comoving distance interval.
 This is done by drawing a light travel distance from a flat distribution and
 converting to a comoving distance.
 */
class SourceUniform1D: public SourceFeature {
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
	void prepareParticle(ParticleState& particle) const;
	void setDescription();
};

/**
 @class SourceDensityGrid
 @brief Random source positions from a density grid
 */
class SourceDensityGrid: public SourceFeature {
	ref_ptr<ScalarGrid> grid;
public:
	SourceDensityGrid(ref_ptr<ScalarGrid> densityGrid);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
 @class SourceDensityGrid1D
 @brief Random source positions from a 1D density grid
 */
class SourceDensityGrid1D: public SourceFeature {
	ref_ptr<ScalarGrid> grid;
public:
	SourceDensityGrid1D(ref_ptr<ScalarGrid> densityGrid);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
 @class SourceIsotropicEmission
 @brief Isotropic emission from a source
 */
class SourceIsotropicEmission: public SourceFeature {
public:
	SourceIsotropicEmission();
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
 @class SourceDirection
 @brief Emission in a discrete direction
 */
class SourceDirection: public SourceFeature {
	Vector3d direction;
public:
	SourceDirection(Vector3d direction = Vector3d(-1, 0, 0));
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
 @class SourceEmissionMap
 @brief Deactivate Candidate if it has zero probability in provided EmissionMap
 */
class SourceEmissionMap: public SourceFeature {
	ref_ptr<EmissionMap> emissionMap;
public:
	SourceEmissionMap(EmissionMap *emissionMap);
	void prepareCandidate(Candidate &candidate) const;
	void setEmissionMap(EmissionMap *emissionMap);
	void setDescription();
};

/**
 @class SourceEmissionCone
 @brief Uniform random emission inside a cone
 */
class SourceEmissionCone: public SourceFeature {
	Vector3d direction;
	double aperture;
public:
	SourceEmissionCone(Vector3d direction, double aperture);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
 @class SourceRedshift
 @brief Discrete redshift (time of emission)
 */
class SourceRedshift: public SourceFeature {
	double z;
public:
	SourceRedshift(double z);
	void prepareCandidate(Candidate &candidate) const;
	void setDescription();
};

/**
 @class SourceUniformRedshift
 @brief Random redshift (time of emission) from uniform distribution
 */
class SourceUniformRedshift: public SourceFeature {
	double zmin, zmax;
public:
	SourceUniformRedshift(double zmin, double zmax);
	void prepareCandidate(Candidate &candidate) const;
	void setDescription();
};

/**
 @class SourceRedshiftEvolution
 @brief Random redshift (time of emission) from (1+z)^m distribution
 */
class SourceRedshiftEvolution: public SourceFeature {
	double zmin, zmax, m;
public:
	SourceRedshiftEvolution(double m, double zmin, double zmax);
	void prepareCandidate(Candidate &candidate) const;
};

/**
 @class SourceRedshift1D
 @brief Redshift according to the distance to 0

 This source property sets the redshift according to the distance to 0.
 It must be added after a position setting source property.
 */
class SourceRedshift1D: public SourceFeature {
public:
	SourceRedshift1D();
	void prepareCandidate(Candidate &candidate) const;
	void setDescription();
};

#ifdef CRPROPA_HAVE_MUPARSER
/**
 @class SourceGenericComposition
 @brief Multiple nuclei with energies described by an expression string
 */
class SourceGenericComposition: public SourceFeature {
public:
	struct Nucleus {
		int id;
		std::vector<double> cdf;
	};

	SourceGenericComposition(double Emin, double Emax, std::string expression, size_t bins = 1024);
	void add(int id, double abundance);
	void add(int A, int Z, double abundance);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();

    const std::vector<double> *getNucleusCDF(int id) const {
        for (size_t i = 0; i<nuclei.size(); i++)
            if (nuclei[i].id == id)
            	return &nuclei[i].cdf;
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

/**  @} */ // end of group SourceFeature
#endif

}// namespace crpropa

#endif // CRPROPA_SOURCE_H
