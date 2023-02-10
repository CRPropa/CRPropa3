#ifndef CRPROPA_ACCELERATION
#define CRPROPA_ACCELERATION

#include <crpropa/Candidate.h>
#include <crpropa/Geometry.h>
#include <crpropa/Module.h>
#include <crpropa/Units.h>
#include <crpropa/Vector3.h>

#include <string>

namespace crpropa {
/** @addtogroup Acceleration
 *  @{
 */

/// @class StepLengthModifier
/// @brief Modifies the steplength of an acceleration module.
class StepLengthModifier : public Referenced {
  public:
	/// Returns an update of the steplength
	/// @param steplength 	Modifies step length, e.g., based on scattering 
	///						model.	
	/// @param candidate 	Additional candidate properties are usually 
	///						included in the calculation of the updated
	///						step length.
	virtual double modify(double steplength, Candidate *candidate) = 0;
};


/// @class AbstractAccelerationModule
/// @brief Core functionallity for acceleration by scattering with scatter
///  centers moving in a velocity field.
/// @details The velocity field is implicity implemented in the derived classes
///  for performance reasons. Models for the dependence of the step length of
///  the scatter process are set via modifiers.
class AbstractAccelerationModule : public Module {
	double stepLength;
	std::vector<ref_ptr<StepLengthModifier>> modifiers;

  public:
	/// The parent's constructor need to be called on initialization!
	AbstractAccelerationModule(double _stepLength = 1. * parsec);
	// add a step length modifier to the model
	void add(StepLengthModifier *modifier);
	// update the candidate
	void process(Candidate *candidate) const;

	/// Returns the velocity vector of the scatter centers in the rest frame of
	/// the candidate. Needs to be implemented in inheriting classes.
	virtual Vector3d scatterCenterVelocity(Candidate *candidate) const = 0;

	/// Scatter the candidate with a center with given scatter center
	/// velocity into a random direction. Assumes that the
	/// candidate is ultra-relativistic (m = 0).
	void scatter(Candidate *candidate,
	             const Vector3d &scatter_center_velocity) const;
};


/// @class SecondOrderFermi
/// @brief  Implements scattering with centers moving in isotropic directions.
///   All scatter centers have the same velocity.
class SecondOrderFermi : public AbstractAccelerationModule {
	double scatterVelocity;
	std::vector<double> angle;
	std::vector<double> angleCDF;

  public:
	/** Constructor
	@param scatterVelocity			velocity of scattering centers
	@param stepLength				average mean free path
	@param sizeOfPitchangleTable	number of precalculated pitch angles
	*/
	SecondOrderFermi(double scatterVelocity = .1 * crpropa::c_light,
	                 double stepLength = 1. * crpropa::parsec,
	                 unsigned int sizeOfPitchangleTable = 10000);
	virtual crpropa::Vector3d
	scatterCenterVelocity(crpropa::Candidate *candidate) const;
};


/// @class DirectedFlowScattering
/// @brief Scattering in a directed flow of scatter centers.
/// @details Two of these region with different
///    velocities can be used to create first order Fermi scenario.
///    Thanks to Aritra Ghosh, Groningn University, for first work in 2017 on
///    the shock acceleration in CRPropa leading to this module.
class DirectedFlowScattering : public AbstractAccelerationModule {
  private:
	crpropa::Vector3d __scatterVelocity;

  public:
  /** Constructor
   * @param scatterCenterVelocity	velocity of scattering centers
   * @param stepLength				average mean free path
  */
	DirectedFlowScattering(crpropa::Vector3d scatterCenterVelocity,
	                       double stepLength = 1. * parsec);
	virtual crpropa::Vector3d
	scatterCenterVelocity(crpropa::Candidate *candidate) const;
};


/// @class DirectedFlowOfScatterCenters
/// @brief In a directed flow, the step length depend on the direction of the
/// particles as headon collisions are more likely than tail=on collisions -
/// propagating against the flow is harder.
class DirectedFlowOfScatterCenters : public StepLengthModifier {
  private:
	Vector3d __scatterVelocity;

  public:
  /** Constructor
   * @param scatterCenterVelocity	velocity of scattering centers
  */
	DirectedFlowOfScatterCenters(const Vector3d &scatterCenterVelocity);
	double modify(double steplength, Candidate *candidate);
};


/// @class QuasiLinearTheory
/// @brief Scales the steplength according to quasi linear theory.
/// @details Following quasi-linear theory [Schlickeiser1989], the mean free
/// path \f$\lambda\f$ of a
///  particle with energy \f$E\f$ and charge \f$Z\f$ in a field with turbulence
///  spectrum \f$\frac{k}{k_{\min}}^{-q}\f$ is
///   \f$ \lambda = {\left(\frac{B}{\delta B}\right)}^2 {\left(R_G\;
///   k_{\min}\right)}^{1-q} R_G \equiv \lambda_0 {\left( \frac{E}{1
///   EeV}\frac{1}{Z} \right)}^{2-q} \f$
///  where \f$R_G = \frac{E}{B Z}\f$ is the gyro-radius of the
///  particles.
/// This class implements the rigidity dependent scaling factor used to modify
/// the base step length.
/// @par
/// \b [Schlickeiser1989] R. Schlickeiser, Cosmic-Ray Transport and
///      Acceleration. II. Cosmic Rays in Moving Cold Media with Application to
///      Diffusive Shock Wave Acceleration,
///      The Astrophysical Journal 336 (1989) 264. doi:10.1086/167010.
class QuasiLinearTheory : public StepLengthModifier {
	private:
	double __referenceEnergy;
	double __turbulenceIndex;
	double __minimumRigidity;

  public:
  /** Constructor
   * @param referenecEnergy	reference energy - break of power spectrum
   * @param turbulenceIndex	power law index of the isotropic magnetic 
   * 						turbulence power spectrum; default is set 
   * 						to Kolmogorov turbulence.
   * @param minimumRigidity	minimal rigidity
  */
	QuasiLinearTheory(double referenecEnergy = 1. * EeV,
	                  double turbulenceIndex = 5. / 3,
	                  double minimumRigidity = 0);
	double modify(double steplength, Candidate *candidate);
};


/// @class ParticleSplitting
/// @brief  Implements particle splitting, i.e. inverse thinning, to speed up
///  the simulation.
/// @details After crossing a surface a given number of times, the particle is
/// split to N partilces with weight 1/N. This eases performance constraints in
/// acceleration simulations due to the power law nature of many acceleration
/// mechanisms.
/// Thanks to Matthew Weiss, Penn State University for the first work on this
/// feature in 2017.
class ParticleSplitting : public Module {
	int numberSplits;
	int crossingThreshold;
	double minWeight;
	ref_ptr<Surface> surface;
	std::string counterid;

	public:
	/** Constructor
	@param surface              The surface to monitor
	@param crossingThreshold   Number of crossings after which a particle is split
	@param numberSplits           Number of particles the candidate is split into
	@param minWeight           Minimum weight to consider. Particles with
	                         	a lower weight are not split again.
	@param counterid            An unique string to identify the particle
	                            property used for counting. Useful if
	                            multiple splitting modules are present.
	*/
	ParticleSplitting(Surface *surface, int crossingThreshold = 50,
	                  int numberSplits = 5, double minWeight = 0.01,
	                  std::string counterid = "ParticleSplittingCounter");

	// update the candidate
	void process(Candidate *candidate) const;
};


/**  @} */ // end of group Acceleration

} // namespace crpropa

#endif
