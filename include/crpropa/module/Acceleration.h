#ifndef ACCELERATION_MODULE
#define ACCELERATION_MODULE

#include <crpropa/Candidate.h>
#include <crpropa/Module.h>
#include <crpropa/Vector3.h>
#include <crpropa/Units.h>

namespace crpropa
{

/// Modifies the steplength of an acceleration module.
class StepLengthModifier : public Referenced
{
	public:
		/// Returns an update of the steplength
		virtual double modify(double steplength, Candidate* candidate) = 0;
};


// Core functionallity for acceleration by scattering with scatter centers
// moving in a velocity field. The velocity field is implicity implemented in
// the derived class for performance reasons.
// Models for the dependence of the step length of
// the scatter process are set via modifiers.
class AbstractAccelerationModule : public Module
{
	double stepLength;
	std::vector<ref_ptr<StepLengthModifier> > modifiers;

	public:
		/// The parent's constructor need to be called on initialization!
		AbstractAccelerationModule(double _stepLength = 1. * parsec);
		// add a step length modifier to the model
		void add(StepLengthModifier *modifier);
		// update the candidate
		void process(Candidate *candidate) const;
		/// Returns the velocity vector of the scatter centers in the Lab Frame.
		/// Needs to be implemented in inheriting classes.
		virtual Vector3d scatterCenterVelocity(Candidate *candidate) const = 0;

		/// Returns new candidate momentum post collision when supplied with initial candidate momentum
		/// [All in the rest frame of the scatter center]. Can be overriden in inheriting classes.
		/// The default is a uniform distribution.
		Vector3d scatterMomentum(const Vector3d &v) const;

		/// Scatter the candidate with a center and the given scatter center
		/// velocity
		/// Scattering the candidate into a random direction. Assumes that the
		/// candidate is ultra-relativistic (m = 0).
		void scatter(Candidate* candidate, const Vector3d& scatter_center_velocity) const;
};






// Implements scattering with centers moving in isotropic directions. Scatter
// centers have the same velocity.
class SecondOrderFermi : public AbstractAccelerationModule
{
	double scatterVelocity;
	std::vector<double> angle;
	std::vector<double> angleCDF;

	public:
		SecondOrderFermi(double _scatterVelocity=.1 * crpropa::c_light, double stepLength = 1. * crpropa::parsec, unsigned int size_of_pitchangle_table=10000);
		virtual crpropa::Vector3d scatterCenterVelocity(crpropa::Candidate *candidate) const;
};


/// Step length depend on the direction of the particle - propagating against
/// the flow is harder.
class DirectedFlowOfScatterCenters: public StepLengthModifier
{
	private:
		Vector3d __scatterVelocity;
	public:
	DirectedFlowOfScatterCenters(const Vector3d &scatterCenterVelocity);
	double modify(double steplength, Candidate* candidate);
};


/*
		A directed flow of scatter centers. Two of these region with different
		velocities can be used to create a shock, i.e. first order Fermi.
		Thanks to Aritra Ghosh, Groningn University, for first work in 2017 on the
		shock acceleration in CRPropa leading to this module.
*/
class DirectedFlowScattering : public AbstractAccelerationModule
{
	private:
		crpropa::Vector3d __scatterVelocity;
	public:
		/// In a directed field of scatter centers, the probability to encounter a
		/// scatter center depends on the direction of the candidate.
		DirectedFlowScattering(crpropa::Vector3d scatterCenterVelocity, double stepLength=1. * parsec);
		virtual crpropa::Vector3d scatterCenterVelocity(crpropa::Candidate *candidate) const;
};

/*
Scales the steplength according to quasi linear theory via
		L = L0 * ( Z * (E_ref / E))**(2-q)
*/
class QuasiLinearTheory : public StepLengthModifier
{
	private:
		double __referenceEnergy;
		double __turbulenceIndex;
		double __minimumRigidity;

	public:
		QuasiLinearTheory(double referenecEnergy=1.*EeV, double turbulence_index=5./3, double minimumRigidity=0);
		double modify(double steplength, Candidate* candidate);
};



} // namespace crpropa

#endif
