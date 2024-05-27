/// CRPRopa public definitions
#ifndef CRPROPA_H
#define CRPROPA_H

#include "crpropa/Candidate.h"
#include "crpropa/Common.h"
#include "crpropa/Cosmology.h"
#include "crpropa/EmissionMap.h"
#include "crpropa/Geometry.h"
#include "crpropa/Grid.h"
#include "crpropa/GridTools.h"
#include "crpropa/Logging.h"
#include "crpropa/Module.h"
#include "crpropa/ModuleList.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/ParticleState.h"
#include "crpropa/PhotonBackground.h"
#include "crpropa/Random.h"
#include "crpropa/Referenced.h"
#include "crpropa/Source.h"
#include "crpropa/Units.h"
#include "crpropa/Variant.h"
#include "crpropa/Vector3.h"
#include "crpropa/Version.h"

#include "crpropa/module/AdiabaticCooling.h"
#include "crpropa/module/Acceleration.h"
#include "crpropa/module/Boundary.h"
#include "crpropa/module/BreakCondition.h"
#include "crpropa/module/CandidateSplitting.h"
#include "crpropa/module/DiffusionSDE.h"
#include "crpropa/module/EMDoublePairProduction.h"
#include "crpropa/module/EMInverseComptonScattering.h"
#include "crpropa/module/EMPairProduction.h"
#include "crpropa/module/EMTripletPairProduction.h"
#include "crpropa/module/ElasticScattering.h"
#include "crpropa/module/ElectronPairProduction.h"
#include "crpropa/module/HDF5Output.h"
#include "crpropa/module/MomentumDiffusion.h"
#include "crpropa/module/NuclearDecay.h"
#include "crpropa/module/Observer.h"
#include "crpropa/module/OutputShell.h"
#include "crpropa/module/ParticleCollector.h"
#include "crpropa/module/PhotoDisintegration.h"
#include "crpropa/module/PhotoPionProduction.h"
#include "crpropa/module/PhotonOutput1D.h"
#include "crpropa/module/PropagationBP.h"
#include "crpropa/module/PropagationCK.h"
#include "crpropa/module/Redshift.h"
#include "crpropa/module/RestrictToRegion.h"
#include "crpropa/module/SimplePropagation.h"
#include "crpropa/module/SynchrotronRadiation.h"
#include "crpropa/module/TextOutput.h"
#include "crpropa/module/Tools.h"

#include "crpropa/magneticField/ArchimedeanSpiralField.h"
#include "crpropa/magneticField/JF12Field.h"
#include "crpropa/magneticField/JF12FieldSolenoidal.h"
#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/magneticField/GalacticMagneticField.h"
#include "crpropa/magneticField/MagneticFieldGrid.h"
#include "crpropa/magneticField/PolarizedSingleModeMagneticField.h"
#include "crpropa/magneticField/PT11Field.h"
#include "crpropa/magneticField/QuimbyMagneticField.h"
#include "crpropa/magneticField/TF17Field.h"
#include "crpropa/magneticField/UF23Field.h"
#include "crpropa/magneticField/CMZField.h"
#include "crpropa/magneticField/turbulentField/GridTurbulence.h"
#include "crpropa/magneticField/turbulentField/HelicalGridTurbulence.h"
#include "crpropa/magneticField/turbulentField/PlaneWaveTurbulence.h"
#include "crpropa/magneticField/turbulentField/SimpleGridTurbulence.h"
#include "crpropa/magneticField/turbulentField/TurbulentField.h"

#include "crpropa/advectionField/AdvectionField.h"

#include "crpropa/massDistribution/Density.h"
#include "crpropa/massDistribution/Nakanishi.h"
#include "crpropa/massDistribution/Cordes.h"
#include "crpropa/massDistribution/Massdistribution.h"
#include "crpropa/massDistribution/Ferriere.h"
#include "crpropa/massDistribution/ConstantDensity.h"

/** \namespace crpropa
 *  @brief CRPropa is a public astrophysical simulation framework for propagating extraterrestrial ultra-high energy particles.
 **/

// Groups of Modules for Doxygen
/**
 * \defgroup Core Core Classes
 * @{ @brief Core classes used to build CRPropa
 * @}
 *
 * \defgroup PhysicsDefinitions Physics Definitions
 * @{ @brief Fundamental physical data, function, units, constants, etc.
 * @}
 *
 * \defgroup Propagation Propagation Modules
 * @{ @brief Modules that propagate a Candidate
 * @}
 *
 * \defgroup EnergyLosses Energy Loss Processes
 * @{ @brief Energy losses of candidates.
 * @}
 *
 * \defgroup PhotonFields Photon Fields
 * @{ @brief Photon fields for particle interactions.
 * @}
 *
 * \defgroup MagneticFields Magnetic Fields
 * @{ @brief Magnetic field models
 * @}
 *
 * \defgroup Observer Observers
 * @{ @brief Observers and ObserverFeatures
 *
 * ObserverFeatures are added to the ObserverModule to check for detection and perform actions on detection.
 * @}
 *
 * \defgroup Condition Conditions
 * @{ @brief Propagation boundaries and breaking conditions
 *
 * Conditions are set of modules similar to the Observer as they can perform an action if fulfilled.
 * @}
 *
 * \defgroup SourceFeatures Sources
 * @{ @brief Source and SourceFeatures, i.e. properties of the partcle injection
 *
 *  Sourcefeatures are added to sources and manipulate the properties of the
 *  emitted candidate.
 * @}
 *
 * \defgroup Output Output Modules
 * @{ @brief File outputs.
 * @}
 *
 * \defgroup MagneticLenses Magnetic Lenses
 * @{ @brief Lensing technique to account for deflections in the galactic
 * magnetic field.
 * @}
 *
 * \defgroup Acceleration Acceleration Processes
 * @{ @brief Modules and techniques to simualte particle acceleration
 *
 *  These particle acceleration features is analogue to an energy loss due to
 *  an interaction with a background photon, but for energy loss/gain in
 *  interaction with scattering centers.
 * @}
 *
 * \defgroup Tools Tools
 * @{ @brief Collection of helper functions and modules.
 * @}
 */






#endif // CRPROPA_H
