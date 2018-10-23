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
#include "crpropa/PhotonPropagation.h"
#include "crpropa/Random.h"
#include "crpropa/Referenced.h"
#include "crpropa/Source.h"
#include "crpropa/Units.h"
#include "crpropa/Variant.h"
#include "crpropa/Vector3.h"
#include "crpropa/Version.h"

#include "crpropa/module/Boundary.h"
#include "crpropa/module/BreakCondition.h"
#include "crpropa/module/DiffusionSDE.h"
#include "crpropa/module/EMCascade.h"
#include "crpropa/module/EMDoublePairProduction.h"
#include "crpropa/module/EMInverseComptonScattering.h"
#include "crpropa/module/EMPairProduction.h"
#include "crpropa/module/EMTripletPairProduction.h"
#include "crpropa/module/ElasticScattering.h"
#include "crpropa/module/ElectronPairProduction.h"
#include "crpropa/module/HDF5Output.h"
#include "crpropa/module/NuclearDecay.h"
#include "crpropa/module/Observer.h"
#include "crpropa/module/OutputShell.h"
#include "crpropa/module/ParticleCollector.h"
#include "crpropa/module/PhotoDisintegration.h"
#include "crpropa/module/PhotoPionProduction.h"
#include "crpropa/module/PhotonEleCa.h"
#include "crpropa/module/PhotonOutput1D.h"
#include "crpropa/module/PropagationCK.h"
#include "crpropa/module/Redshift.h"
#include "crpropa/module/RestrictToRegion.h"
#include "crpropa/module/SimplePropagation.h"
#include "crpropa/module/SynchrotronRadiation.h"
#include "crpropa/module/TextOutput.h"
#include "crpropa/module/Tools.h"
#include "crpropa/module/AdiabaticCooling.h"

#include "crpropa/magneticField/AMRMagneticField.h"
#include "crpropa/magneticField/JF12Field.h"
#include "crpropa/magneticField/JF12FieldSolenoidal.h"
#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/magneticField/MagneticFieldGrid.h"
#include "crpropa/magneticField/PT11Field.h"
#include "crpropa/magneticField/QuimbyMagneticField.h"
#include "crpropa/magneticField/ArchimedeanSpiralField.h"

#include "crpropa/advectionField/AdvectionField.h"

#include "crpropa/Massdistribution/Density.h"
#include "crpropa/Massdistribution/Nakanishi.h"
#include "crpropa/Massdistribution/Cordes.h"
#include "crpropa/Massdistribution/Massdistribution.h"
#include "crpropa/Massdistribution/Ferriere.h"
#include "crpropa/Massdistribution/Pohl.h"
#include "crpropa/Massdistribution/ConstantDensity.h"



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
 * @{ @brief Source and SourceFeatures
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
 * \defgroup Tools Tools
 * @{ @brief Collection of helper functinos and modules.
 * @}
 */






#endif // CRPROPA_H
