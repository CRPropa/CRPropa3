%module mpc

%include stl.i
%include std_set.i
%include std_multiset.i
%include std_map.i
%include std_pair.i
%include std_multimap.i
%include std_vector.i
%include std_string.i
%include stdint.i
%include std_container.i

%{
#include "mpc/module/NuclearDecay.h"
#include "mpc/module/ElectronPairProduction.h"
#include "mpc/module/PhotoDisintegration.h"
#include "mpc/module/PhotoPionProduction.h"
#include "mpc/module/Redshift.h"
#include "mpc/module/BreakCondition.h"
#include "mpc/module/Output.h"
#include "mpc/module/SimplePropagation.h"
#include "mpc/module/DeflectionCK.h"
#include "mpc/module/GlutDisplay.h"

#include "mpc/magneticField/magneticField.h"
#include "mpc/magneticField/uniformMagneticField.h"
#include "mpc/magneticField/magneticFieldGrid.h"
#include "mpc/magneticField/turbulentMagneticFieldGrid.h"
#include "mpc/magneticField/sphMagneticField.h"

#include "mpc/Candidate.h"
#include "mpc/ParticleState.h"
#include "mpc/Module.h"
#include "mpc/ModuleChain.h"
#include "mpc/PhasePoint.h"
#include "mpc/ExplicitRungeKutta.h"
#include "mpc/HepPID.h"
#include "mpc/Random.h"
#include "mpc/Units.h"
#include "mpc/Vector3.h"
#include "mpc/Referenced.h"
#include "mpc/module/common.h"
%}
 
/* Parse the header file to generate wrappers */
%feature("ref")   mpc::Referenced "$this->addReference();"
%feature("unref") mpc::Referenced "$this->removeReference();"
%include "mpc/Referenced.h"
%include "mpc/Units.h"
%include "mpc/HepPID.h"
%include "mpc/Vector3.h"
%include "mpc/Random.h"
%include "mpc/ParticleState.h"
%include "mpc/module/common.h"
%template(CandidateVector) std::vector< mpc::ref_ptr<mpc::Candidate> >;
%template(CandidateRefPtr) mpc::ref_ptr<mpc::Candidate>;
%include "mpc/Candidate.h"

%template(ModuleVector) std::vector< mpc::ref_ptr<mpc::Module> >;
%template(ModuleRefPtr) mpc::ref_ptr<mpc::Module>;
%include "mpc/Module.h"

%template(MagneticFieldVector) std::vector< mpc::ref_ptr<mpc::MagneticField> >;
%template(MagneticFieldRefPtr) mpc::ref_ptr<mpc::MagneticField>;
%include "mpc/magneticField/magneticField.h"
%include "mpc/magneticField/magneticFieldGrid.h"
%include "mpc/magneticField/uniformMagneticField.h"
%include "mpc/magneticField/sphMagneticField.h"
%include "mpc/magneticField/turbulentMagneticFieldGrid.h"

%include "mpc/ExplicitRungeKutta.h"
%include "mpc/PhasePoint.h"

%include "mpc/module/BreakCondition.h"
%include "mpc/module/SimplePropagation.h"
%include "mpc/module/DeflectionCK.h"
%include "mpc/module/Output.h"
%include "mpc/module/NuclearDecay.h"
%include "mpc/module/ElectronPairProduction.h"
%include "mpc/module/PhotoPionProduction.h"
%include "mpc/module/PhotoDisintegration.h"
%include "mpc/module/Redshift.h"
%include "mpc/module/GlutDisplay.h"


%include "mpc/ModuleChain.h"
