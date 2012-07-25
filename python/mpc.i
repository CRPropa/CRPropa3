%module(directors="1") mpc

%include stl.i
%include std_set.i
%include std_multiset.i
%include std_map.i
%include std_pair.i
%include std_multimap.i
%include std_vector.i
%include std_string.i
%include std_list.i
%include stdint.i
%include std_container.i
%include "exception.i"

%{
#include "mpc/module/NuclearDecay.h"
#include "mpc/module/ElectronPairProduction.h"
#include "mpc/module/StochasticInteraction.h"
#include "mpc/module/PhotoDisintegration.h"
#include "mpc/module/PhotoPionProduction.h"
#include "mpc/module/Redshift.h"
#include "mpc/module/BreakCondition.h"
#include "mpc/module/Output.h"
#include "mpc/module/SimplePropagation.h"
#include "mpc/module/DeflectionCK.h"
#include "mpc/module/Tools.h"

#include "mpc/magneticField/MagneticField.h"
#include "mpc/magneticField/UniformMagneticField.h"
#include "mpc/magneticField/PeriodicGrid.h"
#include "mpc/magneticField/MagneticFieldGrid.h"
#include "mpc/magneticField/SPHMagneticField.h"

#include "mpc/Referenced.h"
#include "mpc/Candidate.h"
#include "mpc/ParticleState.h"
#include "mpc/Nucleus.h"
#include "mpc/Module.h"
#include "mpc/ModuleList.h"
#include "mpc/PhasePoint.h"
#include "mpc/ExplicitRungeKutta.h"
#include "mpc/Random.h"
#include "mpc/Units.h"
#include "mpc/Vector3.h"
#include "mpc/Source.h"
#include "mpc/Common.h"
%}


%exception
{
 try
 {
   $action
 }
 catch (const std::runtime_error& e) {
   SWIG_exception(SWIG_RuntimeError, e.what());
 }
 catch (...) { 
   SWIG_exception(SWIG_RuntimeError, "unknown exception");
 } 
}


%ignore operator<<;
%ignore operator>>;

%ignore operator mpc::Source*;
%ignore operator mpc::Candidate*;
%ignore operator mpc::Module*;
%ignore operator mpc::ModuleList*;
%ignore operator mpc::MagneticField*;

%feature("ref")   mpc::Referenced "$this->addReference();"
%feature("unref") mpc::Referenced "$this->removeReference();"


%include "mpc/Vector3.h"
%template(Vector3d) mpc::Vector3<double>;
%template(Vector3f) mpc::Vector3<float>;


%include "mpc/Referenced.h"
%include "mpc/Units.h"
%include "mpc/Nucleus.h"
%include "mpc/Common.h"
%include "mpc/Random.h"
%include "mpc/ParticleState.h"


%template(SourceRefPtr) mpc::ref_ptr<mpc::Source>;
%include "mpc/Source.h"

%template(CandidateVector) std::vector< mpc::ref_ptr<mpc::Candidate> >;
%template(CandidateRefPtr) mpc::ref_ptr<mpc::Candidate>;
%include "mpc/Candidate.h"

%template(ModuleRefPtr) mpc::ref_ptr<mpc::Module>;
%template(stdModuleList) std::list< mpc::ref_ptr<mpc::Module> >;
%include "mpc/Module.h"


%implicitconv mpc::ref_ptr<mpc::MagneticField>;
%template(MagneticFieldRefPtr) mpc::ref_ptr<mpc::MagneticField>;
%include "mpc/magneticField/MagneticField.h"
%include "mpc/magneticField/UniformMagneticField.h"

%include "mpc/magneticField/PeriodicGrid.h"

%implicitconv mpc::ref_ptr<mpc::PeriodicGrid<mpc::Vector3<float> > >;
%template(VectorGridRefPtr) mpc::ref_ptr<mpc::PeriodicGrid<mpc::Vector3<float> > >;
%template(VectorGrid) mpc::PeriodicGrid<mpc::Vector3<float> >;

%implicitconv mpc::ref_ptr<mpc::PeriodicGrid<float> >;
%template(ScalarGridRefPtr) mpc::ref_ptr<mpc::PeriodicGrid<float> >;
%template(ScalarGrid) mpc::PeriodicGrid<float>;

%include "mpc/magneticField/MagneticFieldGrid.h"
%include "mpc/magneticField/SPHMagneticField.h"

%include "mpc/ExplicitRungeKutta.h"
%include "mpc/PhasePoint.h"


%include "mpc/module/BreakCondition.h"
%include "mpc/module/SimplePropagation.h"
%include "mpc/module/DeflectionCK.h"
%include "mpc/module/Output.h"
%include "mpc/module/ElectronPairProduction.h"
%include "mpc/module/StochasticInteraction.h"
%include "mpc/module/NuclearDecay.h"
%include "mpc/module/PhotoPionProduction.h"
%include "mpc/module/PhotoDisintegration.h"
%include "mpc/module/Redshift.h"
%include "mpc/module/Tools.h"


%template(ModuleListRefPtr) mpc::ref_ptr<mpc::ModuleList>;
%include "mpc/ModuleList.h"