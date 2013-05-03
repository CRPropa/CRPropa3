%module(directors="1") mpc
%feature("autodoc", "1"); // automatic docstrings

%{
// workaround for SWIG < 2.0.5 with GCC >= 4.7
#include <cstddef>
using std::ptrdiff_t;
%}

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

#ifdef MPC_HAVE_QUIMBY
%import (module="quimby") quimby.i
#endif

%{
#include "mpc/module/NuclearDecay.h"
#include "mpc/module/ElectronPairProduction.h"
#include "mpc/module/StochasticInteraction.h"
#include "mpc/module/PhotoDisintegration.h"
#include "mpc/module/PhotoPionProduction.h"
#include "mpc/module/Redshift.h"
#include "mpc/module/BreakCondition.h"
#include "mpc/module/Boundary.h"
#include "mpc/module/Observer.h"
#include "mpc/module/Output.h"
#include "mpc/module/OutputROOT.h"
#include "mpc/module/OutputCRPropa2.h"
#include "mpc/module/PhotonOutput.h"
#include "mpc/module/SimplePropagation.h"
#include "mpc/module/DeflectionCK.h"
#include "mpc/module/Tools.h"

#include "mpc/magneticField/MagneticField.h"
#include "mpc/magneticField/MagneticFieldGrid.h"
#include "mpc/magneticField/QuimbyMagneticField.h"
#include "mpc/magneticField/JF12Field.h"

#include "mpc/Referenced.h"
#include "mpc/Candidate.h"
#include "mpc/ParticleState.h"
#include "mpc/ParticleID.h"
#include "mpc/ParticleMass.h"
#include "mpc/Module.h"
#include "mpc/ModuleList.h"
#include "mpc/PhasePoint.h"
#include "mpc/ExplicitRungeKutta.h"
#include "mpc/Random.h"
#include "mpc/Units.h"
#include "mpc/Vector3.h"
#include "mpc/Source.h"
#include "mpc/Common.h"
#include "mpc/Cosmology.h"
#include "mpc/PhotonBackground.h"
#include "mpc/Grid.h"
#include "mpc/GridTools.h"
%}


%exception
{
 try
 {
   $action
 }
 catch (const std::exception& e) {
   SWIG_exception(SWIG_RuntimeError, e.what());
 }
 catch (const char *e) {
   SWIG_exception(SWIG_RuntimeError, e);
 }
 catch (Swig::DirectorException &e) {
   SWIG_exception(SWIG_RuntimeError, e.getMessage());
 }
}

%ignore operator<<;
%ignore operator>>;
%ignore *::operator=;
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
%include "mpc/Common.h"
%include "mpc/Cosmology.h"
%include "mpc/PhotonBackground.h"
%include "mpc/Random.h"
%include "mpc/ParticleState.h"
%include "mpc/ParticleID.h"
%include "mpc/ParticleMass.h"

%template(CandidateVector) std::vector< mpc::ref_ptr<mpc::Candidate> >;
%template(CandidateRefPtr) mpc::ref_ptr<mpc::Candidate>;
%include "mpc/Candidate.h"

%template(ModuleRefPtr) mpc::ref_ptr<mpc::Module>;
%template(stdModuleList) std::list< mpc::ref_ptr<mpc::Module> >;
%feature("director") mpc::Module;      
%include "mpc/Module.h"

%implicitconv mpc::ref_ptr<mpc::MagneticField>;
%template(MagneticFieldRefPtr) mpc::ref_ptr<mpc::MagneticField>;
%include "mpc/magneticField/MagneticField.h"

%include "mpc/Grid.h"
%include "mpc/GridTools.h"

%implicitconv mpc::ref_ptr<mpc::Grid<mpc::Vector3<float> > >;
%template(VectorGridRefPtr) mpc::ref_ptr<mpc::Grid<mpc::Vector3<float> > >;
%template(VectorGrid) mpc::Grid<mpc::Vector3<float> >;

%implicitconv mpc::ref_ptr<mpc::Grid<float> >;
%template(ScalarGridRefPtr) mpc::ref_ptr<mpc::Grid<float> >;
%template(ScalarGrid) mpc::Grid<float>;

%include "mpc/magneticField/MagneticFieldGrid.h"
%include "mpc/magneticField/QuimbyMagneticField.h"
%include "mpc/magneticField/JF12Field.h"

%include "mpc/ExplicitRungeKutta.h"
%include "mpc/PhasePoint.h"

%include "mpc/module/BreakCondition.h"
%include "mpc/module/Boundary.h"
%include "mpc/module/Observer.h"
%include "mpc/module/SimplePropagation.h"
%include "mpc/module/DeflectionCK.h"
%include "mpc/module/Output.h"
%include "mpc/module/OutputROOT.h"
%include "mpc/module/OutputCRPropa2.h"
%include "mpc/module/PhotonOutput.h"
%include "mpc/module/ElectronPairProduction.h"
%include "mpc/module/StochasticInteraction.h"
%include "mpc/module/NuclearDecay.h"
%include "mpc/module/PhotoPionProduction.h"
%include "mpc/module/PhotoDisintegration.h"
%include "mpc/module/Redshift.h"
%include "mpc/module/Tools.h"

%template(SourceRefPtr) mpc::ref_ptr<mpc::Source>;
%include "mpc/Source.h"

%template(ModuleListRefPtr) mpc::ref_ptr<mpc::ModuleList>;
%include "mpc/ModuleList.h"


// nice python print
%pythoncode %{
Module.__str__ = Module.getDescription
def Vector3__str__(self):
  return "(%.4e, %.4e, %.4e)" % (self.x, self.y, self.z)
Vector3d.__str__ = Vector3__str__
Vector3f.__str__ = Vector3__str__
%}
