%module(directors="1") crpropa
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

#ifdef CRPROPA_HAVE_QUIMBY
%import (module="quimby") quimby.i
#endif

%{
#include "crpropa/module/NuclearDecay.h"
#include "crpropa/module/ElectronPairProduction.h"
#include "crpropa/module/PhotoDisintegration.h"
#include "crpropa/module/PhotoPionProduction.h"
#include "crpropa/module/Redshift.h"
#include "crpropa/module/BreakCondition.h"
#include "crpropa/module/Boundary.h"
#include "crpropa/module/Observer.h"
#include "crpropa/module/OutputTXT.h"
#include "crpropa/module/OutputShell.h"
#include "crpropa/module/OutputROOT.h"
#include "crpropa/module/OutputCRPropa2.h"
#include "crpropa/module/PhotonDINT.h"
#include "crpropa/module/PhotonEleCa.h"
#include "crpropa/module/SimplePropagation.h"
#include "crpropa/module/DeflectionCK.h"
#include "crpropa/module/Tools.h"

#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/magneticField/MagneticFieldGrid.h"
#include "crpropa/magneticField/QuimbyMagneticField.h"
#include "crpropa/magneticField/JF12Field.h"

#include "crpropa/Referenced.h"
#include "crpropa/Candidate.h"
#include "crpropa/ParticleState.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Module.h"
#include "crpropa/ModuleList.h"
#include "crpropa/PhasePoint.h"
#include "crpropa/ExplicitRungeKutta.h"
#include "crpropa/Random.h"
#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Source.h"
#include "crpropa/Common.h"
#include "crpropa/Cosmology.h"
#include "crpropa/PhotonBackground.h"
#include "crpropa/Grid.h"
#include "crpropa/GridTools.h"
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
%ignore operator crpropa::Source*;
%ignore operator crpropa::Candidate*;
%ignore operator crpropa::Module*;
%ignore operator crpropa::ModuleList*;
%ignore operator crpropa::MagneticField*;

%feature("ref")   crpropa::Referenced "$this->addReference();"
%feature("unref") crpropa::Referenced "$this->removeReference();"


%include "crpropa/Vector3.h"
%template(Vector3d) crpropa::Vector3<double>;
%template(Vector3f) crpropa::Vector3<float>;

%include "crpropa/Referenced.h"
%include "crpropa/Units.h"
%include "crpropa/Common.h"
%include "crpropa/Cosmology.h"
%include "crpropa/PhotonBackground.h"
%include "crpropa/Random.h"
%include "crpropa/ParticleState.h"
%include "crpropa/ParticleID.h"
%include "crpropa/ParticleMass.h"

%template(CandidateVector) std::vector< crpropa::ref_ptr<crpropa::Candidate> >;
%template(CandidateRefPtr) crpropa::ref_ptr<crpropa::Candidate>;
%include "crpropa/Candidate.h"

%template(ModuleRefPtr) crpropa::ref_ptr<crpropa::Module>;
%template(stdModuleList) std::list< crpropa::ref_ptr<crpropa::Module> >;
%feature("director") crpropa::Module;      
%include "crpropa/Module.h"

%implicitconv crpropa::ref_ptr<crpropa::MagneticField>;
%template(MagneticFieldRefPtr) crpropa::ref_ptr<crpropa::MagneticField>;
%include "crpropa/magneticField/MagneticField.h"

%include "crpropa/Grid.h"
%include "crpropa/GridTools.h"

%implicitconv crpropa::ref_ptr<crpropa::Grid<crpropa::Vector3<float> > >;
%template(VectorGridRefPtr) crpropa::ref_ptr<crpropa::Grid<crpropa::Vector3<float> > >;
%template(VectorGrid) crpropa::Grid<crpropa::Vector3<float> >;

%implicitconv crpropa::ref_ptr<crpropa::Grid<float> >;
%template(ScalarGridRefPtr) crpropa::ref_ptr<crpropa::Grid<float> >;
%template(ScalarGrid) crpropa::Grid<float>;

%include "crpropa/magneticField/MagneticFieldGrid.h"
%include "crpropa/magneticField/QuimbyMagneticField.h"
%include "crpropa/magneticField/JF12Field.h"

%include "crpropa/ExplicitRungeKutta.h"
%include "crpropa/PhasePoint.h"

%include "crpropa/module/BreakCondition.h"
%include "crpropa/module/Boundary.h"
%include "crpropa/module/Observer.h"
%include "crpropa/module/SimplePropagation.h"
%include "crpropa/module/DeflectionCK.h"
%include "crpropa/module/OutputTXT.h"
%include "crpropa/module/OutputShell.h"
%include "crpropa/module/OutputROOT.h"
%include "crpropa/module/OutputCRPropa2.h"
%include "crpropa/module/PhotonDINT.h"
%include "crpropa/module/PhotonEleCa.h"
%include "crpropa/module/ElectronPairProduction.h"
%include "crpropa/module/NuclearDecay.h"
%include "crpropa/module/PhotoPionProduction.h"
%include "crpropa/module/PhotoDisintegration.h"
%include "crpropa/module/Redshift.h"
%include "crpropa/module/Tools.h"

%template(SourceRefPtr) crpropa::ref_ptr<crpropa::Source>;
%include "crpropa/Source.h"

%template(ModuleListRefPtr) crpropa::ref_ptr<crpropa::ModuleList>;
%include "crpropa/ModuleList.h"


// nice python print
%pythoncode %{
Module.__str__ = Module.getDescription
def Vector3__str__(self):
  return "(%.4e, %.4e, %.4e)" % (self.x, self.y, self.z)
Vector3d.__str__ = Vector3__str__
Vector3f.__str__ = Vector3__str__
%}
