/* 2: SWIG and CRPropa headers */

%include "stl.i"
%include "std_set.i"
%include "std_multiset.i"
%include "std_map.i"
%include "std_pair.i"
%include "std_multimap.i"
%include "std_vector.i"
%include "std_string.i"
%include "std_list.i"
%include "stdint.i"
%include "std_container.i"
%include "exception.i"

#ifdef CRPROPA_HAVE_QUIMBY
%import (module="quimby") "quimby/Referenced.h"
%import (module="quimby") "quimby/Vector3.h"
%import (module="quimby") "quimby/MagneticField.h"
//%import (module="quimby") quimby.i
#endif

#ifdef CRPROPA_HAVE_SAGA
%import (module="saga") saga.i
#endif


%{
#include "crpropa/module/NuclearDecay.h"
#include "crpropa/module/ElectronPairProduction.h"
#include "crpropa/module/PhotoPionProduction.h"
#include "crpropa/module/PhotoDisintegration.h"
#include "crpropa/module/ElasticScattering.h"
#include "crpropa/module/Redshift.h"
#include "crpropa/module/BreakCondition.h"
#include "crpropa/module/Boundary.h"
#include "crpropa/module/Observer.h"
#include "crpropa/module/TextOutput.h"
#include "crpropa/module/ParticleCollector.h"
#include "crpropa/module/HDF5Output.h"
#include "crpropa/module/OutputShell.h"
#include "crpropa/module/OutputROOT.h"
#include "crpropa/module/OutputCRPropa2.h"
#include "crpropa/module/EMCascade.h"
#include "crpropa/module/PhotonEleCa.h"
#include "crpropa/module/PhotonOutput1D.h"
#include "crpropa/module/SimplePropagation.h"
#include "crpropa/module/PropagationCK.h"
#include "crpropa/module/EMPairProduction.h"
#include "crpropa/module/EMDoublePairProduction.h"
#include "crpropa/module/EMTripletPairProduction.h"
#include "crpropa/module/EMInverseComptonScattering.h"
#include "crpropa/module/SynchrotronRadiation.h"
#include "crpropa/module/Tools.h"
#include "crpropa/module/DiffusionSDE.h"

#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/magneticField/MagneticFieldGrid.h"
#include "crpropa/magneticField/QuimbyMagneticField.h"
#include "crpropa/magneticField/AMRMagneticField.h"
#include "crpropa/magneticField/JF12Field.h"
#include "crpropa/magneticField/PshirkovField.h"

#include "crpropa/Referenced.h"
#include "crpropa/Candidate.h"
#include "crpropa/EmissionMap.h"
#include "crpropa/ParticleState.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Module.h"
#include "crpropa/ModuleList.h"
#include "crpropa/Random.h"
#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Source.h"
#include "crpropa/Common.h"
#include "crpropa/Cosmology.h"
#include "crpropa/PhotonBackground.h"
#include "crpropa/PhotonPropagation.h"
#include "crpropa/Grid.h"
#include "crpropa/GridTools.h"

#include "crpropa/Version.h"
%}

%{
#include <iostream>
#include <iomanip>
%}

%ignore operator<<;
%ignore operator>>;
%ignore *::operator=;
%ignore operator crpropa::Source*;
%ignore operator crpropa::SourceList*;
%ignore operator crpropa::SourceInterface*;
%ignore operator crpropa::SourceFeature*;
%ignore operator crpropa::Candidate*;
%ignore operator crpropa::Module*;
%ignore operator crpropa::ModuleList*;
%ignore operator crpropa::Observer*;
%ignore operator crpropa::ObserverFeature*;
%ignore operator crpropa::MagneticField*;
%ignore operator crpropa::ParticleCollector*;
%ignore crpropa::TextOutput::load;

%feature("ref")   crpropa::Referenced "$this->addReference();"
%feature("unref") crpropa::Referenced "$this->removeReference();"


%include "crpropa/Vector3.h"


%include "crpropa/Referenced.h"
%include "crpropa/Units.h"
%include "crpropa/Common.h"
%include "crpropa/Cosmology.h"
%include "crpropa/PhotonBackground.h"
%include "crpropa/PhotonPropagation.h"
%include "crpropa/Random.h"
%include "crpropa/ParticleState.h"
%include "crpropa/ParticleID.h"
%include "crpropa/ParticleMass.h"

/* override Candidate::getProperty() */
%ignore crpropa::Candidate::getProperty(const std::string &, std::string &) const;

%nothread; /* disable threading for extend*/
%extend crpropa::Candidate {
    PyObject * getProperty(PyObject * name){

        std::string value;
        std::string input;

        if (PyString_Check( name )){
            input = PyString_AsString( name );
        } else {
            std::cerr << "ERROR: The argument of getProperty() must be a string!" << std::endl;
            return NULL;
        }
        $self->getProperty( input, value ); 

        return PyString_FromString( value.c_str() );
    }
}; 
%thread; /* reenable threading */

%template(CandidateVector) std::vector< crpropa::ref_ptr<crpropa::Candidate> >;
%template(CandidateRefPtr) crpropa::ref_ptr<crpropa::Candidate>;
%include "crpropa/Candidate.h"

%template(ModuleRefPtr) crpropa::ref_ptr<crpropa::Module>;
%template(stdModuleList) std::list< crpropa::ref_ptr<crpropa::Module> >;
%feature("director") crpropa::Module;
%feature("director") crpropa::AbstractCondition;
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

%include "crpropa/EmissionMap.h"
%implicitconv crpropa::ref_ptr<crpropa::EmissionMap>;
%template(EmissionMapRefPtr) crpropa::ref_ptr<crpropa::EmissionMap>;
%implicitconv crpropa::ref_ptr<crpropa::CylindricalProjectionMap>;
%template(CylindricalProjectionMapRefPtr) crpropa::ref_ptr<crpropa::CylindricalProjectionMap>;

%include "crpropa/magneticField/MagneticFieldGrid.h"
%feature("notabstract") QuimbyMagneticFieldAdapter;
%include "crpropa/magneticField/QuimbyMagneticField.h"
%include "crpropa/magneticField/AMRMagneticField.h"
%include "crpropa/magneticField/JF12Field.h"
%include "crpropa/magneticField/PshirkovField.h"
%include "crpropa/module/BreakCondition.h"
%include "crpropa/module/Boundary.h"

%feature("director") crpropa::Observer;
%feature("director") crpropa::ObserverFeature;
%include "crpropa/module/Observer.h"
%include "crpropa/module/SimplePropagation.h"
%include "crpropa/module/PropagationCK.h"
%include "crpropa/module/Output.h"
%include "crpropa/module/DiffusionSDE.h"
%include "crpropa/module/TextOutput.h"
%inline %{
class RangeError {};
%}

%template(ParticleCollectorRefPtr) crpropa::ref_ptr<crpropa::ParticleCollector>;

%include "crpropa/module/ParticleCollector.h"
%include "crpropa/module/HDF5Output.h"
%include "crpropa/module/OutputShell.h"
%include "crpropa/module/OutputROOT.h"
%include "crpropa/module/OutputCRPropa2.h"
%include "crpropa/module/EMCascade.h"
%include "crpropa/module/PhotonEleCa.h"
%include "crpropa/module/PhotonOutput1D.h"
%include "crpropa/module/NuclearDecay.h"
%include "crpropa/module/ElectronPairProduction.h"
%include "crpropa/module/PhotoPionProduction.h"
%include "crpropa/module/PhotoDisintegration.h"
%include "crpropa/module/ElasticScattering.h"
%include "crpropa/module/Redshift.h"
%include "crpropa/module/EMPairProduction.h"
%include "crpropa/module/EMDoublePairProduction.h"
%include "crpropa/module/EMTripletPairProduction.h"
%include "crpropa/module/EMInverseComptonScattering.h"
%include "crpropa/module/SynchrotronRadiation.h"

%template(IntSet) std::set<int>;
%include "crpropa/module/Tools.h"

%template(SourceInterfaceRefPtr) crpropa::ref_ptr<crpropa::SourceInterface>;
%feature("director") crpropa::SourceInterface;
%template(SourceFeatureRefPtr) crpropa::ref_ptr<crpropa::SourceFeature>;
%feature("director") crpropa::SourceFeature;
%include "crpropa/Source.h"

%template(ModuleListRefPtr) crpropa::ref_ptr<crpropa::ModuleList>;
%include "crpropa/ModuleList.h"

%exception crpropa::ParticleCollector::__getitem__ {
  try {
        $action
  }
  catch (RangeError) {
        SWIG_exception(SWIG_IndexError, "Index out of bounds");
        return NULL;
  }

}

%extend crpropa::ParticleCollector {
  crpropa::ref_ptr<crpropa::Candidate> __getitem__(size_t i) {
        if (i >= $self->getCount()) {
                throw RangeError();
        }
        return (*($self))[i];
  }
};

