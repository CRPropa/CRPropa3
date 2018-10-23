/* 2: CRPropa headers and Python extensions */

/* Python slots */
%feature("python:slot", "sq_length", functype="lenfunc") __len__;
%feature("python:slot", "mp_subscript", functype="binaryfunc") __getitem__;
%feature("python:slot", "tp_iter", functype="unaryfunc") __iter__;
#ifdef SWIG_PYTHON3
%feature("python:slot", "tp_iternext", functype="iternextfunc") __next__;
#else
%feature("python:slot", "tp_iternext", functype="iternextfunc") next;
#endif

/* Include headers */

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
#include "CRPropa.h"
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
%ignore operator crpropa::AdvectionField*;
%ignore operator crpropa::ParticleCollector*;
%ignore crpropa::TextOutput::load;

%feature("ref")   crpropa::Referenced "$this->addReference();"
%feature("unref") crpropa::Referenced "$this->removeReference();"


%include "crpropa/Logging.h"
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
%include "crpropa/Version.h"

%import "crpropa/Variant.h"

/* override Candidate::getProperty() */
%ignore crpropa::Candidate::getProperty(const std::string &) const;

%nothread; /* disable threading for extend*/
%extend crpropa::Candidate {
    PyObject * getProperty(PyObject * name){

        std::string input;

        if (PyUnicode_Check(name)){
          #ifdef SWIG_PYTHON3
          // test on PY_MAJOR_VERSION >= 3 wont work with swig
              input = PyUnicode_AsUTF8(name);
          #else
              PyObject *s =  PyUnicode_AsUTF8String(name);
              input = PyString_AsString(s);
          #endif
        }
        #ifndef SWIG_PYTHON3
        else if (PyString_Check(name)){
            input = PyString_AsString(name);
        }
        #endif
        else {
            std::cerr << "ERROR: The argument of getProperty() must be a string/unicode object!" << std::endl;
            return NULL;
        }

        crpropa::Variant value = $self->getProperty(input);

        // implement this conversion here and not in the Variant as
        // __asPythonObject, as extensions cannot be called from extension.
        if (! value.isValid())
        {
          Py_INCREF(Py_None);
          return Py_None;
        }
        else if (value.getTypeInfo() == typeid(bool))
        {
         if(value.toBool())
         {
          Py_RETURN_TRUE;
         }
         else
         {
          Py_RETURN_FALSE;
         }
        }
        // convert all integer types to python long
        else if (value.getTypeInfo() == typeid(char))
        {
          return PyInt_FromLong(value.toInt64());
        }
        else if (value.getTypeInfo() == typeid(unsigned char))
        {
          return PyInt_FromLong(value.toInt64());
        }
        else if (value.getTypeInfo() == typeid(int16_t))
        {
          return PyInt_FromLong(value.toInt64());
        }
        else if (value.getTypeInfo() == typeid(uint16_t))
        {
          return PyInt_FromLong(value.toInt64());
        }
        else if (value.getTypeInfo() == typeid(int32_t))
        {
          return PyInt_FromLong(value.toInt64());
        }
        else if (value.getTypeInfo() == typeid(uint32_t))
        {
          return PyInt_FromLong(value.toInt64());
        }
        else if (value.getTypeInfo() == typeid(int64_t))
        {
          return PyLong_FromLong(value.toInt64());
        }
        else if (value.getTypeInfo() == typeid(uint64_t))
        {
          return PyLong_FromUnsignedLong(value.toInt64());
        }
        // convert float and double to pyfloat which is double precision
        else if (value.getTypeInfo() == typeid(float))
        {
          return PyFloat_FromDouble(value.toDouble());
        }
        else if (value.getTypeInfo() == typeid(double))
        {
          return PyFloat_FromDouble(value.toDouble());
        }
        else if (value.getTypeInfo() == typeid(std::string))
        {
        #ifdef SWIG_PYTHON3
          return PyUnicode_FromString(value.toString().c_str());
        #else
          return PyString_FromString(value.toString().c_str());
        #endif
        }

        std::cerr << "ERROR: Unknown Type" << std::endl;
        return NULL;
    }


    PyObject * setProperty(PyObject * name, PyObject * value){

        std::string input;

        if (PyUnicode_Check(name)){
          #ifdef SWIG_PYTHON3
              input = PyUnicode_AsUTF8(name);
          #else
              input = PyUnicode_AS_DATA(name);
              PyObject *s =  PyUnicode_AsUTF8String(name);
              input = PyString_AsString(s);
          #endif
        }
        #ifndef SWIG_PYTHON3
        else if (PyString_Check( name )){
            input = PyString_AsString( name );
        }
        #endif
        else {
            std::cerr << "ERROR: The argument of setProperty() must be a string/unicode object!" << std::endl;
            return NULL;
        }


        if (value == Py_None)
        {
          $self->setProperty(input, crpropa::Variant());
        Py_RETURN_TRUE;
        }
        else if (PyBool_Check(value))
        {
         if(value == Py_True)
         {
          $self->setProperty(input, true);
         }
         else
         {
          $self->setProperty(input, false);
         }
          Py_RETURN_TRUE;
        }
        else if (PyInt_Check(value))
        {
          $self->setProperty(input, crpropa::Variant::fromInt32(PyInt_AsLong(value)));
          Py_RETURN_TRUE;
        }
        else if (PyLong_Check(value))
        {
          $self->setProperty(input, crpropa::Variant::fromUInt64(PyLong_AsLong(value)));
          Py_RETURN_TRUE;
        }
        else if (PyFloat_Check(value))
        {
          $self->setProperty(input, crpropa::Variant::fromDouble(PyFloat_AsDouble(value)));
          Py_RETURN_TRUE;
        }
        else if (PyUnicode_Check(value)){
        #ifdef SWIG_PYTHON3
          $self->setProperty(input, PyUnicode_AsUTF8(value));
        #else
          PyObject *s =  PyUnicode_AsUTF8String(value);
          $self->setProperty(input, PyString_AsString(s));
        #endif
          Py_RETURN_TRUE;
        }
        #ifndef SWIG_PYTHON3
        else if (PyString_Check( value))
        {
          $self->setProperty(input, PyString_AsString(value));
          Py_RETURN_TRUE;
        }
        #endif
        else
        {
          PyObject *t = PyObject_Str(PyObject_Type(value));
          std::string ot;

          #ifdef SWIG_PYTHON3
            ot = PyUnicode_AsUTF8(t);
          #else
            ot = PyString_AsString(t);
          #endif
          std::cerr << "ERROR: Unknown Type: " << ot << std::endl;
          return NULL;
        }
    }
};
%thread; /* reenable threading */


%template(CandidateVector) std::vector< crpropa::ref_ptr<crpropa::Candidate> >;
%template(CandidateRefPtr) crpropa::ref_ptr<crpropa::Candidate>;
%include "crpropa/Candidate.h"

%feature("director") crpropa::Surface;
%feature("director") crpropa::ClosedSurface;
%include "crpropa/Geometry.h"

%template(ModuleRefPtr) crpropa::ref_ptr<crpropa::Module>;
%template(stdModuleList) std::list< crpropa::ref_ptr<crpropa::Module> >;
%feature("director") crpropa::Module;
%feature("director") crpropa::AbstractCondition;
%include "crpropa/Module.h"

%implicitconv crpropa::ref_ptr<crpropa::MagneticField>;
%template(MagneticFieldRefPtr) crpropa::ref_ptr<crpropa::MagneticField>;
%include "crpropa/magneticField/MagneticField.h"

%implicitconv crpropa::ref_ptr<crpropa::AdvectionField>;
%template(AdvectionFieldRefPtr) crpropa::ref_ptr<crpropa::AdvectionField>;
%include "crpropa/advectionField/AdvectionField.h"

%implicitconv crpropa::ref_ptr<crpropa::Density>;
%template(DensityRefPtr) crpropa::ref_ptr<crpropa::Density>;
%include "crpropa/Massdistribution/Density.h"

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
%include "crpropa/magneticField/JF12FieldSolenoidal.h"
%include "crpropa/magneticField/PT11Field.h"
%include "crpropa/magneticField/ArchimedeanSpiralField.h"
%include "crpropa/module/BreakCondition.h"
%include "crpropa/module/Boundary.h"

%feature("director") crpropa::Observer;
%feature("director") crpropa::ObserverFeature;
%include "crpropa/module/Observer.h"
%include "crpropa/module/SimplePropagation.h"
%include "crpropa/module/PropagationCK.h"

%ignore crpropa::Output::enableProperty(const std::string &property, const Variant& defaultValue, const std::string &comment = "");
%extend crpropa::Output{
  PyObject * enableProperty(const std::string &name, PyObject* defaultValue, const std::string &comment="")
  {

       if (defaultValue == Py_None)
        {
          Py_RETURN_TRUE;
        }
        else if (PyBool_Check(defaultValue))
        {
         if(defaultValue == Py_True)
         {
          $self->enableProperty(name, true, comment);
         }
         else
         {
          $self->enableProperty(name, false, comment);
         }
          Py_RETURN_TRUE;
        }
        else if (PyInt_Check(defaultValue))
        {
          $self->enableProperty(name, crpropa::Variant::fromInt32(PyInt_AsLong(defaultValue)), comment);
          Py_RETURN_TRUE;
        }
        else if (PyLong_Check(defaultValue))
        {
          $self->enableProperty(name, crpropa::Variant::fromInt64(PyLong_AsLong(defaultValue)), comment);
          Py_RETURN_TRUE;
        }
        else if (PyFloat_Check(defaultValue))
        {
          $self->enableProperty(name, crpropa::Variant::fromDouble(PyFloat_AsDouble(defaultValue)), comment);
          Py_RETURN_TRUE;
        }
        else if (PyUnicode_Check(defaultValue)){
        #ifdef SWIG_PYTHON3
          std::string ss = PyUnicode_AsUTF8(defaultValue);
        #else
          PyObject *s =  PyUnicode_AsUTF8String(defaultValue);
          std::string ss = PyString_AsString(s);
        #endif
          $self->enableProperty(name, ss, comment);
          Py_RETURN_TRUE;
        }
        #ifndef SWIG_PYTHON3
        else if (PyString_Check( defaultValue))
        {
          std::string ss = PyString_AsString(defaultValue);
          $self->enableProperty(name, ss, comment);
          Py_RETURN_TRUE;
        }
        #endif
        else
        {
          PyObject *t = PyObject_Str(PyObject_Type(defaultValue));
          std::string ot;

          #ifdef SWIG_PYTHON3
            ot = PyUnicode_AsUTF8(t);
          #else
            ot = PyString_AsString(t);
          #endif
          std::cerr << "ERROR: Unknown Type: " << ot << std::endl;
          return NULL;
        }

  }
}


%include "crpropa/module/Output.h"
%include "crpropa/module/DiffusionSDE.h"
%include "crpropa/module/TextOutput.h"

%include "crpropa/module/HDF5Output.h"
%include "crpropa/module/OutputShell.h"
%include "crpropa/module/EMCascade.h"
%include "crpropa/module/PhotonEleCa.h"
%include "crpropa/module/PhotonOutput1D.h"
%include "crpropa/module/NuclearDecay.h"
%include "crpropa/module/ElectronPairProduction.h"
%include "crpropa/module/PhotoPionProduction.h"
%include "crpropa/module/PhotoDisintegration.h"
%include "crpropa/module/ElasticScattering.h"
%include "crpropa/module/Redshift.h"
%include "crpropa/module/RestrictToRegion.h"
%include "crpropa/module/EMPairProduction.h"
%include "crpropa/module/EMDoublePairProduction.h"
%include "crpropa/module/EMTripletPairProduction.h"
%include "crpropa/module/EMInverseComptonScattering.h"
%include "crpropa/module/SynchrotronRadiation.h"
%include "crpropa/module/AdiabaticCooling.h"

%template(IntSet) std::set<int>;
%include "crpropa/module/Tools.h"

%template(SourceInterfaceRefPtr) crpropa::ref_ptr<crpropa::SourceInterface>;
%feature("director") crpropa::SourceInterface;
%template(SourceFeatureRefPtr) crpropa::ref_ptr<crpropa::SourceFeature>;
%feature("director") crpropa::SourceFeature;
%include "crpropa/Source.h"

%inline %{
class ModuleListIterator {
  public:
        ModuleListIterator(
                crpropa::ModuleList::iterator _cur,
                crpropa::ModuleList::iterator _end) :
                        cur(_cur), end(_end) {}
        ModuleListIterator* __iter__() { return this; }
        crpropa::ModuleList::iterator cur;
        crpropa::ModuleList::iterator end;
  };
%}

%extend ModuleListIterator {
#ifdef SWIG_PYTHON3
  crpropa::ref_ptr<crpropa::Module>& __next__() {
#else
  crpropa::ref_ptr<crpropa::Module>& next() {
#endif
    if ($self->cur != $self->end) {
        return *$self->cur++;
    }
    throw StopIterator();
  }
}

%extend crpropa::ModuleList {
  ModuleListIterator __iter__() {
        return ModuleListIterator($self->begin(), $self->end());
  }
  crpropa::ref_ptr<crpropa::Module> __getitem__(size_t i) {
        if (i >= $self->size()) {
                throw RangeError();
        }
        return (*($self))[i];
  }
  size_t __len__() {
        return $self->size();
  }
};

%template(ModuleListRefPtr) crpropa::ref_ptr<crpropa::ModuleList>;
%include "crpropa/ModuleList.h"

%template(ParticleCollectorRefPtr) crpropa::ref_ptr<crpropa::ParticleCollector>;

%inline %{
class ParticleCollectorIterator {
  public:
        ParticleCollectorIterator(
                crpropa::ParticleCollector::iterator _cur,
                crpropa::ParticleCollector::iterator _end) :
                        cur(_cur), end(_end) {}
        ParticleCollectorIterator* __iter__() { return this; }
        crpropa::ParticleCollector::iterator cur;
        crpropa::ParticleCollector::iterator end;
  };
%}

%extend ParticleCollectorIterator {
#ifdef SWIG_PYTHON3
  crpropa::ref_ptr<crpropa::Candidate>& __next__() {
#else
  crpropa::ref_ptr<crpropa::Candidate>& next() {
#endif
    if ($self->cur != $self->end) {
        return *$self->cur++;
    }
    throw StopIterator();
  }
}

%extend crpropa::ParticleCollector {
  ParticleCollectorIterator __iter__() {
        return ParticleCollectorIterator($self->begin(), $self->end());
  }
  crpropa::ref_ptr<crpropa::Candidate> __getitem__(size_t i) {
        if (i >= $self->size()) {
                throw RangeError();
        }
        return (*($self))[i];
  }
  std::vector< crpropa::ref_ptr<crpropa::Candidate> > __getitem__(PyObject *param) {
        std::vector< crpropa::ref_ptr<crpropa::Candidate> > result;

        if (PySlice_Check(param)) {
                Py_ssize_t len = 0, start = 0, stop = 0, step = 0, slicelength = 0, i = 0;
                len = $self->size();

                #ifdef SWIG_PYTHON3
                    PySlice_GetIndicesEx(param, len, &start, &stop, &step, &slicelength);
                #else
                    PySlice_GetIndicesEx((PySliceObject*)param, len, &start, &stop, &step, &slicelength);
                #endif

                for(crpropa::ParticleCollector::iterator itr = $self->begin(); itr != $self->end(); ++itr){
                        if( i >= start && i < stop){
                                result.push_back(itr->get());
                        }
                        ++i;
                }
                return result;
        } else {
                throw RangeError();
        }
  }
  size_t __len__() {
        return $self->size();
  }
};

%include "crpropa/module/ParticleCollector.h"

%include "crpropa/Massdistribution/Density.h"
%include "crpropa/Massdistribution/Nakanishi.h"
%include "crpropa/Massdistribution/Cordes.h"
%include "crpropa/Massdistribution/Ferriere.h"
%include "crpropa/Massdistribution/Pohl.h"
%include "crpropa/Massdistribution/Massdistribution.h"
%include "crpropa/Massdistribution/ConstantDensity.h"

