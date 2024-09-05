/* 2: CRPropa headers and Python extensions */

/* Python slots */
%feature("python:slot", "sq_length", functype="lenfunc") __len__;
%feature("python:slot", "mp_subscript", functype="binaryfunc") __getitem__;
%feature("python:slot", "tp_iter", functype="unaryfunc") __iter__;
%feature("python:slot", "tp_iternext", functype="iternextfunc") __next__;

/* Include headers */

#ifdef CRPROPA_HAVE_QUIMBY
  %import (module="quimby") "quimby/Referenced.h"
  %import (module="quimby") "quimby/Vector3.h"
  %import (module="quimby") "quimby/MagneticField.h"
  //%import (module="quimby") quimby.i
#endif


%{
#include "CRPropa.h"
  // for usage of namespace in header files, necessary for keyword arguments with units
  using namespace crpropa;
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
%ignore operator crpropa::PhotonField*;
%ignore operator crpropa::AdvectionField*;
%ignore operator crpropa::ParticleCollector*;
%ignore operator crpropa::Density*;
%ignore operator crpropa::CylindricalProjectionMap*;
%ignore operator crpropa::EmissionMap*;
%ignore operator crpropa::Grid< crpropa::Vector3< float > >*;
%ignore operator crpropa::Grid< crpropa::Vector3< double > >*;
%ignore operator crpropa::Grid< float >*;
%ignore operator crpropa::Grid< double >*;
%ignore crpropa::TextOutput::load;

%feature("ref")   crpropa::Referenced "$this->addReference();"
%feature("unref") crpropa::Referenced "$this->removeReference();"


%include "crpropa/Logging.h"

/* ignore public references and replace with attributes for Vector3d and Vector3f*/
%attribute(crpropa::Vector3<double>, double, x, getX, setX);
%attribute(crpropa::Vector3<double>, double, y, getY, setY);
%attribute(crpropa::Vector3<double>, double, z, getZ, setZ);
%attribute(crpropa::Vector3<float>, float, x, getX, setX);
%attribute(crpropa::Vector3<float>, float, y, getY, setY);
%attribute(crpropa::Vector3<float>, float, z, getZ, setZ);

/* implement array interface for numpy compatibility */
%feature("python:slot", "sq_length", functype="lenfunc") crpropa::Vector3::__len__;
%feature("python:slot", "mp_subscript", functype="binaryfunc") crpropa::Vector3::__getitem__;
%feature("python:slot", "mp_ass_subscript", functype="objobjargproc") crpropa::Vector3::__setitem__;
%typemap(directorin,numinputs=1) (const double *v) {
  npy_intp dim = 3;
  $input = PyArray_SimpleNewFromData(1, &dim, NPY_DOUBLE, (void *)$1);
}
%ignore crpropa::Vector3::data;

%include "crpropa/Vector3.h"

%exception crpropa::Vector3::__getitem__ {
  try {
    $action
  } catch (RangeError) {
    SWIG_exception(SWIG_IndexError, "Index out of bounds");
    return NULL;
  }
}

%exception crpropa::Vector3::__setitem__ {
  try {
    $action
  } catch (RangeError) {
    SWIG_exception(SWIG_IndexError, "Index out of bounds");
    return NULL;
  }
}

%extend crpropa::Vector3 {
  size_t __len__() {
    return 3;
  }

  PyObject* __array__() {
    npy_intp shape[1];
    shape[0] = 3;
    PyObject *ro;
    if (sizeof($self->data[0]) == NPY_SIZEOF_FLOAT) {
      ro = PyArray_SimpleNewFromData(1, shape, NPY_FLOAT, $self->data);
    } else if (sizeof($self->data[0]) == NPY_SIZEOF_DOUBLE) {
      ro = PyArray_SimpleNewFromData(1, shape, NPY_DOUBLE, $self->data);
    } else {
      KISS_LOG_ERROR << "crpropa::Vector3 has fixed size of 3 elements!";
    }

    return ro;
  }

  double __getitem__(size_t i) {
    if(i > 2) {
        throw RangeError();
    }

    return $self->data[i];
  }

  int __setitem__(size_t i, T value) {
    if(i > 2) {
      throw RangeError();
    }

    $self->data[i] = value;
    return 0;
  }

  const std::string getDescription() {
    char buffer[256];
    sprintf( buffer, "Vector(%.6G, %.6G, %.6G)", $self->x, $self->y, $self->z );
    return buffer;
  }
}

%feature("python:slot", "tp_str", functype="reprfunc") crpropa::Vector3::getDescription();
%feature("python:slot", "tp_repr", functype="reprfunc") crpropa::Vector3::getDescription();

%template(Vector3d) crpropa::Vector3<double>;
%template(Vector3f) crpropa::Vector3<float>;

%include "crpropa/Referenced.h"
%include "crpropa/Units.h"
%include "crpropa/Common.h"
%include "crpropa/Cosmology.h"
%template(RandomSeed) std::vector<uint32_t>;
%template(RandomSeedThreads) std::vector< std::vector<uint32_t> >;
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
  PyObject * getProperty(PyObject* name) {

    std::string input;

    if (PyUnicode_Check(name)) {
      input = PyUnicode_AsUTF8(name);
    } else if (PyString_Check(name)){
      input = PyString_AsString(name);
    } else {
      std::cerr << "ERROR: The argument of getProperty() must be a string/unicode object!" << std::endl;
      return NULL;
    }

    crpropa::Variant value = $self->getProperty(input);

    // implement this conversion here and not in the Variant as
    // __asPythonObject, as extensions cannot be called from extension.
    if (! value.isValid())  {
      Py_INCREF(Py_None);
      return Py_None;
    } else if (value.getTypeInfo() == typeid(bool)) {
      if(value.toBool()) {
        Py_RETURN_TRUE;
      } else {
        Py_RETURN_FALSE;
      }
    } else if (value.getTypeInfo() == typeid(char)) {
      // convert all integer types to python long
      return PyInt_FromLong(value.toInt64());
    } else if (value.getTypeInfo() == typeid(unsigned char)) {
      return PyInt_FromLong(value.toInt64());
    } else if (value.getTypeInfo() == typeid(int16_t)) {
      return PyInt_FromLong(value.toInt64());
    } else if (value.getTypeInfo() == typeid(uint16_t)) {
      return PyInt_FromLong(value.toInt64());
    } else if (value.getTypeInfo() == typeid(int32_t)) {
      return PyInt_FromLong(value.toInt64());
    } else if (value.getTypeInfo() == typeid(uint32_t)) {
      return PyInt_FromLong(value.toInt64());
    } else if (value.getTypeInfo() == typeid(int64_t)) {
      return PyLong_FromLong(value.toInt64());
    } else if (value.getTypeInfo() == typeid(uint64_t)) {
      return PyLong_FromUnsignedLong(value.toInt64());
    } else if (value.getTypeInfo() == typeid(float)) {
      // convert float and double to pyfloat which is double precision
      return PyFloat_FromDouble(value.toDouble());
    } else if (value.getTypeInfo() == typeid(double)) {
      return PyFloat_FromDouble(value.toDouble());
    } else if (value.getTypeInfo() == typeid(std::string)) {
      return PyUnicode_FromString(value.toString().c_str());
    }

    std::cerr << "ERROR: Unknown Type" << std::endl;
    return NULL;
  }

  PyObject * setProperty(PyObject * name, PyObject * value) {

    std::string input;

    if (PyUnicode_Check(name)){
      input = PyUnicode_AsUTF8(name);
    } else {
      std::cerr << "ERROR: The argument of setProperty() must be a string/unicode object!" << std::endl;
      return NULL;
    }

    if (value == Py_None) {
      $self->setProperty(input, crpropa::Variant());
      Py_RETURN_TRUE;
    } else if (PyBool_Check(value)) {
      if(value == Py_True) {
        $self->setProperty(input, true);
      } else {
        $self->setProperty(input, false);
      }
      Py_RETURN_TRUE;
    } else if (PyInt_Check(value)) {
      $self->setProperty(input, crpropa::Variant::fromInt32(PyInt_AsLong(value)));
      Py_RETURN_TRUE;
    } else if (PyLong_Check(value)) {
      $self->setProperty(input, crpropa::Variant::fromUInt64(PyLong_AsLong(value)));
      Py_RETURN_TRUE;
    } else if (PyFloat_Check(value)) {
      $self->setProperty(input, crpropa::Variant::fromDouble(PyFloat_AsDouble(value)));
      Py_RETURN_TRUE;
    } else if (PyUnicode_Check(value)){
      $self->setProperty(input, PyUnicode_AsUTF8(value));
      Py_RETURN_TRUE;
    } else {
      PyObject *t = PyObject_Str(PyObject_Type(value));
      std::string ot = PyUnicode_AsUTF8(t);
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

%template(OutputRefPtr) crpropa::ref_ptr<Output>;
%feature("director") crpropa::Output;
%ignore crpropa::Output::dumpIndexList(std::vector<int>);
%include "crpropa/module/Output.h"

%implicitconv crpropa::ref_ptr<crpropa::MagneticField>;
%template(MagneticFieldRefPtr) crpropa::ref_ptr<crpropa::MagneticField>;
%feature("director") crpropa::MagneticField;
%include "crpropa/magneticField/MagneticField.h"

%implicitconv crpropa::ref_ptr<crpropa::PhotonField>;
%template(PhotonFieldRefPtr) crpropa::ref_ptr<crpropa::PhotonField>;
%feature("director") crpropa::PhotonField;
%include "crpropa/PhotonBackground.h"

%implicitconv crpropa::ref_ptr<crpropa::AdvectionField>;
%template(AdvectionFieldRefPtr) crpropa::ref_ptr<crpropa::AdvectionField>;
%feature("director") crpropa::AdvectionField;
%include "crpropa/advectionField/AdvectionField.h"

%implicitconv crpropa::ref_ptr<crpropa::Density>;
%template(DensityRefPtr) crpropa::ref_ptr<crpropa::Density>;
%feature("director") crpropa::Density;
%include "crpropa/massDistribution/Density.h"

%include "crpropa/Grid.h"
%include "crpropa/GridTools.h"

%template(Array3d) std::array<double, 3>;
%template(Array3f) std::array<float, 3>;

%implicitconv crpropa::ref_ptr<crpropa::Grid<crpropa::Vector3<float> > >;
%template(Grid3fRefPtr) crpropa::ref_ptr<crpropa::Grid<crpropa::Vector3<float> > >;
%template(Grid3f) crpropa::Grid<crpropa::Vector3<float> >;

%implicitconv crpropa::ref_ptr<crpropa::Grid<crpropa::Vector3<double> > >;
%template(Grid3dRefPtr) crpropa::ref_ptr<crpropa::Grid<crpropa::Vector3<double> > >;
%template(Grid3d) crpropa::Grid<crpropa::Vector3<double> >;

%implicitconv crpropa::ref_ptr<crpropa::Grid<float> >;
%template(Grid1fRefPtr) crpropa::ref_ptr<crpropa::Grid<float> >;
%template(Grid1f) crpropa::Grid<float>;

%implicitconv crpropa::ref_ptr<crpropa::Grid<double> >;
%template(Grid1dRefPtr) crpropa::ref_ptr<crpropa::Grid<double> >;
%template(Grid1d) crpropa::Grid<double>;

%implicitconv std::pair<std::vector<int>, std::vector<float> >;
%template(PairIntFloat) std::pair<int, float>;
%template(PairVector) std::vector<std::pair<int, float> >;

%include "crpropa/EmissionMap.h"
%implicitconv crpropa::ref_ptr<crpropa::EmissionMap>;
%template(EmissionMapRefPtr) crpropa::ref_ptr<crpropa::EmissionMap>;
%implicitconv crpropa::ref_ptr<crpropa::CylindricalProjectionMap>;
%template(CylindricalProjectionMapRefPtr) crpropa::ref_ptr<crpropa::CylindricalProjectionMap>;

%include "crpropa/magneticField/MagneticFieldGrid.h"
%include "crpropa/magneticField/GalacticMagneticField.h"
%feature("notabstract") QuimbyMagneticFieldAdapter;
%include "crpropa/magneticField/QuimbyMagneticField.h"
%include "crpropa/magneticField/JF12Field.h"
%include "crpropa/magneticField/JF12FieldSolenoidal.h"
%include "crpropa/magneticField/PolarizedSingleModeMagneticField.h"
%include "crpropa/magneticField/PT11Field.h"
%include "crpropa/magneticField/TF17Field.h"
%include "crpropa/magneticField/UF23Field.h"
%include "crpropa/magneticField/ArchimedeanSpiralField.h"
%include "crpropa/magneticField/CMZField.h"
%include "crpropa/magneticField/turbulentField/TurbulentField.h"
%include "crpropa/magneticField/turbulentField/GridTurbulence.h"
%include "crpropa/magneticField/turbulentField/SimpleGridTurbulence.h"
%include "crpropa/magneticField/turbulentField/HelicalGridTurbulence.h"
%include "crpropa/magneticField/turbulentField/PlaneWaveTurbulence.h"
%include "crpropa/module/BreakCondition.h"
%include "crpropa/module/Boundary.h"

%feature("director") crpropa::Observer;
%feature("director") crpropa::ObserverFeature;
%include "crpropa/module/Observer.h"
%include "crpropa/module/SimplePropagation.h"
%include "crpropa/module/PropagationCK.h"
%include "crpropa/module/PropagationBP.h"

%ignore crpropa::Output::enableProperty(const std::string &property, const Variant& defaultValue, const std::string &comment = "");
%extend crpropa::Output{
  PyObject* enableProperty(const std::string &name, PyObject* defaultValue, const std::string &comment="") {

    if (defaultValue == Py_None) {
      Py_RETURN_TRUE;
    } else if (PyBool_Check(defaultValue)) {
      if(defaultValue == Py_True) {
        $self->enableProperty(name, true, comment);
      } else {
        $self->enableProperty(name, false, comment);
      }
      Py_RETURN_TRUE;
    } else if (PyInt_Check(defaultValue)) {
      $self->enableProperty(name, crpropa::Variant::fromInt32(PyInt_AsLong(defaultValue)), comment);
      Py_RETURN_TRUE;
    } else if (PyLong_Check(defaultValue)) {
      $self->enableProperty(name, crpropa::Variant::fromInt64(PyLong_AsLong(defaultValue)), comment);
      Py_RETURN_TRUE;
    } else if (PyFloat_Check(defaultValue)) {
      $self->enableProperty(name, crpropa::Variant::fromDouble(PyFloat_AsDouble(defaultValue)), comment);
      Py_RETURN_TRUE;
    } else if (PyUnicode_Check(defaultValue)){
      std::string ss = PyUnicode_AsUTF8(defaultValue);
      $self->enableProperty(name, ss, comment);
      Py_RETURN_TRUE;
    } else {
      PyObject* t = PyObject_Str(PyObject_Type(defaultValue));
      std::string ot = PyUnicode_AsUTF8(t);
      std::cerr << "ERROR: Unknown Type: " << ot << std::endl;
      return NULL;
    }

  }
}

%include "crpropa/module/DiffusionSDE.h"
%include "crpropa/module/TextOutput.h"
%include "crpropa/module/HDF5Output.h"
%include "crpropa/module/OutputShell.h"
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
%include "crpropa/module/MomentumDiffusion.h"
%include "crpropa/module/CandidateSplitting.h"

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
          cur(_cur), end(_end) {
          }
      ModuleListIterator* __iter__() {
        return this;
      }
      crpropa::ModuleList::iterator cur;
      crpropa::ModuleList::iterator end;
    };
%}

%extend ModuleListIterator {
  crpropa::ref_ptr<crpropa::Module>& __next__() {
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
        cur(_cur), end(_end) {
        }
    ParticleCollectorIterator* __iter__() {
      return this;
    }
    crpropa::ParticleCollector::iterator cur;
    crpropa::ParticleCollector::iterator end;
  };
%}

%extend ParticleCollectorIterator {
  crpropa::ref_ptr<crpropa::Candidate>& __next__() {
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

  std::vector<crpropa::ref_ptr<crpropa::Candidate>> __getitem__(PyObject *param) {
    std::vector< crpropa::ref_ptr<crpropa::Candidate> > result;

    if (PySlice_Check(param)) {
      Py_ssize_t len = 0, start = 0, stop = 0, step = 0, slicelength = 0, i = 0;
      len = $self->size();

      PySlice_GetIndicesEx(param, len, &start, &stop, &step, &slicelength);

      for(crpropa::ParticleCollector::iterator itr = $self->begin(); itr != $self->end(); ++itr) {
        if(i >= start && i < stop) {
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
%include "crpropa/massDistribution/Density.h"
%include "crpropa/massDistribution/Nakanishi.h"
%include "crpropa/massDistribution/Cordes.h"
%include "crpropa/massDistribution/Ferriere.h"
%include "crpropa/massDistribution/Massdistribution.h"
%include "crpropa/massDistribution/ConstantDensity.h"


%template(StepLengthModifierRefPtr) crpropa::ref_ptr<crpropa::StepLengthModifier>;
%feature("director") crpropa::StepLengthModifier;
%include "crpropa/module/Acceleration.h"
