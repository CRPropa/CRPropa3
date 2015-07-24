/* CRPropa3 SWIG interface (for Python) */

/* Content:
 *
 * 1. SWIG settings and workarounds
 * 2. SWIG and CRPropa headers
 * 3. Pretty print for Python
 * 4. Magnetic Lens
 * 5. Particle Maps Container
 *
 */


/* 1. SWIG settings and workarounds */

%module(directors="1", threads="1") crpropa

%feature("autodoc", "1"); // automatic docstrings

%{
// workaround for SWIG < 2.0.5 with GCC >= 4.7
#include <cstddef>
using std::ptrdiff_t;
%}

%exception
{
 try
 {
   $action
 }
 catch (Swig::DirectorException &e) {
   SWIG_exception(SWIG_RuntimeError, e.getMessage());
 }
 catch (const std::exception& e) {
   SWIG_exception(SWIG_RuntimeError, e.what());
 }
 catch (const char *e) {
   SWIG_exception(SWIG_RuntimeError, e);
 }
}


/* 2: SWIG and CRPropa headers */

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
#include "crpropa/module/PhotonDINT1D.h"
#include "crpropa/module/PhotonEleCa.h"
#include "crpropa/module/PhotonOutput1D.h"
#include "crpropa/module/SimplePropagation.h"
#include "crpropa/module/PropagationCK.h"
#include "crpropa/module/Tools.h"

#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/magneticField/MagneticFieldGrid.h"
#include "crpropa/magneticField/QuimbyMagneticField.h"
#include "crpropa/magneticField/AMRMagneticField.h"
#include "crpropa/magneticField/JF12Field.h"
#include "crpropa/magneticField/TurbulentMagneticField.h"

#include "crpropa/Referenced.h"
#include "crpropa/Candidate.h"
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
%ignore operator crpropa::MagneticField*;

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
%include "crpropa/magneticField/AMRMagneticField.h"
%include "crpropa/magneticField/JF12Field.h"
%include "crpropa/magneticField/TurbulentMagneticField.h"

%include "crpropa/module/BreakCondition.h"
%include "crpropa/module/Boundary.h"
%include "crpropa/module/Observer.h"
%include "crpropa/module/SimplePropagation.h"
%include "crpropa/module/PropagationCK.h"
%include "crpropa/module/OutputTXT.h"
%include "crpropa/module/OutputShell.h"
%include "crpropa/module/OutputROOT.h"
%include "crpropa/module/OutputCRPropa2.h"
%include "crpropa/module/PhotonDINT.h"
%include "crpropa/module/PhotonDINT1D.h"
%include "crpropa/module/PhotonEleCa.h"
%include "crpropa/module/PhotonOutput1D.h"
%include "crpropa/module/ElectronPairProduction.h"
%include "crpropa/module/NuclearDecay.h"
%include "crpropa/module/PhotoPionProduction.h"
%include "crpropa/module/PhotoDisintegration.h"
%include "crpropa/module/Redshift.h"
%include "crpropa/module/Tools.h"

%template(SourceInterfaceRefPtr) crpropa::ref_ptr<crpropa::SourceInterface>;
%feature("director") crpropa::SourceInterface;
%template(SourceFeatureRefPtr) crpropa::ref_ptr<crpropa::SourceFeature>;
%feature("director") crpropa::SourceFeature;
%include "crpropa/Source.h"

%template(ModuleListRefPtr) crpropa::ref_ptr<crpropa::ModuleList>;
%include "crpropa/ModuleList.h"

%include "crpropa/Version.h"
%pythoncode %{
    __version__ = g_GIT_DESC 
%}

/* 3. Pretty print for Python */
/*  __repr__ << getDescription */

%define __REPR__( classname ) 
%feature("python:slot", "tp_str", functype="reprfunc") classname::repr();
%feature("python:slot", "tp_repr", functype="reprfunc") classname::repr();

%extend classname {
        const std::string repr() {
            return $self->getDescription();
        }
}

%enddef

%define VECTOR3__REPR__( classname ) 
%feature("python:slot", "tp_str", functype="reprfunc") classname::repr();
%feature("python:slot", "tp_repr", functype="reprfunc") classname::repr();

%extend classname {
        const std::string repr() {
            char buffer[1024];
            sprintf( buffer, "Vector( %.5e, %.5e, %.5e )", $self->x, $self->y, $self->z );
            return buffer;
        }
}

%enddef

__REPR__( crpropa::ParticleState );
__REPR__( crpropa::Candidate );
__REPR__( crpropa::Module );
__REPR__( crpropa::ModuleList );

__REPR__( crpropa::Source );

__REPR__( crpropa::SourceList );
__REPR__( crpropa::SourceParticleType );
__REPR__( crpropa::SourceMultipleParticleTypes );
__REPR__( crpropa::SourceEnergy );
__REPR__( crpropa::SourcePowerLawSpectrum );
__REPR__( crpropa::SourceMultiplePositions );
__REPR__( crpropa::SourceComposition );
__REPR__( crpropa::SourcePosition );
__REPR__( crpropa::SourceUniform1D );
__REPR__( crpropa::SourceUniformBox );
__REPR__( crpropa::SourceUniformShell );
__REPR__( crpropa::SourceUniformSphere );
__REPR__( crpropa::SourceDensityGrid );
__REPR__( crpropa::SourceDensityGrid1D );
__REPR__( crpropa::SourceDirection );
__REPR__( crpropa::SourceIsotropicEmission );
__REPR__( crpropa::SourceEmissionCone );
__REPR__( crpropa::SourceRedshift );
__REPR__( crpropa::SourceRedshift1D );
__REPR__( crpropa::SourceUniformRedshift );

__REPR__( crpropa::Observer );
__REPR__( crpropa::ObserverPoint );
__REPR__( crpropa::ObserverSmallSphere );
__REPR__( crpropa::ObserverLargeSphere );
__REPR__( crpropa::ObserverRedshiftWindow );
__REPR__( crpropa::ObserverNucleusVeto );
__REPR__( crpropa::ObserverNeutrinoVeto );
__REPR__( crpropa::ObserverPhotonVeto );

VECTOR3__REPR__( crpropa::Vector3 );

%template(Vector3d) crpropa::Vector3<double>;
%template(Vector3f) crpropa::Vector3<float>;

%pythoncode %{
    DeflectionCK = PropagationCK  # legacy name
%}


/* 4. Magnetic Lens */

#ifdef WITHNUMPY
%{
/* Include numpy array interface, if available */
  #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
  #include "numpy/arrayobject.h"
  #include "numpy/ufuncobject.h"
%}
#endif

/* Initialize numpy array interface, if available */
#ifdef WITHNUMPY
%init %{
import_array();
import_ufunc();
%}

%pythoncode %{
import numpy
__WITHNUMPY = True
%}

#else
%pythoncode %{
__WITHNUMPY = False 
%}
#endif


#ifdef WITH_GALACTIC_LENSES

%include typemaps.i

%template(IntVector) std::vector<int>;
%template(DoubleVector) std::vector<double>;

%{
#include "crpropa/magneticLens/ModelMatrix.h"
#include "crpropa/magneticLens/Pixelization.h"
#include "crpropa/magneticLens/MagneticLens.h"
#include "crpropa/magneticLens/ParticleMapsContainer.h"
%}

%include "crpropa/magneticLens/ModelMatrix.h"
%apply double &INOUT {double &longitude, double &latitude};
%typemap(in,numinputs=0) double& longitude (double temp) "$1 = &temp;"
%typemap(in,numinputs=0) double& latitude (double temp) "$1 = &temp;"
%typemap(argout) double& longitude {
  %append_output(PyFloat_FromDouble(*$1));
}
%typemap(argout) double& latitude {
  %append_output(PyFloat_FromDouble(*$1));
}



%include "crpropa/magneticLens/Pixelization.h"

%ignore crpropa::Pixelization::nPix( uint8_t order );

%pythoncode %{
class Pixelization(Pixelization):
    def nPix(self, order=None):
      if order == None:
        return Pixelization_nPix(self.getOrder())
      else:
        return Pixelization_nPix(order)
%}



%apply double &INOUT {double &phi, double &theta};
%ignore MagneticLens::transformModelVector(double *,double) const;
%include "crpropa/magneticLens/MagneticLens.h"
%template(LenspartVector) std::vector< crpropa::LensPart *>;

#ifdef WITHNUMPY
%extend crpropa::MagneticLens{
    PyObject * transformModelVector_numpyArray(PyObject *input, double rigidity)
    {
      PyArrayObject *arr = NULL;
      PyArray_Descr *dtype = NULL;
      int ndim = 0;
      npy_intp dims[NPY_MAXDIMS];
      if (PyArray_GetArrayParamsFromObject(input, NULL, 1, &dtype, &ndim, dims, &arr, NULL) < 0)
      {
        return NULL; 
      }

      if (arr == NULL) 
      {
        return NULL;
      }

      double *dataPointer = (double*) PyArray_DATA(arr);
      $self->transformModelVector(dataPointer, rigidity);
      return input;
    }
};
#else
%extend crpropa::MagneticLens{
    PyObject * transformModelVector_numpyArray(PyObject *input, double rigidity)
    {
      std::cerr << "ERROR: PARSEC was compiled without numpy support!" << std::endl;
      return NULL;
    }
};
#endif

%pythoncode %{
class MagneticLens(MagneticLens):
        transformModelVector = MagneticLens.transformModelVector_numpyArray
%}


/* 5. Particle Maps Container */

%ignore ParticleMapsContainer::getMap;
%ignore ParticleMapsContainer::getParticleIds;
%ignore ParticleMapsContainer::getEnergies;
%ignore ParticleMapsContainer::getRandomParticles;
%include "crpropa/magneticLens/ParticleMapsContainer.h"

#ifdef WITHNUMPY
%extend crpropa::ParticleMapsContainer{
  PyObject *getMap_numpyArray(const int particleId, double energy)
  {
      double* data = $self->getMap(particleId, energy);
      npy_intp npix = $self->getNumberOfPixels();
      npy_intp dims[1] = {npix};
      return PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, (void*)data);
  }

  PyObject *getParticleIds_numpyArray()
  {
      std::vector<int> v = $self->getParticleIds();
      npy_intp size = v.size();
      PyObject *out = PyArray_SimpleNew(1, &size, NPY_INT);
      memcpy(PyArray_DATA((PyArrayObject *) out), &v[0], v.size() * sizeof(int));
      return out; 
  }

  PyObject *getEnergies_numpyArray(const int pid)
  {
      std::vector<double> v = $self->getEnergies(pid);
      npy_intp size = v.size();
      PyObject *out = PyArray_SimpleNew(1, &size, NPY_DOUBLE);
      memcpy(PyArray_DATA((PyArrayObject *) out), &v[0], v.size() * sizeof(double));
      return out; 
  }

  PyObject *getRandomParticles_numpyArray(size_t N)
  {
      vector<int> particleId;
			vector<double> energy;
      vector<double> galacticLongitudes;
			vector<double> galacticLatitudes;
      $self->getRandomParticles(N, particleId, energy, galacticLongitudes,
          galacticLatitudes);
      
      npy_intp size = N;
      PyObject *oId = PyArray_SimpleNew(1, &size, NPY_INT);
      PyObject *oEnergy = PyArray_SimpleNew(1, &size, NPY_DOUBLE);
      PyObject *oLon = PyArray_SimpleNew(1, &size, NPY_DOUBLE);
      PyObject *oLat = PyArray_SimpleNew(1, &size, NPY_DOUBLE);

      memcpy(PyArray_DATA((PyArrayObject *) oId), &particleId[0],
          particleId.size() * sizeof(int));
      memcpy(PyArray_DATA((PyArrayObject *) oEnergy), &energy[0], energy.size()
          * sizeof(double));
      memcpy(PyArray_DATA((PyArrayObject *) oLon), &galacticLongitudes[0],
          galacticLongitudes.size() * sizeof(double));
      memcpy(PyArray_DATA((PyArrayObject *) oLat), &galacticLatitudes[0],
          galacticLatitudes.size() * sizeof(double));

      PyObject *returnList = PyList_New(4);
      PyList_SET_ITEM(returnList, 0, oId);
      PyList_SET_ITEM(returnList, 1, oEnergy);
      PyList_SET_ITEM(returnList, 2, oLon);
      PyList_SET_ITEM(returnList, 3, oLat);
      return returnList;
  }

};
#else // with numpy
%extend crpropa::ParticleMapsContainer{
  PyObject *getMap_numpyArray(const int particleId, double energy)
  {
      std::cerr << "ERROR: PARSEC was compiled without numpy support!" << std::endl;
      return NULL;
  }
};
%extend crpropa::ParticleMapsContainer{
  PyObject *getParticleIds_numpyArray(const int particleId, double energy)
  {
      std::cerr << "ERROR: PARSEC was compiled without numpy support!" << std::endl;
      return NULL;
  }
};
%extend crpropa::ParticleMapsContainer{
  PyObject *getEnergies_numpyArray(const int particleId, double energy)
  {
      std::cerr << "ERROR: PARSEC was compiled without numpy support!" << std::endl;
      return NULL;
  }
};
%extend crpropa::ParticleMapsContainer{
  PyObject *getRandomParticles_numpyArray(const int particleId, double energy)
  {
      std::cerr << "ERROR: PARSEC was compiled without numpy support!" << std::endl;
      return NULL;
  }
};
#endif // with numpy

%pythoncode %{
class ParticleMapsContainer( ParticleMapsContainer ):
    getParticleMap = ParticleMapsContainer.getMap_numpyArray
    getParticleIds = ParticleMapsContainer.getParticleIds_numpyArray
    getEnergies = ParticleMapsContainer.getEnergies_numpyArray
    getRandomParticles = ParticleMapsContainer.getRandomParticles_numpyArray
%}

#endif // WITH_GALACTIC_LENSES_

