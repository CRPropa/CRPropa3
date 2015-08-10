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

