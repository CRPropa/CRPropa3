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
        Py_RETURN_NONE; 
      }

      if (arr == NULL) 
      {
        Py_RETURN_NONE;
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
      Py_RETURN_NONE;
    }
};
#endif



/* 5. Particle Maps Container */

%ignore ParticleMapsContainer::getMap;
%ignore ParticleMapsContainer::getParticleIds;
%ignore ParticleMapsContainer::getEnergies;
%ignore ParticleMapsContainer::getRandomParticles;
%include "crpropa/magneticLens/ParticleMapsContainer.h"

#ifdef WITHNUMPY
%extend crpropa::ParticleMapsContainer{

		PyObject *addParticles(PyObject *particleIds, 
        PyObject *energies, 
        PyObject *galacticLongitudes, 
        PyObject *galacticLatitudes, 
        PyObject *weights)
    {
      //ToDo: Handle strided arrays 
      
      //ToDo: Check that input objects are arrays  PyArray_Check
      if (!PyArray_Check(particleIds))
      {
        std::cerr << "ParticleMapsContainer::addParticles -  require array as input for particleIds\n";
        Py_RETURN_NONE;
      }
      if (!PyArray_Check(energies))
      {
        std::cerr << "ParticleMapsContainer::addParticles -  require array as input for energy\n";
        Py_RETURN_NONE;
      }
      if (!PyArray_Check(galacticLongitudes))
      {
        std::cerr << "ParticleMapsContainer::addParticles -  require array as input for galacticLongitudes\n";
        Py_RETURN_NONE;
      }
      if (!PyArray_Check(galacticLatitudes))
      {
        std::cerr << "ParticleMapsContainer::addParticles -  require array as input for galacticLatitudes\n";
        Py_RETURN_NONE;
      }
      if (!PyArray_Check(weights))
      {
        std::cerr << "ParticleMapsContainer::addParticles -  require array as input for weights\n";
        Py_RETURN_NONE;
      }
      
      int ndim = 0;
      npy_intp dims[NPY_MAXDIMS];

      PyArrayObject *particleIds_arr = NULL;
      PyArray_Descr *particleIds_dtype = NULL;
      if (PyArray_GetArrayParamsFromObject(particleIds, NULL, 1,
            &particleIds_dtype, &ndim, dims, &particleIds_arr, NULL) < 0)
        { Py_RETURN_NONE; };
      if (particleIds_arr == NULL)
        Py_RETURN_NONE;
      
      npy_intp *D = PyArray_DIMS(particleIds_arr);
      int arraySize = D[0];

      PyArrayObject *energies_arr = NULL;
      PyArray_Descr *energies_dtype = NULL;
      
      if (PyArray_GetArrayParamsFromObject(energies, NULL, 1,
            &energies_dtype, &ndim, dims, &energies_arr,
            NULL) < 0)
        { Py_RETURN_NONE; };
      if (energies_arr == NULL)
        Py_RETURN_NONE;

      PyArrayObject *galacticLongitudes_arr = NULL;
      PyArray_Descr *galacticLongitudes_dtype = NULL;
      if (PyArray_GetArrayParamsFromObject(galacticLongitudes, NULL, 1,
            &galacticLongitudes_dtype, &ndim, dims, &galacticLongitudes_arr,
            NULL) < 0)
        { Py_RETURN_NONE; };
      if (galacticLongitudes_arr == NULL)
        Py_RETURN_NONE;

      PyArrayObject *galacticLatitudes_arr = NULL;
      PyArray_Descr *galacticLatitudes_dtype = NULL;
      if (PyArray_GetArrayParamsFromObject(galacticLatitudes, NULL, 1,
            &galacticLatitudes_dtype, &ndim, dims, &galacticLatitudes_arr,
            NULL) < 0)
        { Py_RETURN_NONE; };
      if (galacticLatitudes_arr == NULL)
        Py_RETURN_NONE;

      PyArrayObject *weights_arr = NULL;
      PyArray_Descr *weights_dtype = NULL;
      if (PyArray_GetArrayParamsFromObject(weights, NULL, 1,
            &weights_dtype, &ndim, dims, &weights_arr,
            NULL) < 0)
        { Py_RETURN_NONE; };
      if (weights_arr == NULL)
        Py_RETURN_NONE;

      int *particleIds_dp = (int*) PyArray_DATA(particleIds_arr);
      double *energies_dp = (double*) PyArray_DATA(energies_arr);
      double *galacticLongitudes_dp = (double*) PyArray_DATA(galacticLongitudes_arr);
      double *galacticLatitudes_dp = (double*) PyArray_DATA(galacticLatitudes_arr);
      double *weights_dp= (double*) PyArray_DATA(weights_arr);

      for(size_t i =0; i < arraySize; i++ )
      {
        $self->addParticle(particleIds_dp[i], energies_dp[i],
            galacticLongitudes_dp[i], galacticLatitudes_dp[i], weights_dp[i]);

      }
      Py_RETURN_TRUE;
    }
      

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
      Py_RETURN_NONE;
  }
};
%extend crpropa::ParticleMapsContainer{
  PyObject *getParticleIds_numpyArray()
  {
      std::cerr << "ERROR: PARSEC was compiled without numpy support!" << std::endl;
      Py_RETURN_NONE;
  }
};
%extend crpropa::ParticleMapsContainer{
  PyObject *getEnergies_numpyArray(const int pid)
  {
      std::cerr << "ERROR: PARSEC was compiled without numpy support!" << std::endl;
      Py_RETURN_NONE;
  }
};
%extend crpropa::ParticleMapsContainer{
  PyObject *getRandomParticles_numpyArray(size_t N)
  {
      std::cerr << "ERROR: PARSEC was compiled without numpy support!" << std::endl;
      Py_RETURN_NONE;
  }
};
#endif // with numpy

#endif // WITH_GALACTIC_LENSES_

