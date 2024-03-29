/* 1. SWIG settings and workarounds */

%module(directors="1", threads="1", allprotected="1") crpropa

%feature("director:except") {
  if( $error != NULL ) {
    PyObject *ptype, *pvalue, *ptraceback;
    PyErr_Fetch( &ptype, &pvalue, &ptraceback );
    PyErr_Restore( ptype, pvalue, ptraceback );
    PyErr_Print();
    Py_Exit(1);
  }
}


#ifdef WITH_DOXYGEN
  %include "docstrings_from_doxy.i"
#else
  %feature("autodoc", "1"); // automatic docstrings
#endif

%{
// workaround for SWIG < 2.0.5 with GCC >= 4.7
#include <cstddef>
  using std::ptrdiff_t;
%}

/* SWIG headers */

%include "stl.i"
%include "std_set.i"
%include "std_multiset.i"
%include "std_map.i"
%include "std_pair.i"
%include "std_multimap.i"
%include "std_vector.i"
%include "std_array.i"
%include "std_string.i"
%include "std_list.i"
%include "stdint.i"
%include "std_container.i"
%include "exception.i"
%include "std_iostream.i"
%include "attribute.i"

/* SWIG Exceptions */

%inline %{
  class RangeError {};
  class StopIterator {};
%}

%exception {
  try {
    $action
  } catch (Swig::DirectorException &e) {
    SWIG_exception(SWIG_RuntimeError, e.getMessage());
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  } catch (const char *e) {
    SWIG_exception(SWIG_RuntimeError, e);
  }
}

/* Exceptions for Python lists and iterators */
%exception __next__ {
  try {
    $action
  } catch (StopIterator) {
    PyErr_SetString(PyExc_StopIteration, "End of iterator");
    return NULL;
  }
}

%exception __getitem__ {
  try {
    $action
  } catch (RangeError) {
    SWIG_exception(SWIG_IndexError, "Index out of bounds");
    return NULL;
  }
};


/* Include numpy array interface, if available */
%{
  #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
  #include "numpy/arrayobject.h"
  #include "numpy/ufuncobject.h"
%}

%init %{
  import_array();
  import_ufunc();
%}

%pythoncode %{
  import numpy
%}


/* Hide some known warnings */
#pragma SWIG nowarn=312,325,361,389,401,503,508,509
