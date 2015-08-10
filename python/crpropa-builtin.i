/* CRPropa3 SWIG interface (for Python) */

/* Content:
 *
 * 1. SWIG settings and workarounds
 * 2. SWIG and CRPropa headers
 * 3. Pretty print for Python
 * 4. Magnetic Lens and Particle Maps Container
 *
 */


/* 1. SWIG settings and workarounds */
%include "1_swig.i"

/* 2: SWIG and CRPropa headers */
%include "2_headers.i"

%include "crpropa/Version.h"
%pythoncode %{
    __version__ = g_GIT_DESC 
%}

/* 3. Pretty print for Python */

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
            sprintf( buffer, "Vector( %.3g, %.3g, %.3g )", $self->x, $self->y, $self->z );
            return buffer;
        }
}

%template(Vector3d) crpropa::Vector3<double>;
%template(Vector3f) crpropa::Vector3<float>;

%enddef

%include "3_repr.i"

/* 4. Magnetic Lens */
%include "4_lens.i"

