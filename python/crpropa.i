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

/* 3. Pretty print for Python */

%define __REPR__( classname ) 

%pythoncode %{
globals()["classname"["classname".find('::')+2:]].__repr__ = globals()["classname"["classname".find('::')+2:]].getDescription
%}

%enddef

%define VECTOR3__REPR__( classname ) 

%template(Vector3d) classname<double>;
%template(Vector3f) classname<float>;

%pythoncode %{

def Vector3__repr__(self):
    return "Vector(%.3g, %.3g, %.3g)" % (self.x, self.y, self.z)
Vector3d.__repr__ = Vector3__repr__
Vector3f.__repr__ = Vector3__repr__

%}

%enddef


%include "3_repr.i"

/* 4. Magnetic Lens */
%include "4_lens.i"

