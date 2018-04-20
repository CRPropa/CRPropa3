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
    return "Vector(%.6G, %.6G, %.6G)" % (self.x, self.y, self.z)
Vector3d.__repr__ = Vector3__repr__
Vector3f.__repr__ = Vector3__repr__

%}

%enddef

%include "3_repr.i"

/* 4. Magnetic Lens */
%include "4_lens.i"

#ifdef WITH_GALACTIC_LENSES

%ignore Pixelization::nPix();

%pythoncode %{

def Pixelization_nonStaticnPix(self, order=None):
  if order == None:
    return Pixelization_nPix(self.getOrder())
  else:
    return Pixelization_nPix(order)
Pixelization.nPix = Pixelization_nonStaticnPix

MagneticLens.transformModelVector = MagneticLens.transformModelVector_numpyArray

ParticleMapsContainer.getMap = ParticleMapsContainer.getMap_numpyArray
ParticleMapsContainer.getParticleIds = ParticleMapsContainer.getParticleIds_numpyArray
ParticleMapsContainer.getEnergies = ParticleMapsContainer.getEnergies_numpyArray
ParticleMapsContainer.getRandomParticles = ParticleMapsContainer.getRandomParticles_numpyArray

%}

#endif // WITH_GALACTIC_LENSES_

