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
%feature("python:slot", "tp_str", functype="reprfunc") classname::repr();
%feature("python:slot", "tp_repr", functype="reprfunc") classname::repr();

%extend classname {
        const std::string repr() {
            return $self->getDescription();
        }
}
%enddef

%template(Vector3d) crpropa::Vector3<double>;
%template(Vector3f) crpropa::Vector3<float>;


/* Division of vector fix #34 */
%feature("python:slot", "nb_divide", functype="binaryfunc") *::operator/;

%include "3_repr.i"

/* 4. Magnetic Lens */
%include "4_lens.i"

#ifdef WITH_GALACTIC_LENSES

%pythoncode %{

class Pixelization(Pixelization):
    def nPix(self, order=None):
      if order == None:
        return Pixelization_nPix(self.getOrder())
      else:
        return Pixelization_nPix(order)

class MagneticLens(MagneticLens):
    transformModelVector = MagneticLens.transformModelVector_numpyArray

class ParticleMapsContainer( ParticleMapsContainer ):
    getMap = ParticleMapsContainer.getMap_numpyArray
    getParticleIds = ParticleMapsContainer.getParticleIds_numpyArray
    getEnergies = ParticleMapsContainer.getEnergies_numpyArray
    getRandomParticles = ParticleMapsContainer.getRandomParticles_numpyArray

%}

%include "5_hepid.i"
#endif // WITH_GALACTIC_LENSES_

