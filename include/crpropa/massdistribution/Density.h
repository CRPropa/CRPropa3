#ifndef CRPROPA_DENSITY_H
#define CRPROPA_DENSITY_H

#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Referenced.h"

namespace crpropa {
/**
 @class Density
 @brief Abstract base class for targetdensity
 */
class Density: public Referenced {
 public:
  virtual ~Density() {
  }
  virtual double getDensity(const Vector3d &position) const { 	// sum of all densities
    return 0;
  };
  virtual double getHIDensity(const Vector3d &position) const {
  	return 0;
  };
  virtual double getHIIDensity(const Vector3d &position) const {
  	return 0;
  };
  virtual double getH2Density(const Vector3d &position) const {		
  	return 0;
  };
  virtual double getNucleonDensity(const Vector3d &position) const{	// sum of nucleons (H2 with factor 2)
  	return 0;
  };
  
  virtual bool getisforHI() {
  	return false;
  }
  virtual bool getisforHII() {
  	return false;
  }
  virtual bool getisforH2() {
  	return false;
  }
};


}//namespace crpropa

#endif //CRPROPA_DENSITY_H
