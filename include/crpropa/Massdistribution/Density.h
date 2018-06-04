#ifndef CRPROPA_DENSITY_H
#define CRPROPA_DENSITY_H

#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Referenced.h"

namespace crpropa {
/**
 @class Density
 @brief Abstract base class for Targetdensity
 */
class Density: public Referenced {
bool isforHI;
bool isforHII;
bool isforH2;
 public:
  virtual ~Density() {
  }
  virtual double getDensity(const Vector3d &position) const {
    return 0;
  };
  virtual double getDensity_HE(const Vector3d &position) const {
    return getDensity(position)*0.11;		// H-He ratio
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
	
	bool getisforHI(){
		return isforHI;
	};
	bool getisforHII(){
		return isforHII;
	};
	bool getisforH2(){	
		return isforH2;
	};

};

}//namespace crpropa

#endif //CRPROPA_DENSITY_H
