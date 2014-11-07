#ifndef CRPROPA_FITSMAGNETICFIELD_H
#define CRPROPA_FITSMAGNETICFIELD_H

#ifdef CRPROPA_HAVE_CFITSIO

#include <iostream> 
#include <string>
#include <cstdio>

#include "crpropa/Units.h"
#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/Vector3.h"

#include "fitsio.h"

namespace crpropa {
  
  /**
     @class FITSMagneticField
     @brief Read MagneticField from a given FITS file
  */
  class FITSMagneticField: public MagneticField {
    
  private:
    
    double cfLength;
    double cfDensity;
    double cfMagneticField;
    
  public:        
    
    FITSMagneticField(){};
    FITSMagneticField(const std::string&);
    Vector3d getField(const Vector3d&) const;

    //dtr
  };
  
} // namespace crpropa

#endif // CRPROPA_HAVE_CFITSIO
#endif // CRPROPA_FITSMAGNETICFIELD_H
