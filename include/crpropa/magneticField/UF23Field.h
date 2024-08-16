#ifndef _UF23Field_h_
#define _UF23Field_h_

#include <vector>
#include "crpropa/magneticField/MagneticField.h"

namespace crpropa {
/**
 * \addtogroup MagneticFields
 * @{
 */


/**
 @class UF23Field
 @brief UF23Field Galactic magnetic field model

 Implements the eight coherent magnetic field models of UF23
 See: M. Unger and G.R. Farrar, Astrophys.J. 970 (2024) 95, arXiv:2311.12120

 Assumes a galactocentric coordinate system with the Galactic center
 at the origin, the x-axis pointing in the opposite direction of the
 Sun, and the z-axis pointing towards Galactic North.

 */

class UF23Field : public MagneticField  {
public:
  /// model variations (see Tab.2 of UF23 paper)
  enum ModelType {
    base,
    neCL,
    expX,
    spur,
    cre10,
    synCG,
    twistX,
    nebCor
  };


public:
  /**
     @brief constructor
     @param mt model type  (see Tab.2 of UF23 paper)
     @param maxRadiusInKpc maximum radius of field in kpc
  */
  UF23Field(const ModelType mt);
  /// no default constructor
  UF23Field() = delete;

  Vector3d getField(const Vector3d& pos) const;

private:

  /**
     @brief calculate coherent magnetic field at a given position
     @param posInKpc position with components given in kpc
     @return coherent field in microgauss
  */
  Vector3d operator()(const Vector3d& posInKpc) const;

  /// model parameters, see Table 3 of UF23 paper
  enum EPar {
    eDiskB1 = 0,
    eDiskB2,
    eDiskB3,
    eDiskH,
    eDiskPhase1,
    eDiskPhase2,
    eDiskPhase3,
    eDiskPitch,
    eDiskW,
    ePoloidalA,
    ePoloidalB,
    ePoloidalP,
    ePoloidalR,
    ePoloidalW,
    ePoloidalZ,
    ePoloidalXi,
    eSpurCenter,
    eSpurLength,
    eSpurWidth,
    eStriation,
    eToroidalBN,
    eToroidalBS,
    eToroidalR,
    eToroidalW,
    eToroidalZ,
    eTwistingTime,
    eNpar
  };

  /// model type given in constructor
  const ModelType fModelType;
  /// maximum galacto-centric radius beyond which B=0
  const double fMaxRadiusSquared;

  // parameters are stored in array
  double fParameters[eNpar] = { 0 };
  // references to parameters for convience
  double& fDiskB1       = fParameters[eDiskB1];
  double& fDiskB2       = fParameters[eDiskB2];
  double& fDiskB3       = fParameters[eDiskB3];
  double& fDiskH        = fParameters[eDiskH];
  double& fDiskPhase1   = fParameters[eDiskPhase1];
  double& fDiskPhase2   = fParameters[eDiskPhase2];
  double& fDiskPhase3   = fParameters[eDiskPhase3];
  double& fDiskPitch    = fParameters[eDiskPitch];
  double& fDiskW        = fParameters[eDiskW];
  double& fPoloidalA    = fParameters[ePoloidalA];
  double& fPoloidalB    = fParameters[ePoloidalB];
  double& fPoloidalP    = fParameters[ePoloidalP];
  double& fPoloidalR    = fParameters[ePoloidalR];
  double& fPoloidalW    = fParameters[ePoloidalW];
  double& fPoloidalZ    = fParameters[ePoloidalZ];
  double& fPoloidalXi   = fParameters[ePoloidalXi];
  double& fSpurCenter   = fParameters[eSpurCenter];
  double& fSpurLength   = fParameters[eSpurLength];
  double& fSpurWidth    = fParameters[eSpurWidth];
  double& fStriation    = fParameters[eStriation];
  double& fToroidalBN   = fParameters[eToroidalBN];
  double& fToroidalBS   = fParameters[eToroidalBS];
  double& fToroidalR    = fParameters[eToroidalR];
  double& fToroidalW    = fParameters[eToroidalW];
  double& fToroidalZ    = fParameters[eToroidalZ];
  double& fTwistingTime = fParameters[eTwistingTime];

  // some pre-calculated derived parameter values
  double fSinPitch  = 0;
  double fCosPitch  = 0;
  double fTanPitch  = 0;

  /// major field components
  Vector3d getDiskField(const Vector3d& pos) const;
  Vector3d getHaloField(const Vector3d& pos) const;

  /// sub-components depending on model type
  /// -- Sec. 5.2.2
  Vector3d getSpiralField(const double x, const double y, const double z) const;
  /// -- Sec. 5.2.3
  Vector3d getSpurField(const double x, const double y, const double z) const;
  /// -- Sec. 5.3.1
  Vector3d getToroidalHaloField(const double x, const double y, const double z) const;
  /// -- Sec. 5.3.2
  Vector3d getPoloidalHaloField(const double x, const double y, const double z) const;
  /// -- Sec. 5.3.3
  Vector3d getTwistedHaloField(const double x, const double y, const double z) const;

};
/** @} */
} // namespace crpropa
#endif
