#ifndef CRPROPA_CMZFIELD_H
#define CRPROPA_CMZFIELD_H

#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/Grid.h"
#include "kiss/logger.h"
#include "crpropa/Common.h"

namespace crpropa {

/**
 * \addtogroup MagneticFields
 * @{
 */


/**
 * @class CMZField
 * @brief Magnetic Field Model in the Galactic Center from M. Guenduez et al.
 * 
 * Poloidal Model (Model C) is taken from Katia Ferriere and Philippe Terral 2014 "Analytical models of X-shape magnetic fields in galactic halos" [arXiv:1312.1974]
 * Azimuthal Model is taken from M.Guenduez, J.B. Tjus, K.Ferriere, R.-J. Dettmar (2019) [arXiv:1906.05211]
 */

class CMZField: public MagneticField {
protected:
    bool useMCField;
    bool useICField;
    bool useNTFField;
    bool useRadioArc;

    /** Transform observational parameter a1 in model parameter a. Used for the poloidal model*/
    double getA(double a1) const;
    /** Transform observational parameter a2 in model parameter L. Used for the poloidal model*/
    double getL(double a2) const;
    
public:
    CMZField();

    bool getUseMCField() const;
    bool getUseICField() const;
    bool getUseNTFField() const;
    bool getUseRadioArc() const;

    void setUseMCField(bool use);
    void setUseICField(bool use);
    void setUseNTFField(bool use);
    void setUseRadioArc(bool use);

    /** Magnetic field in the poloidal model. Used for inter-cloud (IC), non-thermal-filaments (NTF) and for the radio arc.
    @param position position in galactic coordinates with Earth at (-8.5kpc, 0,0)
    @param mid  midpoint of the object
    @param B1   normalized magnetic field strength
    @param a    fitting parameter for the radial scale height
    @param L    fitting parameter for the z scale height
    @return     magnetic field vector */
    Vector3d BPol(const Vector3d& position, const Vector3d& mid, double B1, double a, double L) const;
    
    /** Magnetic field in the azimuthal model. Used for molecular clouds (MC)
    @param position position in galactic coordinates with Earth at (-8.5kpc, 0,0)
    @param mid  midpoint of the object
    @param B1   normalized magnetic field strength
    @param eta  ratio between radial and azimuthal component
    @param R    Radius of the cloud
    @return     magnetic field vector*/
    Vector3d BAz(const Vector3d& position, const Vector3d& mid, double B1, double eta, double R) const;
     
    Vector3d getMCField(const Vector3d& pos) const;
    Vector3d getICField(const Vector3d& pos) const;
    Vector3d getNTFField(const Vector3d& pos) const;
    Vector3d getRadioArcField(const Vector3d& pos) const;

    Vector3d getField(const Vector3d& pos) const;
};
/** @} */

} // namespace crpropa

#endif // CRPROPA_CMZFIELD_H
