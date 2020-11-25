#ifndef CRPROPA_CMZFIELD_H
#define CRPROPA_CMZFIELD_H

#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/Grid.h"

namespace crpropa {

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
    double scale(double x, double d) const;
    double Br(double r, double z, double B1, double a, double L) const;
    double Bz(double r, double z, double B1, double a, double L) const;
    double BrAz(double r, double phi, double z, double m, double B1, double eta, double R,double rr) const;
    double BphiAz(double r, double phi, double z, double m, double B1, double eta, double R, double rr) const;
    double BxAz(double x, double y, double z, double m, double B1, double eta, double R, double rr) const;
    double ByAz(double x, double y, double z, double m, double B1, double eta, double R, double rr) const;
    double BzAz(double x,double y, double z, double m, double B1, double eta, double R, double rr) const;
    double ByPol(double x,double y, double z, double B1, double a1, double a2, double eta) const;
    double BxPol(double x,double y, double z, double B1, double a1, double a2, double eta) const;
    double BzPol(double x, double y, double z, double B1, double a1, double a2, double eta) const;
    

public:

	CMZField();

    bool getUseMCField() const;
    bool getUseICField() const;
    bool getUseNTFField() const;
    bool getUseRadioArc() const;

    void setUseMCField(bool use);
    void setUseICField(bool use);
    void setUseNTFField(bool use);
    void setUseRaidoArc(bool use);

    Vector3d getMCField(const Vector3d& pos) const;
    Vector3d getICField(const Vector3d& pos) const;
    Vector3d getNTFField(const Vector3d& pos) const;
    Vector3d getRadioArcField(const Vector3d& pos) const;

    Vector3d getField(const Vector3d& pos) const;
};

} // namespace crpropa

#endif // CRPROPA_CMZFIELD_H
