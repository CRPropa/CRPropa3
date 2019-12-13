#ifndef CRPROPA_TD13FIELD_H
#define CRPROPA_TD13FIELD_H

#include <crpropa/magneticField/MagneticField.h>
#include <crpropa/Vector3.h>

#include <vector>

using namespace std;
using namespace crpropa;

/* Implementation of turbulent magnetic field from Tautz and Dosch 2013
 * "On numerical turbulence generation for test-particle simulations"
 * in Physics of Plasmas
 * doi: 10.1063/1.4789861
 * see https://ui.adsabs.harvard.edu/abs/2013PhPl...20b2302T/abstract
 */

class TD13Field : public MagneticField {
private:
	double B_0; // Magnetic field strength at reference level
	double Lmax; // Magnetic field maximum wave length
	double Lmin; // Magnetic field minimum wave length
    double q, s; // Magnetic field spectrum indexes, default for Kolmogorov spectrum
    int Nm; // number of Fourier modes
    int seed; // random generator seed, if 0 -> not initialized

    vector<double> k_n;
    vector<double> Ak_n;
    vector<double> eta_n; // = cos(theta_n)
    vector<double> sqrt_eta_n; // = sqrt(1-eta**2)
    vector<double> cos_phi_n;
    vector<double> sin_phi_n;
    vector<double> phase_n;
    vector<Vector3d> Xi_n;

public:
    /* Default values for a Kolmogorov spectrum */
    TD13Field(double B_0, double Lmin, double Lmax, double s=5./3., double q=0, int Nm=64, int seed=42);

    Vector3d getField(const Vector3d &pos) const;

    void setB0(double B);
    void setLmin(double L);
    void setLmax(double L);
    void setSpec(double ind1, double ind2);
    void setNm(int N);

    double getB0() const;
    double getLmin() const;
    double getLmax() const;
    double getLc() const;
    double get_q() const;
    double get_s() const;
    int getNm() const;
};

#endif // CRPROPA_TD13FIELD_H
