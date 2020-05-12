#include "crpropa/magneticField/TD13Field.h"

#include <crpropa/Random.h>
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace crpropa;

TD13Field::TD13Field() {
}

void TD13Field::initTurbulence(double B_0, double Lmin, double Lmax, double s, double q, int Nm, int seed) {
    setNm(Nm);
    Random random;
	if (seed != 0)
		random.seed(seed); // use given seed

    double kmin = 1./Lmax; // MF max wave length
    double kmax = 1./Lmin; // MF min wave length
    double alpha = pow(kmax/kmin , 1./(Nm-1.)); //  for log distrib norm
    double sum_Gkn_delta_kn = 0;
    
    for(int n=0; n<Nm; n++) {
        double alpha_n = random.rand(2*M_PI);
        phi_n.push_back( random.rand(2*M_PI) ); 
        eta_n.push_back( random.randUniform(-1.,1.) ); // cos(theta_n) draw uniformaly between [-1;1]
        phase_n.push_back( random.rand(2*M_PI) ); 

        // compute the random wave vector
        Vector3d randVector;
        randVector.x = eta_n[n] * cos(phi_n[n]) * sin(alpha_n) - sin(phi_n[n]) * cos(alpha_n);
        randVector.y = eta_n[n] * sin(phi_n[n]) * sin(alpha_n) + cos(phi_n[n]) * cos(alpha_n);
        randVector.z = -sqrt(1 - eta_n[n]*eta_n[n]) * sin(alpha_n);
        Xi_n.push_back(randVector);

        k_n.push_back( pow(alpha, n) * kmin);
        double delta_kn = k_n[n]*(alpha-1);
        double G_kn = pow(k_n[n],q) / pow( 1+k_n[n]*k_n[n], (s+q)/2. );
        Ak_n.push_back( G_kn*delta_kn ); // first part eq. 4 Tautz and Dosch 2013: A^2(k) = G(k)*delta_k
        sum_Gkn_delta_kn += Ak_n[n]; // for the second part of the equation 
    }

    double epsilon = sqrt(2.); // correction factor (see Tautz and Dosch 2013, section II.B.3 )
    for(int n=0; n<Nm; n++) { // second part eq. 4 Tautz and Dosch 2013
        Ak_n[n] = epsilon * B_0 * sqrt( Ak_n[n] / sum_Gkn_delta_kn );
    }

    return;
}

Vector3d TD13Field::getField(const Vector3d &pos) const {
    Vector3d EGMFdir;
    double z_prime, cos_kn_z;
    for(int n=0; n<Nm; n++) { 
        z_prime  = sqrt(1 - eta_n[n]*eta_n[n]) * cos(phi_n[n]) * pos.x;
        z_prime += sqrt(1 - eta_n[n]*eta_n[n]) * sin(phi_n[n]) * pos.y;
        z_prime += eta_n[n] * pos.z;
        cos_kn_z = cos( k_n[n] * z_prime + phase_n[n] );

        EGMFdir.x += Ak_n[n] * cos_kn_z * Xi_n[n].x;
        EGMFdir.y += Ak_n[n] * cos_kn_z * Xi_n[n].y;
        EGMFdir.z += Ak_n[n] * cos_kn_z * Xi_n[n].z;
    }
    return EGMFdir;
}

void TD13Field::setB0(double B){
    B_0 = B;
    return;
}

void TD13Field::setLmin(double L){
    Lmin = L;
    return;
}

void TD13Field::setLmax(double L){
    Lmax = L;
    return;
}

void TD13Field::setSpec(double ind1, double ind2){
    s = ind1;
    q = ind2;
    return;
}

void TD13Field::setNm(int N){
    Nm = N;
    return;
}

double TD13Field::getB0() const {
    return B_0;
}

double TD13Field::getLmin() const {
    return Lmin;
}

double TD13Field::getLmax() const {
    return Lmax;
}

double TD13Field::getLc() const {
    // According to Harari et Al JHEP03(2002)045
    double Lc;
    Lc = Lmax/2.;
    Lc*= (s-1.)/s;
    Lc*= 1 - pow(Lmin/Lmax,s);
    Lc/= 1 - pow(Lmin/Lmax,s-1);
    return Lc;
}

double TD13Field::get_q() const {
    return q;
}

double TD13Field::get_s() const {
    return s;
}

int TD13Field::getNm() const {
    return Nm;
}
