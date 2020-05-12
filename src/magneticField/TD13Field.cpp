#include "crpropa/magneticField/TD13Field.h"

#include <crpropa/Random.h>
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace crpropa;

TD13Field::TD13Field(double B_0, double Lmin, double Lmax, double s, double q, int Nm, int seed) {
    Nmodes = Nm;
    spec_s = s;
    spec_q = q;
    spec_Lmin = Lmin;
    spec_Lmax = Lmax;

    Random random;
	if (seed != 0)
		random.seed(seed); // use given seed

    double kmin = 1./Lmax; // MF max wave length
    double kmax = 1./Lmin; // MF min wave length
    double alpha = pow(kmax/kmin , 1./(Nm-1.)); //  for log distrib norm
    double sum_Gkn_delta_kn = 0;
    vector<double> Ak_n;
    
    for(int n=0; n<Nm; n++) {
        double phi_n = random.rand(2*M_PI); 
        cos_phi_n.push_back( cos(phi_n) );
        sin_phi_n.push_back( sin(phi_n));

        double eta = random.randUniform(-1.,1.); // cos(theta_n) draw uniformaly between [-1;1]
        double sqrt_eta2 = sqrt(1-eta*eta); // = sqrt(1-eta_n**2)
        eta_n.push_back(eta); 
        sqrt_eta_n.push_back( sqrt_eta2 );
        phase_n.push_back( random.rand(2*M_PI) ); 

        // compute the random wave vector
        double alpha_n = random.rand(2*M_PI);
        Vector3d randVector;
        randVector.x = eta * cos(phi_n) * sin(alpha_n) - sin(phi_n) * cos(alpha_n);
        randVector.y = eta * sin(phi_n) * sin(alpha_n) + cos(phi_n) * cos(alpha_n);
        randVector.z = -sqrt_eta2 * sin(alpha_n);
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
        // Xi_n = Xi_n * Ak_n
        Xi_n[n].x = Xi_n[n].x * Ak_n[n]; 
        Xi_n[n].y = Xi_n[n].y * Ak_n[n]; 
        Xi_n[n].z = Xi_n[n].z * Ak_n[n]; 
    }

    return;
}

Vector3d TD13Field::getField(const Vector3d &pos) const {
    Vector3d EGMFdir;
    double z_prime, cos_kn_z;
    for(int n=0; n<Nmodes; n++) { 
        z_prime  = sqrt_eta_n[n] * cos_phi_n[n] * pos.x;
        z_prime += sqrt_eta_n[n] * sin_phi_n[n] * pos.y;
        z_prime += eta_n[n] * pos.z;
        cos_kn_z = cos( k_n[n] * z_prime + phase_n[n] );

        EGMFdir.x += cos_kn_z * Xi_n[n].x;
        EGMFdir.y += cos_kn_z * Xi_n[n].y;
        EGMFdir.z += cos_kn_z * Xi_n[n].z;
    }
    return EGMFdir;
}

double TD13Field::getLc() const {
    // According to Harari et Al JHEP03(2002)045
    double Lc;
    Lc = spec_Lmax/2.;
    Lc*= (spec_s-1.)/spec_s;
    Lc*= 1 - pow(spec_Lmin/spec_Lmax, spec_s);
    Lc/= 1 - pow(spec_Lmin/spec_Lmax, spec_s-1);
    return Lc;
}
=======
#include "crpropa/Units.h"
#include "crpropa/GridTools.h"
#include "crpropa/Random.h"

#include <iostream>

namespace crpropa {

  TD13Field::TD13Field(double kmin, double kmax, double gamma, double Nm) {

    Random random;

    // initialize everything
    this->gamma = gamma;
    this->Nm = Nm;

    xi = std::vector<Vector3d>(Nm, Vector3d(0.));
    kappa = std::vector<Vector3d>(Nm, Vector3d(0.));
    phi = std::vector<double>(Nm, 0.);
    costheta = std::vector<double>(Nm, 0.);
    beta = std::vector<double>(Nm, 0.);
    Ak = std::vector<double>(Nm, 0.);

    k = logspace(log10(kmin), log10(kmax), Nm);

    // compute Ak
    double q = 0; // TODO: what is q
    double s = gamma;
    double delta_k0 = (k[1] - k[0]) / k[1]; // multiply this by k[i] to get delta_k[i]

    double Ak2_sum = 0; // sum of Ak^2 over all k
    for (int i=0; i<Nm; i++) {
      double Gk = pow(k[i], q) / pow(1 + k[i]*k[i], (s+q)/2);
      Ak[i] = Gk * delta_k0 * k[i];
      Ak2_sum += Ak[i];
    }
    for (int i=0; i<Nm; i++) {
      Ak[i] = sqrt(Ak[i] / Ak2_sum);
    }

    // generate direction, phase, and polarizarion for each wavemode
    for (int i=0; i<Nm; i++) {
      double phi = random.randUniform(-M_PI, M_PI);
      double costheta = random.randUniform(-1., 1.);
      double sintheta = sqrt(1 - costheta*costheta);

      double alpha = random.randUniform(0, 2*M_PI);
      double beta = random.randUniform(0, 2*M_PI);

      Vector3d kappa = Vector3d ( sintheta * cos(phi), sintheta*sin(phi), costheta );
      Vector3d xi = Vector3d ( costheta*cos(phi)*cos(alpha) + sin(phi)*sin(alpha),
		     costheta*sin(phi)*cos(alpha) - cos(phi)*sin(alpha),
		     -sintheta*cos(alpha) );

      this->xi[i] = xi;
      this->kappa[i] = kappa;
      this->phi[i] = phi;
      this->costheta[i] = costheta;
      this->beta[i] = beta;
    }
}

Vector3d TD13Field::getField(const Vector3d& pos) const {
  Vector3d B(0.);

  for (int i=0; i<Nm; i++) {
    double z_ = pos.dot(kappa[i]);
    B += xi[i] * Ak[i] * cos(k[i] * z_ + beta[i]);
  }
  return B;
}

} // namespace crpropa
