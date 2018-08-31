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

#include <immintrin.h>
#include <sleef.h>

#define USE_SIMD

namespace crpropa {

void print_m256d (__m256d v) {
  const double *testptr = (double *) &v;
  std::cout << testptr[0] << " " << testptr[1] << " " << testptr[2] << " " << testptr[3] << std::endl;
}

std::vector<double> logspace(double start, double stop, size_t N) {

  double delta = stop - start;
  std::vector<double> values = std::vector<double>(N, 0.);
  for (int i=0; i<N; i++) {
    values[i] = pow(10, ((double) i) / ((double) (N-1)) * delta + start);
  }
  return values;
}

//see https://stackoverflow.com/questions/49941645/get-sum-of-values-stored-in-m256d-with-sse-avx
double hsum_double_avx(__m256d v) {
  __m128d vlow  = _mm256_castpd256_pd128(v);
  __m128d vhigh = _mm256_extractf128_pd(v, 1); // high 128
  vlow  = _mm_add_pd(vlow, vhigh);     // reduce down to 128

  __m128d high64 = _mm_unpackhi_pd(vlow, vlow);
  return  _mm_cvtsd_f64(_mm_add_sd(vlow, high64));  // reduce to scalar
}

// code for hsum_float_avx taken from:
// https://stackoverflow.com/questions/13219146/how-to-sum-m256-horizontally

// x = ( x7, x6, x5, x4, x3, x2, x1, x0 )
float hsum_float_avx(__m256 x) {
  // hiQuad = ( x7, x6, x5, x4 )
  const __m128 hiQuad = _mm256_extractf128_ps(x, 1);
  // loQuad = ( x3, x2, x1, x0 )
  const __m128 loQuad = _mm256_castps256_ps128(x);
  // sumQuad = ( x3 + x7, x2 + x6, x1 + x5, x0 + x4 )
  const __m128 sumQuad = _mm_add_ps(loQuad, hiQuad);
  // loDual = ( -, -, x1 + x5, x0 + x4 )
  const __m128 loDual = sumQuad;
  // hiDual = ( -, -, x3 + x7, x2 + x6 )
  const __m128 hiDual = _mm_movehl_ps(sumQuad, sumQuad);
  // sumDual = ( -, -, x1 + x3 + x5 + x7, x0 + x2 + x4 + x6 )
  const __m128 sumDual = _mm_add_ps(loDual, hiDual);
  // lo = ( -, -, -, x0 + x2 + x4 + x6 )
  const __m128 lo = sumDual;
  // hi = ( -, -, -, x1 + x3 + x5 + x7 )
  const __m128 hi = _mm_shuffle_ps(sumDual, sumDual, 0x1);
  // sum = ( -, -, -, x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7 )
  const __m128 sum = _mm_add_ss(lo, hi);
  return _mm_cvtss_f32(sum);
}

  TD13Field::TD13Field(double Brms, double kmin, double kmax, double gamma, double bendoverScale, int Nm, int seed) {

    // NOTE: the use of the turbulence bend-over scale in the TD13 paper is quite confusing to
    // me. The paper states that k = l_0 * <k tilde> would be used throughout, yet
    // surely they do not mean to say that l_0 * <k tilde> should be used for the k in the
    // scalar product in eq. 2? In this implementation, I've only multiplied in the l_0
    // in the computation of the Gk, not the actual <k>s used for planar wave evaluation,
    // since this would yield obviously wrong results...

    if (kmin > kmax) {
      throw std::runtime_error("TD13Field: kmin > kmax");
    }

    if (Nm <= 1) {
      throw std::runtime_error("TD13Field: Nm <= 1. We need at least two wavemodes in order to generate the k distribution properly, and besides -- *what are you doing?!*");
    }

    if (kmin <= 0) {
      throw std::runtime_error("TD13Field: kmin <= 0");
    } 

    Random random;
    if (seed != 0) { // copied from initTurbulence
      random.seed(seed);
    }

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
    //on second thought, this is probably unnecessary since it's just a factor and will get
    //normalized out anyways.

    double Ak2_sum = 0; // sum of Ak^2 over all k
    //for this loop, the Ak array actually contains Gk*delta_k (ie non-normalized Ak^2)
    for (int i=0; i<Nm; i++) {
      double k = this->k[i] * bendoverScale;
      double Gk = pow(k, q) / pow(1 + k*k, (s+q)/2);
      Ak[i] = Gk * delta_k0 * k  *k*k; //DEBUG volume correction factor
      Ak2_sum += Ak[i];
    }
    //only in this loop are the actual Ak computed and stored
    //(this two-step process is necessary in order to normalize the values properly)
    for (int i=0; i<Nm; i++) {
      Ak[i] = sqrt(Ak[i] / Ak2_sum * 2) * Brms;
    }

    // generate direction, phase, and polarization for each wavemode
    for (int i=0; i<Nm; i++) {
      double phi = random.randUniform(-M_PI, M_PI);
      double costheta = random.randUniform(-1., 1.);
      //// DEBUG set these to zero for aligned FFT
      //phi = 0.;
      //costheta = 0.;
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
    //copy data into AVX-compatible arrays
    avx_Nm = ( (Nm + 8 - 1)/8 ) * 8; //round up to next larger multiple of 8: align is 256 = 8 * sizeof(float) bit
    //std::cout << avx_Nm <<std::endl;
    //std::cout << itotal << std::endl;
    //std::cout << (itotal*avx_Nm + 7) << std::endl;
    avx_data = std::vector<float>(itotal*avx_Nm + 7, 0.);

    //get the first 256-bit aligned element
    size_t size = avx_data.size()*sizeof(float);
    void *pointer = avx_data.data();
    align_offset = (float *) std::align(32, 32, pointer, size) - avx_data.data();

    //copy
    for (int i=0; i<Nm; i++) {
      avx_data[i + align_offset + avx_Nm*ixi0] = xi[i].x;
      avx_data[i + align_offset + avx_Nm*ixi1] = xi[i].y;
      avx_data[i + align_offset + avx_Nm*ixi2] = xi[i].z;

      avx_data[i + align_offset + avx_Nm*ikappa0] = kappa[i].x;
      avx_data[i + align_offset + avx_Nm*ikappa1] = kappa[i].y;
      avx_data[i + align_offset + avx_Nm*ikappa2] = kappa[i].z;

      avx_data[i + align_offset + avx_Nm*iAk] = Ak[i];
      avx_data[i + align_offset + avx_Nm*ik] = k[i];
      avx_data[i + align_offset + avx_Nm*ibeta] = beta[i];
    }
}

Vector3d TD13Field::getField(const Vector3d& pos) const {

#ifndef USE_SIMD
  Vector3d B(0.);
  
  for (int i=0; i<Nm; i++) {
    double z_ = pos.dot(kappa[i]);
    B += xi[i] * Ak[i] * cos(k[i] * z_ + beta[i]);
  }

  return B;

#else
  __m256 acc0 = _mm256_setzero_ps();
  __m256 acc1 = _mm256_setzero_ps();
  __m256 acc2 = _mm256_setzero_ps();

  __m256 pos0 = _mm256_set1_ps(pos.x);
  __m256 pos1 = _mm256_set1_ps(pos.y);
  __m256 pos2 = _mm256_set1_ps(pos.z);

  __m256 test;

  for (int i=0; i<avx_Nm; i+=8) {

    //load data from memory into AVX registers
    __m256 xi0 = _mm256_load_ps(avx_data.data() + i + align_offset + avx_Nm*ixi0);
    __m256 xi1 = _mm256_load_ps(avx_data.data() + i + align_offset + avx_Nm*ixi1);
    __m256 xi2 = _mm256_load_ps(avx_data.data() + i + align_offset + avx_Nm*ixi2);

    __m256 kappa0 = _mm256_load_ps(avx_data.data() + i + align_offset + avx_Nm*ikappa0);
    __m256 kappa1 = _mm256_load_ps(avx_data.data() + i + align_offset + avx_Nm*ikappa1);
    __m256 kappa2 = _mm256_load_ps(avx_data.data() + i + align_offset + avx_Nm*ikappa2);

    __m256 Ak = _mm256_load_ps(avx_data.data() + i + align_offset + avx_Nm*iAk);
    __m256 k = _mm256_load_ps(avx_data.data() + i + align_offset + avx_Nm*ik);
    __m256 beta = _mm256_load_ps(avx_data.data() + i + align_offset + avx_Nm*ibeta);

    //do the computation
    __m256 z = _mm256_add_ps(_mm256_mul_ps(pos0, kappa0),
			      _mm256_add_ps(_mm256_mul_ps(pos1, kappa1),
					    _mm256_mul_ps(pos2, kappa2)
					    )
			      );

    __m256 cos_arg = _mm256_add_ps(_mm256_mul_ps(k, z), beta);
    __m256 mag = _mm256_mul_ps(Ak, Sleef_cosf8_u35(cos_arg));

    acc0 = _mm256_add_ps(_mm256_mul_ps(mag, xi0), acc0);
    acc1 = _mm256_add_ps(_mm256_mul_ps(mag, xi1), acc1);
    acc2 = _mm256_add_ps(_mm256_mul_ps(mag, xi2), acc2);
  }
  
  return Vector3d(hsum_float_avx(acc0),
                  hsum_float_avx(acc1),
                  hsum_float_avx(acc2)
                  );
#endif
}

} // namespace crpropa
