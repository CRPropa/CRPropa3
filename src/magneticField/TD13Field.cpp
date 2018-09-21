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
#include "kiss/logger.h"

#include <iostream>

#if defined(CRPROPA_HAVE_SLEEF) && defined(__SSE__) && defined(__SSE2__) && defined(__SSE3__)
#define FAST_TD13

#include <immintrin.h>
#include <sleef.h>
#endif

namespace crpropa {

std::vector<double> logspace(double start, double stop, size_t N) {

  double delta = stop - start;
  std::vector<double> values = std::vector<double>(N, 0.);
  for (int i=0; i<N; i++) {
    values[i] = pow(10, ((double) i) / ((double) (N-1)) * delta + start);
  }
  return values;
}

#ifdef FAST_TD13
// code from:
// https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-float-vector-sum-on-x86
float hsum_float_sse3(__m128 v) {
    __m128 shuf = _mm_movehdup_ps(v);        // broadcast elements 3,1 to 2,0
    __m128 sums = _mm_add_ps(v, shuf);
    shuf        = _mm_movehl_ps(shuf, sums); // high half -> low half
    sums        = _mm_add_ss(sums, shuf);
    return        _mm_cvtss_f32(sums);
}

  //bool check_sse() {
  //  int32_t result[4];
  //  Sleef_x86CpuID(result, 1, 0);
  //  return (result[3] & (1 << 25)) && (result[3] & (1 << 26)) && (result[2] & (1 << 0))
  //    && (result[2] & (1 << 28)) //DEBUG check for avx, to test if this is working.
	//    ;
  //}
#endif // defined(FAST_TD13)

  TD13Field::TD13Field(double Brms, double kmin, double kmax, double gamma, double bendoverScale, int Nm, int seed) {

    // NOTE: the use of the turbulence bend-over scale in the TD13 paper is quite confusing to
    // me. The paper states that k = l_0 * <k tilde> would be used throughout, yet
    // surely they do not mean to say that l_0 * <k tilde> should be used for the k in the
    // scalar product in eq. 2? In this implementation, I've only multiplied in the l_0
    // in the computation of the Gk, not the actual <k>s used for planar wave evaluation,
    // since this would yield obviously wrong results...

#ifdef FAST_TD13
    KISS_LOG_INFO << "TD13Field: Using SIMD TD13 implementation" << std::endl;

    // In principle, we could dynamically dispatch to the non-SIMD version in
    // this case. However, this complicates the code, incurs runtime overhead,
    // and is unlikely to happen since SSE3 is quite well supported.
    // TODO: this is currently uncommented b/c sleef seems to fail to provide
    // the cpuid function
    //if (!check_sse()) {
    //  throw std::runtime_error("TD13Field: This code was compiled with SIMD support (SSE1-3), but it is running on a CPU that does not support these instructions. Please set USE_SIMD to OFF in CMake and recompile CRPropa.");
    //}
#endif

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
      // phi, costheta, and sintheta are for drawing vectors with
      // uniform distribution on the unit sphere.
      // This is similar to Random::randVector(): their t is our phi,
      // z is costheta, and r is sintheta. Our kappa is equivalent to
      // the return value of randVector(); however, TD13 then reuse
      // these values to generate a random vector perpendicular to kappa.
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

    // copy data into AVX-compatible arrays
    //
    // What is going on here:
    // SIMD load instructions require data to be aligned in memory to
    // the SIMD register bit size (128 bit for SSE, 256 bit for
    // AVX). Specifying memory alignment in C++ felt somewhat icky to
    // me after doing initial research, so here we're using the simple
    // solution: Allocating a vector that is a little bit longer than
    // required, determining the first aligned index, and then only
    // using the array starting at that index. (Here called align_offset.)
    //
    // To cut down on the complexity of storing one such offset for
    // each of the attribute arrays, we're using only a single array to
    // store all of the wavemode attributes, with offset indices to
    // each. This way, only a single align_offset has to be stored.
    //
    // This code was originally written for AVX, hence the
    // naming. I've also kept the larger 256-bit alignment required
    // for AVX, even though SSE only needs 128-bit alignment, to make
    // it simpler to go back to AVX. (And because there's really no
    // disadvantage, except for allocating a few bytes more.)

    avx_Nm = ( (Nm + 8 - 1)/8 ) * 8; //round up to next larger multiple of 8: align is 256 = 8 * sizeof(float) bit
    avx_data = std::vector<float>(itotal*avx_Nm + 7, 0.);

    //get the first 256-bit aligned element
    size_t size = avx_data.size()*sizeof(float);
    void *pointer = avx_data.data();
    align_offset = (float *) std::align(32, 32, pointer, size) - avx_data.data();

    //copy
    for (int i=0; i<Nm; i++) {
      avx_data[i + align_offset + avx_Nm*iAxi0] = Ak[i] * xi[i].x;
      avx_data[i + align_offset + avx_Nm*iAxi1] = Ak[i] * xi[i].y;
      avx_data[i + align_offset + avx_Nm*iAxi2] = Ak[i] * xi[i].z;

      avx_data[i + align_offset + avx_Nm*ikkappa0] = k[i] * kappa[i].x;
      avx_data[i + align_offset + avx_Nm*ikkappa1] = k[i] * kappa[i].y;
      avx_data[i + align_offset + avx_Nm*ikkappa2] = k[i] * kappa[i].z;

      avx_data[i + align_offset + avx_Nm*ibeta] = beta[i];
    }
}

Vector3d TD13Field::getField(const Vector3d& pos) const {

#ifndef FAST_TD13
  Vector3d B(0.);
  
  for (int i=0; i<Nm; i++) {
    double z_ = pos.dot(kappa[i]);
    B += xi[i] * Ak[i] * cos(k[i] * z_ + beta[i]);
  }

  return B;

#else // CRPROPA_USE_SIMD
  __m128 acc0 = _mm_setzero_ps();
  __m128 acc1 = _mm_setzero_ps();
  __m128 acc2 = _mm_setzero_ps();

  __m128 pos0 = _mm_set1_ps(pos.x);
  __m128 pos1 = _mm_set1_ps(pos.y);
  __m128 pos2 = _mm_set1_ps(pos.z);

  __m128 test;

  for (int i=0; i<avx_Nm; i+=4) {

    //load data from memory into AVX registers
    __m128 Axi0 = _mm_load_ps(avx_data.data() + i + align_offset + avx_Nm*iAxi0);
    __m128 Axi1 = _mm_load_ps(avx_data.data() + i + align_offset + avx_Nm*iAxi1);
    __m128 Axi2 = _mm_load_ps(avx_data.data() + i + align_offset + avx_Nm*iAxi2);

    __m128 kkappa0 = _mm_load_ps(avx_data.data() + i + align_offset + avx_Nm*ikkappa0);
    __m128 kkappa1 = _mm_load_ps(avx_data.data() + i + align_offset + avx_Nm*ikkappa1);
    __m128 kkappa2 = _mm_load_ps(avx_data.data() + i + align_offset + avx_Nm*ikkappa2);

    __m128 beta = _mm_load_ps(avx_data.data() + i + align_offset + avx_Nm*ibeta);

    //do the computation
    __m128 z = _mm_add_ps(_mm_mul_ps(pos0, kkappa0),
			      _mm_add_ps(_mm_mul_ps(pos1, kkappa1),
					    _mm_mul_ps(pos2, kkappa2)
					    )
			      );

    __m128 cos_arg = _mm_add_ps(z, beta);
    __m128 mag = Sleef_cosf4_u35(cos_arg);

    acc0 = _mm_add_ps(_mm_mul_ps(mag, Axi0), acc0);
    acc1 = _mm_add_ps(_mm_mul_ps(mag, Axi1), acc1);
    acc2 = _mm_add_ps(_mm_mul_ps(mag, Axi2), acc2);
  }
  
  return Vector3d(hsum_float_sse3(acc0),
                  hsum_float_sse3(acc1),
                  hsum_float_sse3(acc2)
                  );
#endif // FAST_TD13
}

} // namespace crpropa
