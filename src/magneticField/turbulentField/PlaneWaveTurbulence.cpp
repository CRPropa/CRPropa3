// This file contains an implementation of a vectorized cosine, which
// is based in part on the implementations in the library "SLEEF" by
// Naoki Shibata. SLEEF was used under the Boost Software License,
// Version 1.0. The original source file contained the following
// copyright notice:
//
//   //          Copyright Naoki Shibata 2010 - 2018.
//   // Distributed under the Boost Software License, Version 1.0.
//   //    (See accompanying file LICENSE.txt or copy at
//   //          http://www.boost.org/LICENSE_1_0.txt)
//
// SLEEF was used under the following license, which is not necessarily the
// license that applies to this file:
//
//         Boost Software License - Version 1.0 - August 17th, 2003
//
//         Permission is hereby granted, free of charge, to any person or
//         organization obtaining a copy of the software and accompanying
//         documentation covered by this license (the "Software") to use,
//         reproduce, display, distribute, execute, and transmit the Software,
//         and to prepare derivative works of the Software, and to permit
//         third-parties to whom the Software is furnished to do so, all subject
//         to the following:
//
//         The copyright notices in the Software and this entire statement,
//         including the above license grant, this restriction and the following
//         disclaimer, must be included in all copies of the Software, in whole
//         or in part, and all derivative works of the Software, unless such
//         copies or derivative works are solely in the form of
//         machine-executable object code generated by a source language
//         processor.
//
//         THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
//         EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//         MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND
//         NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE
//         DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER
//         LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
//         OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//         THE SOFTWARE.

#include "crpropa/magneticField/turbulentField/PlaneWaveTurbulence.h"
#include "crpropa/GridTools.h"
#include "crpropa/Random.h"
#include "crpropa/Units.h"

#include "kiss/logger.h"

#include <iostream>

#ifdef FAST_WAVES
#include <immintrin.h>
#endif

namespace crpropa {
#ifdef FAST_WAVES
// see
// https://stackoverflow.com/questions/49941645/get-sum-of-values-stored-in-m256d-with-sse-avx
double hsum_double_avx(__m256d v) {
	__m128d vlow = _mm256_castpd256_pd128(v);
	__m128d vhigh = _mm256_extractf128_pd(v, 1); // high 128
	vlow = _mm_add_pd(vlow, vhigh);              // reduce down to 128

	__m128d high64 = _mm_unpackhi_pd(vlow, vlow);
	return _mm_cvtsd_f64(_mm_add_sd(vlow, high64)); // reduce to scalar
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
#endif // defined(FAST_WAVES)

PlaneWaveTurbulence::PlaneWaveTurbulence(const TurbulenceSpectrum &spectrum,
                                         int Nm, int seed)
    : TurbulentField(spectrum), Nm(Nm) {

#ifdef FAST_WAVES
	KISS_LOG_INFO << "PlaneWaveTurbulence: Using SIMD TD13 implementation"
	              << std::endl;

	// In principle, we could dynamically dispatch to the non-SIMD version in
	// this case. However, this complicates the code, incurs runtime overhead,
	// and is unlikely to happen since SSE3 is quite well supported.
	// TODO: this is currently uncommented b/c sleef seems to fail to provide
	// the cpuid function
	// if (!check_sse()) {
	//  throw std::runtime_error("TD13Field: This code was compiled with SIMD
	//  support (SSE1-3), but it is running on a CPU that does not support these
	//  instructions. Please set USE_SIMD to OFF in CMake and recompile
	//  CRPropa.");
	//}
#endif

	if (Nm <= 1) {
		throw std::runtime_error(
		    "PlaneWaveTurbulence: Nm <= 1. Specify at least two wavemodes in "
		    "order to generate the k distribution properly.");
	}

	Random random;
	if (seed != 0)
		random.seed(seed);

	double kmax = 2 * M_PI / spectrum.getLmin();
	double kmin = 2 * M_PI / spectrum.getLmax();

	xi = std::vector<Vector3d>(Nm, Vector3d(0.));
	kappa = std::vector<Vector3d>(Nm, Vector3d(0.));
	phi = std::vector<double>(Nm, 0.);
	costheta = std::vector<double>(Nm, 0.);
	beta = std::vector<double>(Nm, 0.);
	Ak = std::vector<double>(Nm, 0.);
	k = std::vector<double>(Nm, 0.);

	double delta = log10(kmax / kmin);
	for (int i = 0; i < Nm; i++) {
		k[i] = pow(10, log10(kmin) + ((double)i) / ((double)(Nm - 1)) * delta);
	}

	// compute Ak
	double delta_k0 =
	    (k[1] - k[0]) / k[1]; // multiply this by k[i] to get delta_k[i]
	// on second thought, this is probably unnecessary since it's just a factor
	// and will get normalized out anyways.

	double Ak2_sum = 0; // sum of Ak^2 over all k
	// for this loop, the Ak array actually contains Gk*delta_k (ie
	// non-normalized Ak^2)
	for (int i = 0; i < Nm; i++) {
		double k = this->k[i];
		double Gk = spectrum.energySpectrum(k);
		Ak[i] = Gk * delta_k0 * k;
		Ak2_sum += Ak[i];

		// phi, costheta, and sintheta are for drawing vectors with
		// uniform distribution on the unit sphere.
		// This is similar to Random::randVector(): their t is our phi,
		// z is costheta, and r is sintheta. Our kappa is equivalent to
		// the return value of randVector(); however, TD13 then reuse
		// these values to generate a random vector perpendicular to kappa.
		double phi = random.randUniform(-M_PI, M_PI);
		double costheta = random.randUniform(-1., 1.);
		double sintheta = sqrt(1 - costheta * costheta);

		double alpha = random.randUniform(0, 2 * M_PI);
		double beta = random.randUniform(0, 2 * M_PI);

		Vector3d kappa =
		    Vector3d(sintheta * cos(phi), sintheta * sin(phi), costheta);
		Vector3d xi =
		    Vector3d(costheta * cos(phi) * cos(alpha) + sin(phi) * sin(alpha),
		             costheta * sin(phi) * cos(alpha) - cos(phi) * sin(alpha),
		             -sintheta * cos(alpha));

		this->xi[i] = xi;
		this->kappa[i] = kappa;
		this->phi[i] = phi;
		this->costheta[i] = costheta;
		this->beta[i] = beta;
	}

	// only in this loop are the actual Ak computed and stored
	// (this two-step process is necessary in order to normalize the values
	// properly)
	for (int i = 0; i < Nm; i++) {
		Ak[i] = sqrt(2 * Ak[i] / Ak2_sum) * spectrum.getBrms();
	}

#ifdef FAST_WAVES
	// * copy data into AVX-compatible arrays *
	// AVX requires all data to be aligned to 256 bit, or 32 bytes, which is the
	// same as 4 double precision floating point numbers. Since support for
	// alignments this big seems to be somewhat tentative in C++ allocators,
	// we're aligning them manually by allocating a normal double array, and
	// then computing the offset to the first value with the correct alignment.
	// This is a little bit of work, so instead of doing it separately for each
	// of the individual data arrays, we're doing it once for one big array that
	// all of the component arrays get packed into.
	//
	// The other thing to keep in mind is that AVX always reads in units of 256
	// bits, or 4 doubles. This means that our number of wavemodes must be
	// divisible by 4. If it isn't, we simply pad it out with zeros. Since the
	// final step of the computation of each wavemode is multiplication by the
	// amplitude, which will be set to 0, these padding wavemodes won't affect
	// the result.
	avx_Nm = ((Nm + 4 - 1) / 4) * 4; // round up to next larger multiple of 4:
	                                 // align is 256 = 4 * sizeof(double) bit
	avx_data = std::vector<double>(itotal * avx_Nm + 3, 0.);

	// get the first 256-bit aligned element
	size_t size = avx_data.size() * sizeof(double);
	void *pointer = avx_data.data();
	align_offset =
	    (double *)std::align(32, 32, pointer, size) - avx_data.data();

	// copy
	for (int i = 0; i < Nm; i++) {
		avx_data[i + align_offset + avx_Nm * iAxi0] = Ak[i] * xi[i].x;
		avx_data[i + align_offset + avx_Nm * iAxi1] = Ak[i] * xi[i].y;
		avx_data[i + align_offset + avx_Nm * iAxi2] = Ak[i] * xi[i].z;

		// the cosine implementation computes cos(pi*x), so we'll divide out the
		// pi here
		avx_data[i + align_offset + avx_Nm * ikkappa0] =
		    k[i] / M_PI * kappa[i].x;
		avx_data[i + align_offset + avx_Nm * ikkappa1] =
		    k[i] / M_PI * kappa[i].y;
		avx_data[i + align_offset + avx_Nm * ikkappa2] =
		    k[i] / M_PI * kappa[i].z;

		// we also need to divide beta by pi, since that goes into the argument
		// as well
		avx_data[i + align_offset + avx_Nm * ibeta] = beta[i] / M_PI;
	}
#endif // FAST_WAVES
}

Vector3d PlaneWaveTurbulence::getField(const Vector3d &pos) const {

#ifndef FAST_WAVES
	Vector3d B(0.);
	for (int i = 0; i < Nm; i++) {
		double z_ = pos.dot(kappa[i]);
		B += xi[i] * Ak[i] * cos(k[i] * z_ + beta[i]);
	}
	return B;

#else  // FAST_WAVES

	// Initialize accumulators
	//
	// There is one accumulator per component of the result vector.
	// Note that each accumulator contains four numbers. At the end of
	// the loop, each of these number will contain the sum of every
	// fourth wavemodes, starting at a different offset. In the end, all
	// of the accumulator's numbers are added together (using
	// hsum_double_avx), resulting in the total sum.

	__m256d acc0 = _mm256_setzero_pd();
	__m256d acc1 = _mm256_setzero_pd();
	__m256d acc2 = _mm256_setzero_pd();

	// broadcast position into AVX registers
	__m256d pos0 = _mm256_set1_pd(pos.x);
	__m256d pos1 = _mm256_set1_pd(pos.y);
	__m256d pos2 = _mm256_set1_pd(pos.z);

	for (int i = 0; i < avx_Nm; i += 4) {

		// load data from memory into AVX registers
		__m256d Axi0 =
		    _mm256_load_pd(avx_data.data() + i + align_offset + avx_Nm * iAxi0);
		__m256d Axi1 =
		    _mm256_load_pd(avx_data.data() + i + align_offset + avx_Nm * iAxi1);
		__m256d Axi2 =
		    _mm256_load_pd(avx_data.data() + i + align_offset + avx_Nm * iAxi2);

		__m256d kkappa0 = _mm256_load_pd(avx_data.data() + i + align_offset +
		                                 avx_Nm * ikkappa0);
		__m256d kkappa1 = _mm256_load_pd(avx_data.data() + i + align_offset +
		                                 avx_Nm * ikkappa1);
		__m256d kkappa2 = _mm256_load_pd(avx_data.data() + i + align_offset +
		                                 avx_Nm * ikkappa2);

		__m256d beta =
		    _mm256_load_pd(avx_data.data() + i + align_offset + avx_Nm * ibeta);

		// Do the computation

		// this is the scalar product between k*kappa and pos
		__m256d z = _mm256_add_pd(_mm256_mul_pd(pos0, kkappa0),
		                          _mm256_add_pd(_mm256_mul_pd(pos1, kkappa1),
		                                        _mm256_mul_pd(pos2, kkappa2)));

		// here, the phase is added on. this is the argument of the cosine.
		__m256d cos_arg = _mm256_add_pd(z, beta);

		// ********
		// * Computing the cosine
		// *
		// * argument reduction
		// step 1: compute round(x), and store it in q
		__m256d q = _mm256_round_pd(
		    cos_arg, (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC));

		// we now compute s, which will be the input parameter to our polynomial
		// approximation of cos(pi/2*x) between 0 and 1 the andnot_pd is just a
		// fast way of taking the absolute value
		__m256d s = _mm256_sub_pd(cos_arg, q);

		// the following based on the int extraction process described here:
		// https://stackoverflow.com/questions/41144668/how-to-efficiently-perform-double-int64-conversions-with-sse-avx/41223013
		// we assume -2^51 <= q < 2^51 for this, which is unproblematic, as
		// double precision has decayed far enough at that point that the cosine
		// would be useless anyway.

		// we now want to check whether q is even or odd, because the cosine is
		// negative for odd qs, so we'll have to flip the final result. on an
		// int, this is as simple as checking the 0th bit.
		// => manipulate the double in such a way that we can do this.
		// so, we add 2^52, such that the last digit of the mantissa is actually
		// in the ones position. since q may be negative, we'll also add 2^51 to
		// make sure it's positive. note that 2^51 is even and thus leaves
		// evenness invariant, which is the only thing we care about here.

		q = _mm256_add_pd(q, _mm256_set1_pd(0x0018000000000000));

		// unfortunately, integer comparisons were only introduced in avx2, so
		// we'll have to make do with a floating point comparison to check
		// whether the last bit is set. however, masking out all but the last
		// bit will result in a denormal float, which may either result in
		// performance problems or just be rounded down to zero, neither of
		// which is what we want here. To fix this, we'll mask in not only bit
		// 0, but also the exponent (and sign, but that doesn't matter) of q.
		// Luckily, the exponent of q is guaranteed to have the fixed value of
		// 1075 (corresponding to 2^52) after our addition.

		__m256d invert = _mm256_and_pd(
		    q, _mm256_castsi256_pd(_mm256_set1_epi64x(0xfff0000000000001)));

		// if we did have a one in bit 0, our result will be equal to 2^52 + 1
		invert = _mm256_cmp_pd(
		    invert, _mm256_castsi256_pd(_mm256_set1_epi64x(0x4330000000000001)),
		    _CMP_EQ_OQ);

		// finally, we need to turn invert and right_invert into masks for the
		// sign bit on each final double, ie
		invert = _mm256_and_pd(invert, _mm256_set1_pd(-0.0));

		// TODO: clamp floats between 0 and 1? This would ensure that we never
		// see inf's, but maybe we want that, so that things dont just fail
		// silently...

		// * end of argument reduction
		// *******

		// ******
		// * evaluate the cosine using a polynomial approximation
		// * the coefficients for this were generated using sleefs gencoef.c
		// * These coefficients are probably far from optimal.
		// * However, they should be sufficient for this case.
		s = _mm256_mul_pd(s, s);

		__m256d u = _mm256_set1_pd(+0.2211852080653743946e+0);

		u = _mm256_add_pd(_mm256_mul_pd(u, s),
		                  _mm256_set1_pd(-0.1332560668688523853e+1));
		u = _mm256_add_pd(_mm256_mul_pd(u, s),
		                  _mm256_set1_pd(+0.4058509506474178075e+1));
		u = _mm256_add_pd(_mm256_mul_pd(u, s),
		                  _mm256_set1_pd(-0.4934797516664651162e+1));
		u = _mm256_add_pd(_mm256_mul_pd(u, s), _mm256_set1_pd(1.));

		// then, flip the sign of each double for which invert is not zero.
		// since invert has only zero bits except for a possible one in bit 63,
		// we can xor it onto our result to selectively invert the 63st (sign)
		// bit in each double where invert is set.
		u = _mm256_xor_pd(u, invert);

		// * end computation of cosine
		// **********

		// Finally, Ak*xi is multiplied on. Since this is a vector, the
		// multiplication needs to be done for each of the three
		// components, so it happens separately.
		acc0 = _mm256_add_pd(_mm256_mul_pd(u, Axi0), acc0);
		acc1 = _mm256_add_pd(_mm256_mul_pd(u, Axi1), acc1);
		acc2 = _mm256_add_pd(_mm256_mul_pd(u, Axi2), acc2);
	}

	return Vector3d(hsum_double_avx(acc0), hsum_double_avx(acc1),
	                hsum_double_avx(acc2));
#endif // FAST_WAVES
}

} // namespace crpropa
