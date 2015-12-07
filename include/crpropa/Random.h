// Random.h
// Mersenne Twister random number generator -- a C++ class Random
// Based on code by Makoto Matsumoto, Takuji Nishimura, and Shawn Cokus
// Richard J. Wagner  v1.0  15 May 2003  rjwagner@writeme.com

// The Mersenne Twister is an algorithm for generating random numbers.  It
// was designed with consideration of the flaws in various other generators.
// The period, 2^19937-1, and the order of equidistribution, 623 dimensions,
// are far greater.  The generator is also fast; it avoids multiplication and
// division, and it benefits from caches and pipelines.  For more information
// see the inventors' web page at http://www.math.keio.ac.jp/~matumoto/emt.html

// Reference
// M. Matsumoto and T. Nishimura, "Mersenne Twister: A 623-Dimensionally
// Equidistributed Uniform Pseudo-Random Number Generator", ACM Transactions on
// Modeling and Computer Simulation, Vol. 8, No. 1, January 1998, pp 3-30.

// Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
// Copyright (C) 2000 - 2003, Richard J. Wagner
// All rights reserved.                          
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
//   1. Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//
//   2. Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//   3. The names of its contributors may not be used to endorse or promote 
//      products derived from this software without specific prior written 
//      permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// The original code included the following notice:
//
//     When you use this, send an email to: matumoto@math.keio.ac.jp
//     with an appropriate reference to your work.
//
// It would be nice to CC: rjwagner@writeme.com and Cokus@math.washington.edu
// when you write.

// Parts of this file are modified beginning in 29.10.09 for adaption in PXL.
// Parts of this file are modified beginning in 10.02.12 for adaption in CRPropa.

#ifndef RANDOM_H
#define RANDOM_H

// Not thread safe (unless auto-initialization is avoided and each thread has
// its own Random object)
#include "crpropa/Vector3.h"

#include <iostream>
#include <limits>
#include <stdio.h>
#include <time.h>
#include <cmath>
#include <stdexcept>
#include <algorithm>

#include <stdint.h>

//necessary for win32
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace crpropa {

/**
 @class Random
 @brief Random number generator.
 */
class Random {
public:
#ifndef uint32
	typedef unsigned long uint32; // unsigned integer type, at least 32 bits
#endif
	enum {N = 624}; // length of state vector
	enum {SAVE = N + 1}; // length of array for save()

protected:
	enum {M = 397}; // period parameter
	uint32 state[N];// internal state
	uint32 *pNext;// next value to get from state
	int left;// number of values left before reload needed

//Methods
public:
	/// initialize with a simple uint32
	Random( const uint32& oneSeed );
	// initialize with an array
	Random( uint32 *const bigSeed, uint32 const seedLength = N );
	/// auto-initialize with /dev/urandom or time() and clock()
	/// Do NOT use for CRYPTOGRAPHY without securely hashing several returned
	/// values together, otherwise the generator state can be learned after
	/// reading 624 consecutive values.
	Random();
	// Access to 32-bit random numbers
	double rand();///< real number in [0,1]
	double rand( const double& n );///< real number in [0,n]
	double randExc();///< real number in [0,1)
	double randExc( const double& n );///< real number in [0,n)
	double randDblExc();///< real number in (0,1)
	double randDblExc( const double& n );///< real number in (0,n)
	/// Pull a 32-bit integer from the generator state
	/// Every other access function simply transforms the numbers extracted here
	uint32 randInt();///< integer in [0,2^32-1]
	uint32 randInt( const uint32& n );///< integer in [0,n] for n < 2^32

	uint64_t randInt64(); ///< integer in [0, 2**64 -1]. PROBABLY NOT SECURE TO USE
	uint64_t randInt64(const uint64_t &n); ///< integer in [0, n] for n < 2**64 -1. PROBABLY NOT SECURE TO USE
 
	double operator()() {return rand();} ///< same as rand()

	/// Access to 53-bit random numbers (capacity of IEEE double precision)
	double rand53();// real number in [0,1)
	///Exponential distribution in (0,inf)
	double randExponential();
	/// Normal distributed random number
	double randNorm( const double& mean = 0.0, const double& variance = 1.0 );
	/// Uniform distribution in [min, max]
	double randUniform(double min, double max);
	/// Rayleigh distributed random number	
	double randRayleigh(double sigma);
	/// Fisher distributed random number
	double randFisher(double k);

	/// Draw a random bin from a (unnormalized) cumulative distribution function, without leading zero. 
	size_t randBin(const std::vector<float> &cdf);
	size_t randBin(const std::vector<double> &cdf);

	/// Random point on a unit-sphere
	Vector3d randVector();
	/// Random vector with given angular separation around mean direction
	Vector3d randVectorAroundMean(const Vector3d &meanDirection, double angle);
	/// Fisher distributed random vector
	Vector3d randFisherVector(const Vector3d &meanDirection, double kappa);
	/// Uniform distributed random vector inside a cone
	Vector3d randConeVector(const Vector3d &meanDirection, double angularRadius);

	/// Power-law distribution of a given differential spectral index
	double randPowerLaw(double index, double min, double max);
	/// Broken power-law distribution 
	double randBrokenPowerLaw(double index1, double index2, double breakpoint, double min, double max );

	/// Seed the generator with a simple uint32
	void seed( const uint32 oneSeed );
	/// Seed the generator with an array of uint32's
	/// There are 2^19937-1 possible initial states.  This function allows
	/// all of those to be accessed by providing at least 19937 bits (with a
	/// default seed length of N = 624 uint32's).  Any bits above the lower 32
	/// in each element are discarded.
	/// Just call seed() if you want to get array from /dev/urandom
	void seed( uint32 *const bigSeed, const uint32 seedLength = N );
	/// Seed the generator with an array from /dev/urandom if available
	/// Otherwise use a hash of time() and clock() values
	void seed();

	// Saving and loading generator state
	void save( uint32* saveArray ) const;// to array of size SAVE
	void load( uint32 *const loadArray );// from such array
	friend std::ostream& operator<<( std::ostream& os, const Random& mtrand );
	friend std::istream& operator>>( std::istream& is, Random& mtrand );

	static Random &instance();
	static void seedThreads(const uint32 oneSeed);
	
protected:
	/// Initialize generator state with seed
	/// See Knuth TAOCP Vol 2, 3rd Ed, p.106 for multiplier.
	/// In previous versions, most significant bits (MSBs) of the seed affect
	/// only MSBs of the state array.  Modified 9 Jan 2002 by Makoto Matsumoto.
	void initialize( const uint32 oneSeed );

	/// Generate N new values in state
	/// Made clearer and faster by Matthew Bellew (matthew.bellew@home.com)
	void reload();
	uint32 hiBit( const uint32& u ) const {return u & 0x80000000UL;}
	uint32 loBit( const uint32& u ) const {return u & 0x00000001UL;}
	uint32 loBits( const uint32& u ) const {return u & 0x7fffffffUL;}
	uint32 mixBits( const uint32& u, const uint32& v ) const
	{	return hiBit(u) | loBits(v);}

#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4146 )
#endif
	uint32 twist( const uint32& m, const uint32& s0, const uint32& s1 ) const
	{	return m ^ (mixBits(s0,s1)>>1) ^ (-loBit(s1) & 0x9908b0dfUL);}

#ifdef _MSC_VER
#pragma warning( pop )
#endif

	/// Get a uint32 from t and c
	/// Better than uint32(x) in case x is floating point in [0,1]
	/// Based on code by Lawrence Kirby (fred@genesis.demon.co.uk)
	static uint32 hash( time_t t, clock_t c );

};

} //namespace crpropa

#endif  // RANDOM_H
