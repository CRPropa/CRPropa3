// shealpix.h, 12.04.10 by Tobias Winchen, <tobias.winchen@rwth-aachen.de>
//
// This defines the class sHEALPIX (small HEALPIX) which provides a very
// basic Hierarchical Equal Area isoLatitude Pixelization (HEALPIX) of a
// sphere like http://healpix.jpl.nasa.gov (on which it is based).  
// sHEALPIX so far only provides for pixelizations in the RING scheme:
//	- assoziation of a dircetion (theta, phi) with the pixel number 
//	- vice versa
//
// shealpix tries to be compatible to bhealpix_base
//
// sHEALPIX is published unter the GNU GPL v3 as it is based on HEALPIX
// which is published under the GNU GPL v2.
//
// sHEALPIX is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// sHEALPIX is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with sHEALPIX; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//

#ifndef SHEALPIX_HH
#define SHEALPIX_HH

#include <cmath>
#include <stdint.h>

#include <iostream>


namespace shealpix
{

/*! Class representing a 3D cartesian vector in healpix, copied from
 * vec3.h. */
class vec3
{
public:
	double x; /*!< x-coordinate */
	double y; /*!< y-coordinate */
	double z; /*!< z-coordinate */

	/*! Default constructor. Does not initialize \a x, \a y, and \a z. */
	vec3 () {}
	/*! Creates a vector with the coordinates \a xc, \a yc, and \a zc. */
	vec3 (double xc, double yc, double zc) : x(xc), y(yc), z(zc) {}

	/*! Creates a unit vector from a z coordinate and an azimuthal angle. */
	void set_z_phi (double z_, double phi)
	{
		using namespace std;
		double sintheta = sqrt((1.-z_)*(1.+z_));
		x = sintheta*cos(phi);
		y = sintheta*sin(phi);
		z = z_;
	}

/*! Normalizes the vector to length 1. */
	void Normalize ()
	{
		using namespace std;
		double l = 1.0/sqrt (x*x + y*y + z*z);
		x*=l; y*=l; z*=l;
	}

	/*! Returns the length of the vector. */
	double Length () const
	{ return sqrt (x*x + y*y + z*z); }

	/*! Returns the squared length of the vector. */
	double SquaredLength () const
		{ return (x*x + y*y + z*z); }
	/*! Returns the vector with the signs of all coordinates flipped. */
	const vec3 operator- () const
		{ return vec3 (-x, -y, -z); }
	/*! Flips the signs of all coordinates. */
	void Flip ()
		{ x=-x; y=-y; z=-z; }
	/*! Subtracts \a vec from the vector. */
	const vec3 operator- (const vec3 &vec) const
		{ return vec3 (x-vec.x, y-vec.y, z-vec.z); }
	/*! Adds \a vec to the vector. */
	const vec3 operator+ (const vec3 &vec) const
		{ return vec3 (x+vec.x, y+vec.y, z+vec.z); }
	/*! Returns the vector scaled by \a fact. */
	const vec3 operator* (double fact) const
		{ return vec3 (x*fact, y*fact, z*fact); }
	/*! Returns the vector scaled by \a 1/fact. */
	const vec3 operator/ (double fact) const
		{ double xfact = 1./fact; return vec3 (x*xfact, y*xfact, z*xfact); }
	/*! Scales the vector by \a fact. */
	vec3 &operator*= (double fact)
		{ x*=fact; y*=fact; z*=fact; return *this; }
};

/*! Returns the dot product of \a v1 and \a v2.
    \relates vec3 */
inline double dotprod(const vec3 &v1, const vec3 &v2)
  { return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z; }

/*! Returns the cross product of \a a and \a b.
    \relates vec3 */
inline vec3 crossprod(const vec3 &a, const vec3 &b)
  { return vec3 (a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x); }

/*! Writes \a v to \a os.
    \relates vec3 */
inline std::ostream &operator<< (std::ostream &os, const vec3 &v)
  {
  os << v.x << ", " << v.y << ", " << v.z << std::endl;
  return os;
  }








class sHEALPIX
{

	public:

		sHEALPIX( ) : _order(0),  _npix(0), _ncap(0), _nside(0), _npface(0), _fact1(0), _fact2(0)
		{
		}

		sHEALPIX(int order)
		{
			_order = order;
			// copied from healpix_base.h
			_nside  = 1<<_order; // _nside = 2**order
			_npface = _nside << _order; // _nside * 2**order == nside**2
			_ncap   = (_npface-_nside)<<1; //(nside**2-nside) * 2**1
			_npix   = 12*_npface;
			_fact2  = 4./_npix;
			_fact1  = (_nside<<1)*_fact2; // 2*nside * fact2
		}

		void pix2ang_z_phi (int pix, double &z, double &phi) const
		{
			if (pix<_ncap) // North Polar cap
			{
				int iring = int(0.5*(1+isqrt(1+2*pix))); //counted from North pole
				int iphi  = (pix+1) - 2*iring*(iring-1);

				z = 1.0 - (iring*iring)*_fact2;
				phi = (iphi-0.5) * M_PI/2/iring;
			}
			else if (pix<(_npix-_ncap)) // Equatorial region
			{
				int ip  = pix - _ncap;
				int iring = ip/(4*_nside) + _nside; // counted from North pole
				int iphi  = ip%(4*_nside) + 1;
				// 1 if iring+nside is odd, 1/2 otherwise
				double fodd = ((iring+_nside)&1) ? 1 : 0.5;

				int nl2 = 2*_nside;
				z = (nl2-iring)*_fact1;
				phi = (iphi-fodd) * M_PI/nl2;
			}
			else // South Polar cap
			{
				int ip = _npix - pix;
				int iring = int(0.5*(1+isqrt(2*ip-1))); //counted from South pole
				int iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));

				z = -1.0 + (iring*iring)*_fact2;
				phi = (iphi-0.5) * M_PI/2/iring;
			}
		}

		int Order() const
		{
			return _order;
		}

		int Npix() const
		{
			return _npix;
		}


		int ang2pix_z_phi (double z, double phi) const
		{
			double za = std::fabs(z);
			double tt = fmodulo(phi,2*M_PI) * 2/M_PI; // in [0,4)

			if (za<=2./3.) // Equatorial region
			{
				double temp1 = _nside*(0.5+tt);
				double temp2 = _nside*z*0.75;
				int jp = int(temp1-temp2); // index of  ascending edge line
				int jm = int(temp1+temp2); // index of descending edge line

				// ring number counted from z=2/3
				int ir = _nside+ 1 + jp - jm; // in {1,2n+1}
				int kshift = 1-(ir&1); // kshift=1 if ir even, 0 otherwise

				int ip = (jp+jm- _nside+kshift+1)/2; // in {0,4n-1}
				ip = imodulo(ip,4* _nside);

				return _ncap + (ir-1)*4* _nside+ ip;
			}
		else  // North & South polar caps
		{
			double tp = tt-int(tt);
			double tmp = _nside*sqrt(3*(1-za));

			int jp = int(tp*tmp); // increasing edge line index
			int jm = int((1.0-tp)*tmp); // decreasing edge line index

			int ir = jp+jm+1; // ring number counted from the closest pole
			int ip = int(tt*ir); // in {0,4*ir-1}
			ip = imodulo(ip,4*ir);

			if (z>0)
			{
				return 2*ir*(ir-1) + ip;
			}
			else
			{
				return _npix - 2*ir*(ir+1) + ip;
			}
		}

		}

		vec3 pix2vec (int pix) const
		{
			double z, phi;
			pix2ang_z_phi (pix,z,phi);
			vec3 res;
			res.set_z_phi (z, phi);
			return res;
		}

		int vec2pix (const vec3 &vec) const
		{ return ang2pix_z_phi (vec.z/vec.Length(), safe_atan2(vec.y,vec.x)); }


	private:
		int _order; // The Healpix Order (Sometimes refered to as _NSide)
		int _npix;  // The number of pixels
		int _ncap;
		int _nside;
		int _npface;
		double _fact1;
		double _fact2;

		// Returns the remainder of the division \a v1/v2.
		// The result is non-negative.
		// \a v1 can be positive or negative; \a v2 must be positive. 
		// FROM HEALPIX, cxxsupport/cxxutils.h
		inline double fmodulo (double v1, double v2) const
		{
			return (v1>=0) ? ((v1<v2) ? v1 : std::fmod(v1,v2)) : (fmod(v1,v2)+v2);
		}

		//! Returns the remainder of the division \a v1/v2.
		/*! The result is non-negative.
    \a v1 can be positive or negative; \a v2 must be positive. */
		// FROM HEALPIX, cxxsupport/cxxutils.h
		template<typename I> inline I imodulo (I v1, I v2) const
		{ return (v1>=0) ? ((v1<v2) ? v1 : (v1%v2)) : ((v1%v2)+v2); }

		//! Returns the integer \a n, which fulfills \a n*n<=arg<(n+1)*(n+1).
		// FROM HEALPIX, cxxsupport/cxxutils.h
		template<typename I> inline unsigned int isqrt (I arg) const
		{
		using namespace std;
		if (sizeof(I)<=4)
			return unsigned (sqrt(arg+0.5));
		long double arg2 = arg;
		return unsigned (sqrt(arg2+0.5));
		}
	protected:
		//! Returns \a atan2(y,x) if \a x!=0 or \a y!=0; else returns 0.
		inline double safe_atan2 (double y, double x) const
		{
			using namespace std;
			return ((x==0.) && (y==0.)) ? 0.0 : atan2(y,x);
		}

};

} // namespace sHEALPIX

#endif // SHEALPIX_HH
