//----------------------------------------------------------------------
// This file is part of PARSEC (http://physik.rwth-aachen.de/parsec)
// a parametrized simulation engine for cosmic rays.
//
// Copyright (C) 2011  Martin Erdmann, Peter Schiffer, Tobias Winchen
//                     RWTH Aachen University, Germany
// Contact: winchen@physik.rwth-aachen.de
//
//  This program is free software: you can redistribute it and/or
//  modify it under the terms of the GNU General Public License as
//  published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------


#ifndef PIXELIZATION_HH 
#define PIXELIZATION_HH

#include "parsec/shealpix.h"
#include <cmath>
#include <stdint.h>

namespace parsec
{

/// Helpers to makes work with Healpix smooth 
const uint8_t _nOrder_max = 13;
const uint32_t _nPix[] =
{
		48,
		192,
		768,
		3072,
		12288,
		49152,
		196608,
		786432,
		3145728,
		12582912,
		50331648,
		201326592,
		805306368
};

/// Every communication with healpix is done through this class to avoid
/// bugs with missmatching coordinates (and make python hooks easier)
class Pixelization : private shealpix::sHEALPIX
{
public:
	Pixelization() : shealpix::sHEALPIX( )
	{
	}

	/// Constructor creating Pixelization with healpix order 6 (about
	/// 50000 pixels)
	Pixelization(uint8_t order) : shealpix::sHEALPIX(order)
	{
	}

	/// Returns the number of the pixel which includes the direction (phi,theta) 
	/// phi in [-pi, pi], theta in [-pi/2, pi/2]
	uint32_t direction2Pix(double longitude, double latitude) const;

	/// Returns the number of pixels of the pixelization
	uint32_t nPix() const
	{
		return Npix();
	}

	/// Returns the number of pixels given by healpix order
	static uint32_t nPix(uint8_t order);

	/// Returns the number of pixels of the pixelization
	int getNumberOfPixels()
	{
		return Npix();
	}

	/// Returns the order, a given pixel number corresponds to. 0 if no
	/// match!
	static uint8_t pix2Order(uint32_t pix);

	/// Gives the center of pixel i in longitude [rad] and latitude [rad]
	void pix2Direction(uint32_t i, double &longitude, double &latitude) const;

	/// Calculate the angle [rad] between the vectors pointing to pixels i and j
	double angularDistance(uint32_t i, uint32_t j) const;

	/// Returns the maximum possible pixelization order
	uint8_t getMaxOrder() const
	{
		return _nOrder_max;
	}

  /// Returns healpix order
  uint8_t getOrder() const
  {
    return Order();
  }


private:
	void spherCo2Vec(double phi, double theta, shealpix::vec3 &V) const;
	void vec2SphereCo(double &phi , double &theta, const shealpix::vec3 &V) const;
};


} // namespace
#endif // PIXELIZATION_HH 
