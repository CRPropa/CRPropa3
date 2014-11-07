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

#include "crpropa/magneticLens/Pixelization.h"
#include "crpropa/Random.h"

namespace crpropa 
{

	healpix::T_Healpix_Base<int64_t> Pixelization::_healpix_nest = healpix::T_Healpix_Base<int64_t>(29, healpix::NEST);


uint8_t Pixelization::pix2Order(uint32_t pix)
{
	for (uint8_t i = 0; i < _nOrder_max; i++)
	{
		if (pix == _nPix[i])
			return i + 1;
	}
	return 0;
}

uint32_t Pixelization::nPix(uint8_t order)
{
	if (order > _nOrder_max)
	{
		return 0;
	}
	else
	{
		return _nPix[order - 1];
	}
}

uint32_t Pixelization::direction2Pix(double longitude, double latitude) const
{
	healpix::vec3 v;
	spherCo2Vec(longitude, latitude, v);
	try
	{
		uint32_t i = (uint32_t) _healpix->vec2pix(v);
		return i;
	}
	catch (healpix::PlanckError &e)
	{
		std::cerr << "Healpix error triggered from direction2Pix(" << longitude << ", " << latitude  << ")\n";
		std::cerr << " v = " << v.x <<", " << v.y << ", " <<  v.z << std::endl;
		std::cerr << "\n The original exception reads:\n";
		std::cerr << e.what() << std::endl;
		throw;
	}
}

void Pixelization::pix2Direction(uint32_t i, double &longitude,
		double &latitude) const
{
	healpix::vec3 v;
	try{
		v = _healpix->pix2vec(i);
	}
	catch (healpix::PlanckError &e)
	{
		std::cerr << "Healpix error triggered from pix2Direction(" << i << ", &longitude, &latitude " << ")\n";
		std::cerr << "The original exception reads:\n";
		std::cerr << e.what() << std::endl;
		throw;
	}

	vec2SphereCo(longitude, latitude, v);
}

void Pixelization::spherCo2Vec(double phi, double theta,
		healpix::vec3 &V) const
{
	V.x = cos(phi) * cos(theta);
	V.y = sin(phi) * cos(theta);
	V.z = sin(theta);
}

void Pixelization::vec2SphereCo(double &phi, double &theta,
		const healpix::vec3 &V) const
{
	theta = asin(V.z);
	phi = healpix::safe_atan2(V.y, V.x);
}


double Pixelization::angularDistance(uint32_t i, uint32_t j) const
{
	healpix::vec3 v1, v2;
	v1 = _healpix->pix2vec(i);
	v2 = _healpix->pix2vec(j);
	double s = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
	// Failsafe for numerical inaccuracies
	return ((s > 1) ? 0 : ((s < -1) ? M_PI : acos(s)));
}

	void Pixelization::getRandomDirectionInPixel(uint32_t i, double &longitude, double &latitude) 
	{
		
    int64_t inest = _healpix->ring2nest(i);

    int64_t nUp = 29 - _healpix->Order();
    int64_t iUp = inest * pow(4, nUp);
    iUp += Random::instance().randInt(pow(4, nUp));

		healpix::vec3 v = _healpix_nest.pix2vec(iUp);
		
		vec2SphereCo(longitude, latitude, v);
	}
} // namespace
