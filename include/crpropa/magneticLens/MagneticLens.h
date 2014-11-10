//----------------------------------------------------------------------
// This file is part of PARSEC (http://physik.rwth-aachen.de/parsec)
// a parametrized simulation engine for cosmic rays.
//
// Copyright (C) 2011	Martin Erdmann, Peter Schiffer, Tobias Winchen
//										 RWTH Aachen University, Germany
// Contact: winchen@physik.rwth-aachen.de
//
//	This program is free software: you can redistribute it and/or
//	modify it under the terms of the GNU General Public License as
//	published by the Free Software Foundation, either version 3 of
//	the License, or (at your option) any later version.
//
//	This program is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with this program. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#ifndef MAGNETICLENS_HH 
#define MAGNETICLENS_HH 

#include "crpropa/magneticLens/ModelMatrix.h"
#include "crpropa/magneticLens/Pixelization.h"
#include "crpropa/Units.h"
#include "crpropa/Vector3.h"

#include <vector>
#include <string>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdint.h>


namespace crpropa 
{

/// Holds one matrix for the lens and information about the rigidity range
class LensPart
{
	string _filename;
	double _rigidityMin;
	double _rigidityMax;
	ModelMatrixType M;
	double _maximumSumOfColumns;
	bool _maximumSumOfColumns_calculated;

public:
	LensPart()
	{
	}
	/// File containing the matrix to be used in the range rigidityMin,
	/// rigidityMax in Joule 
	LensPart(const std::string &filename, double rigidityMin, double rigidityMax) :
			_filename(filename), _rigidityMin(rigidityMin), _rigidityMax(rigidityMax), _maximumSumOfColumns_calculated(
					false), _maximumSumOfColumns(0)
	{
	}

	~LensPart()
	{
	}

	/// Loads the matrix from file
	void loadMatrixFromFile()
	{
		deserialize(_filename, M);
	}

	/// Returns the filename of the matrix
	const std::string& getFilename()
	{
		return _filename;
	}

	/// Calculates the maximum of the sums of columns for the matrix
	double getMaximumOfSumsOfColumns()
	{
		if (!_maximumSumOfColumns_calculated)
		{ // lazy calculation of maximum
			_maximumSumOfColumns = maximumOfSumsOfColumns(M);
			_maximumSumOfColumns_calculated = true;
		}
		return _maximumSumOfColumns;
	}
 
	/// Returns the minimum of the rigidity range for the lenspart in eV 
	double getMinimumRigidity()
	{
		return _rigidityMin / eV;
	}

	/// Returns the maximum of the rigidity range for the lenspart in eV
	double getMaximumRigidity()
	{
		return _rigidityMax / eV;
	}

	/// Returns the modelmatrix
	ModelMatrixType& getMatrix()
	{
		return M;
	}

	/// Sets the modelmatrix
	void setMatrix(const ModelMatrixType& m)
	{
		M = m;
	}

	
};

/// Function to calculate the mean deflection [rad] of the matrix M, given a pixelization
//double calculateMeanDeflection(const ModelMatrix &M,
//		const Pixelization &pixelization)
//{
//	double totalDeflection = 0;
//	double weightSum = 0;
//	for (const_i2_t it1 = M.begin2(); it1 != (M.end2()); it1++)
//	{
//		for (const_i1_t it2 = it1.begin(); it2 != it1.end(); it2++)
//		{
//			totalDeflection+= pixelization.angularDistance(it2.index1(),
//					it2.index2()) * (*it2) ;
//			weightSum+= (*it2);
//		}
//	}
//	return totalDeflection / weightSum;
//}

typedef std::vector<LensPart*>::iterator LensPartIter;
typedef std::vector<LensPart*>::const_iterator const_LensPartIter;

/// The lens for the galactic magnetic field.
/// Note that the energies refer to protons (Z=1). To be used with other particles with a different charge number please select the rigidity accordingly.
class MagneticLens
{
	
	void updateRigidityBounds(double rigidityMin, double rigidityMax);
	
	/// Loads part of a lens (one matrix) from file to use it in given rigidity range.
	void loadLensPart(const string &filename, double rigidityMin,
			double rigidityMax);

	// Stores the individual lenses
	std::vector<LensPart*> _lensParts;
	Pixelization* _pixelization;
	// Checks Matrix, raises Errors if not ok - also generate
	// _pixelization if called first time
	void _checkMatrix(const ModelMatrixType &M);
	// minimum / maximum rigidity that is covered by the lens [Joule]
	double _minimumRigidity;
	double _maximumRigidity;
	static bool _randomSeeded;
	double _norm;

public:
	/// Default constructor
	MagneticLens() :
			_pixelization(NULL), _minimumRigidity(DBL_MAX), _maximumRigidity(DBL_MIN), _norm(1)
	{
	}

	/// Constructs lens with predefined healpix order
	MagneticLens(uint8_t healpixorder) :
			_pixelization(NULL), _minimumRigidity(DBL_MAX), _maximumRigidity(DBL_MIN)
	{
		_pixelization = new Pixelization(healpixorder);
	}

	/// Construct lens and load lens from file
	MagneticLens(const string &filename) :
			_pixelization(NULL), _minimumRigidity(DBL_MAX), _maximumRigidity(DBL_MIN)
	{
		loadLens(filename);
	}

	/// Returns the pixelization used
	const Pixelization& getPixelization() const
	{
		return (*_pixelization);
	}

	/// Default destructor
	~MagneticLens()
	{
		if (_pixelization)
			delete _pixelization;
		for (std::vector<LensPart*>::iterator iter = _lensParts.begin(); 
				iter != _lensParts.end(); iter++)
		{
			delete (*iter);
		}
		_lensParts.clear();
	}

	/// Try to transform the comsic ray to a new direction.
	/// Returns false and does not change phi and theta if the cosmic ray is
	/// lost due to conservation of cosmic ray flux.
	/// Rigidity is given in Joule, phi and theta in rad
	bool transformCosmicRay(double rigidity, double& phi, double& theta);

	/// Tries transform a cosmic ray with momentum vector p 
	bool transformCosmicRay(double rigidity, Vector3d &p);

	/// transforms the model array assuming that model points to an array of the
	/// correct size. Rigidity is given in Joule
	void transformModelVector(double* model, double rigidity) const;

	/// Loads M as part of a lens and use it in given rigidity range with
	/// rigidities given in Joule 
	void setLensPart(const ModelMatrixType &M, double rigidityMin, double rigidityMax);

	/// Loads a lens from a given file, containing lines like
	/// lensefile.MLDAT rigidityMin rigidityMax
	/// rigidities are given in logarithmic units [log10(E / eV)]
	void loadLens(const string &filename);

	/// Normalizes the lens parts to the maximum of sums of columns of
	/// every lenspart. By doing this, the lens won't distort the spectrum
	void normalizeLens();

	/// Normalizes the lens parts individually. Normalized this way, the
	/// lens generally distorts the spectrum of the sources, but deflects
	/// the UHECR more efficiently.
	void normalizeLensparts();

	/// Checks if rigidity [Joule] is covered by lens
	bool rigidityCovered(double rigidity) const;

	/// Normalizes all matrix columns - the lens will then create fake
	/// anisotropies, but won't drop particles 
	void normalizeMatrixColumns();

	/// Returns minimum rigidity covered by lens, in eV
	double getMinimumRigidity() const
	{
		return _minimumRigidity / eV;
	}
	/// Returns maximum rigidity covered by lens, in eV
	double getMaximumRigidity() const
	{
		return _maximumRigidity / eV;
	}

	//	returns the norm used for the lenses
	double getNorm()
	{
		return _norm;
	}

	/// Returns iterator to the lens part with rigidity Joule 
	LensPart* getLensPart(double rigidity) const;

	/// Returns all lens parts
	const std::vector<LensPart*>& getLensParts() const
	{
		return _lensParts;
	}
};


	


} // namespace

#endif // MAGNETICLENS_HH 
