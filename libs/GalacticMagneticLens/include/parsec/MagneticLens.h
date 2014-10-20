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

#include "parsec/ModelMatrix.h"
#include "parsec/Pixelization.h"

#include <boost/numeric/ublas/operation.hpp>

#include <vector>
#include <string>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdint.h>


namespace parsec
{

/// Holds one matrix for the lens and information about the rigidity range
class LensPart
{
	string _filename;
	double _rigidityMin;
	double _rigidityMax;
	ModelMatrix M;
	double _maximumSumOfColumns;
	bool _maximumSumOfColumns_calculated;

public:
	LensPart()
	{
	}
	/// File containing the matrix to be used in the range rigidityMin,
	/// rigidityMax with logarithmic units	[log10(E / eV)]
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
		M.deserialize(_filename);
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
			_maximumSumOfColumns = M.getMaximumOfSumsOfColumns();
			_maximumSumOfColumns_calculated = true;
		}
		return _maximumSumOfColumns;
	}
 
	/// Returns the minimum of the rigidity range for the lenspart in [log10(E/ eV)]
	double getMinimumRigidity()
	{
		return _rigidityMin;
	}

	/// Returns the maximum of the rigidity range for the lenspart in [log10(E/ eV)]
	double getMaximumRigidity()
	{
		return _rigidityMax;
	}

	/// Returns the modelmatrix
	ModelMatrix& getMatrix()
	{
		return M;
	}

	/// Sets the modelmatrix
	void setMatrix(const ModelMatrix& m)
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
	/// Loads part of a lens (one matrix) from file to use it in given rigidity range.
	void loadLensPart(const string &filename, double rigidityMin,
			double rigidityMax);

	// Stores the individual lenses
	std::vector<LensPart*> _lensParts;
	Pixelization* _pixelization;
	// Checks Matrix, raises Errors if not ok - also generate
	// _pixelization if called first time
	void _checkMatrix(const ModelMatrix &M);
	// minimum / maximum rigidity [log10(E / eV)] that is covered by the lens
	double _minimumRigidity;
	double _maximumRigidity;
	static bool _randomSeeded;
	double _norm;

public:
	/// Default constructor
	MagneticLens() :
			_pixelization(NULL), _minimumRigidity(-1), _maximumRigidity(-1), _norm(1)
	{
	}

	/// Constructs lens with predefined healpix order
	MagneticLens(uint8_t healpixorder) :
			_pixelization(NULL), _minimumRigidity(-1), _maximumRigidity(-1)
	{
		_pixelization = new Pixelization(healpixorder);
	}

	/// Construct lens and load lens from file
	MagneticLens(const string &filename) :
			_pixelization(NULL), _minimumRigidity(-1), _maximumRigidity(-1)
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
	/// Rigidity is given in EeV, phi and theta in rad
	bool transformCosmicRay(double rigidity, double& phi, double& theta);

	/// transforms the model array assuming that model points to an array of the
	/// correct size
	void transformModelVector(double* model, double rigidity) const;

	/// Loads M as part of a lens and use it in given rigidity range with
	/// given in logarithmic units [log10(E / eV)]
	void setLensPart(const ModelMatrix &M, double rigidityMin, double rigidityMax);

	/// Loads a lens from a given file, containing lines like
	/// lensefile.MLDAT rigidityMin rigidityMax
	void loadLens(const string &filename);

	/// Normalizes the lens parts to the maximum of sums of columns of
	/// every lenspart. By doing this, the lens won't distort the spectrum
	void normalizeLens();

	/// Normalizes the lens parts individually. Normalized this way, the
	/// lens generally distorts the spectrum of the sources, but deflects
	/// the UHECR more efficiently.
	void normalizeLensparts();

	/// Checks if rigidity [log10(E / eV)] is covered by lens
	bool rigidityCovered(double rigidity) const;

	/// Normalizes all matrix columns - the lens will then create fake
	/// anisotropies, but won't drop particles 
	void normalizeMatrixColumns();

	/// Returns minimum rigidity covered by lens, in EeV
	double getMinimumRigidity() const
	{
		return pow(10, _minimumRigidity - 18);
	}
	/// Returns maximum rigidity covered by lens, in EeV
	double getMaximumRigidity() const
	{
		return pow(10, _maximumRigidity - 18);
	}

	//	returns the norm used for the lenses
	double getNorm()
	{
		return _norm;
	}

	/// Returns iterator to the lens part with rigidity [log10(E / eV)]
	LensPart* getLensPart(double rigidity) const;

	/// Returns all lens parts
	const std::vector<LensPart*>& getLensParts() const
	{
		return _lensParts;
	}

	//calculates the mean deflection at given rigidity [EeV] 
	//double getMeanDeflection(double rigidity) const
	//{
	//	LensPart* lp = getLensPart(log10(rigidity) + 18);
	//	return calculateMeanDeflection(lp->getMatrix(), *_pixelization);
	//}

};


	


} // namespace

#endif // MAGNETICLENS_HH 
