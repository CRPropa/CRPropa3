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

#ifndef MODELMATRIX_HH
#define MODELMATRIX_HH

#define BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR 1
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

//#include "crpropa/Random.h"

#include <ctime>
#include <ctime>
#include <cmath>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cfloat>
#include <stdexcept>
#include <vector>
#include <stdint.h>

using namespace std;

namespace crpropa 
{

/// ModelMatrixType specifies the used Matrix Type
/// Changes here don't break compability, as matrices are stored in a
/// selfmade format. Compressed_matrix however has been tested as the
/// fastest type
typedef boost::numeric::ublas::compressed_matrix<double,
		boost::numeric::ublas::column_major> ModelMatrixType;
typedef ModelMatrixType::iterator1 i1_t;
typedef ModelMatrixType::iterator2 i2_t;

typedef ModelMatrixType::const_iterator1 const_i1_t;
typedef ModelMatrixType::const_iterator2 const_i2_t;

typedef boost::numeric::ublas::compressed_vector<double> ModelVectorType;
typedef ModelVectorType::iterator MVi_t;

/// Adaptors for direct array access
typedef boost::numeric::ublas::shallow_array_adaptor<double> shallow_adaptor_double;
typedef boost::numeric::ublas::vector<double, shallow_adaptor_double> shallow_vector_double;

/// The class holds a Magnetic Field Model as a Matrix, with m_i,j is
/// the probability that a particle from pixel j reaches pixel i
class ModelMatrix: public ModelMatrixType
{
public:
	ModelMatrix() :
			ModelMatrixType()
	{
	}

	/// Creates a modelmatrix M with size1 x size2 elements as sparse matrix
	/// with nnz non-zero elements
	ModelMatrix(uint32_t size1, uint32_t size2, uint32_t nnz) :
			ModelMatrixType(size1, size2, nnz)
	{
	}

	/// Writes the ModelMatrix to disk as binary files with the format:
	/// Int (number of non zero elements), Int (size1), Int (size2)
	/// (Int, Int, Double) : (column, row, value) triples ...
	void serialize(const string &filename);

	/// Reads a matrix from file
	void deserialize(const string &filename);

	/// Allow construction from matrix
	ModelMatrix& operator=(const ModelMatrixType & source)
	{
		if (this != &source)
		{
			this->ModelMatrixType::operator=(source);
		}
		return *this;
	}

	/// Normalizes each row m_i, so that \f$ \Vert m_i \Vert_1 = 1 \f$
	/// Needed to ensure that each pixel on earth contains the same amount of backtracked particles.
	void normalizeRows();

	/// Normalizes each column j of the matrix so that, \f$ \Vert m_j \Vert_1 = 1 \f$ 
	void normalizeColumns();

	/// Calculate the maximum of the unity norm of the column vectors of the matrix \f$\max_j(\Vert m_j \Vert_1) \f$
	double getMaximumOfSumsOfColumns() const;
	
	// get sum of column j
	double getSumOfColumn(size_t j) const;

	/// Divides matrix by the factor f
	void normalizeMatrix(double factor);

	/// Sets element \f$ m_{i,j} \f$ to value
	void setElement(size_t i, size_t j, double value)
	{
		this->ModelMatrixType::operator()(i, j) = value;
	}

	/// Returns element \f$ m_{i,j} \f$
	double getElement(size_t i, size_t j) const
	{
		return this->ModelMatrixType::operator()(i, j);
	}


};


} // namespace parsec

#endif // MODELMATRIX_HH
