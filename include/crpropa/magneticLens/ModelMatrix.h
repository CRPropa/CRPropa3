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

#include <Eigen/SparseCore>

namespace crpropa 
{

	typedef Eigen::SparseMatrix<double> ModelMatrixType;
	typedef Eigen::SparseVector<double> ModelVectorType;

	/// Writes the ModelMatrix to disk as binary files with the format:
	/// Int (number of non zero elements), Int (size1), Int (size2)
	/// (Int, Int, Double) : (column, row, value) triples ...
	void serialize(const string &filename, const ModelMatrixType &matrix);

	/// Reads a matrix from file
	void deserialize(const string &filename, ModelMatrixType &matrix);

	/// Normalizes each column j of the matrix so that, \f$ \Vert m_j \Vert_1 = 1 \f$ 
	void normalizeColumns(ModelMatrixType &matrix);

	/// Calculate the maximum of the unity norm of the column vectors of the matrix \f$\max_j(\Vert m_j \Vert_1) \f$
	double maximumOfSumsOfColumns(const ModelMatrixType &matrix);
	
	double norm_1(const ModelVectorType &v);

	void normalizeMatrix(ModelMatrixType& matrix, double norm);

	// matrix vector product with update: model = matrix * model
	void prod_up(const ModelMatrixType& matrix, double* model);
} // namespace parsec

#endif // MODELMATRIX_HH
