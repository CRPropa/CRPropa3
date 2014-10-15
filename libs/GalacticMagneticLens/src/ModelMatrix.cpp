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

#include "parsec/ModelMatrix.h"
#include <ctime>
namespace parsec
{

void ModelMatrix::serialize(const string &filename)
{
	ofstream outfile(filename.c_str(), ios::binary);
	if (!outfile)
	{
		throw runtime_error("Can't write file: " + filename);
	}
	uint32_t nnz = 0;
	uint32_t C = 0;
	double val;

	// first write dummy and jump back later
	outfile.write((char*) &nnz, sizeof(uint32_t));
	C = (uint32_t) (this->size1());
	outfile.write((char*) &C, sizeof(uint32_t));
	C = (uint32_t) (this->size2());
	outfile.write((char*) &C, sizeof(uint32_t));

	// serialize non zero elements
	for (i2_t i2 = this->begin2(); i2 != this->end2(); ++i2)
	{
		for (i1_t i1 = i2.begin(); i1 != i2.end(); ++i1)
		{
			val = *i1;
			if (val > DBL_EPSILON)
			{
				C = (uint32_t) i1.index1();
				outfile.write((char*) &C, sizeof(uint32_t));
				C = (uint32_t) i1.index2();
				outfile.write((char*) &C, sizeof(uint32_t));
				outfile.write((char*) &val, sizeof(val));
				nnz++;
				if (outfile.fail())
				{
					throw runtime_error("Error writing file: " + filename);
				}
			}
		}
	}
	// jump back and write nnz
	outfile.seekp(0);
	outfile.write((char*) &nnz, sizeof(uint32_t));
	if (!outfile)
	{
		throw runtime_error("Error writing file: " + filename);
	}
	outfile.close();
}

void ModelMatrix::deserialize(const string &filename)
{
	ifstream infile(filename.c_str(), ios::binary);
	if (!infile)
	{
		throw runtime_error("Can't read file: " + filename);
	}

	uint32_t nnz, size1, size2;
	double val;
	infile.read((char*) &nnz, sizeof(uint32_t));
	infile.read((char*) &size1, sizeof(uint32_t));
	infile.read((char*) &size2, sizeof(uint32_t));
	boost::numeric::ublas::compressed_matrix<double,
			boost::numeric::ublas::column_major> M(size1, size2, nnz);

	for (size_t i = 0; i < nnz; i++)
	{
		infile.read((char*) &size1, sizeof(uint32_t));
		infile.read((char*) &size2, sizeof(uint32_t));
		infile.read((char*) &val, sizeof(double));
		//M(size1,size2) = val;
		M.push_back(size1, size2, val);
	}
	this->ModelMatrixType::operator=(M);
}

void ModelMatrix::normalizeRows()
{
	boost::numeric::ublas::diagonal_matrix<double,
			boost::numeric::ublas::row_major> A(this->size1(), this->size2());
	double rn;
	for (size_t i = 0; i < this->size1(); i++)
	{
		rn = norm_1(row((*this), i));
		A(i, i) = 1 / rn;
	}
	*this = prod(A, (*this));
}

void ModelMatrix::normalizeColumns(){
	boost::numeric::ublas::diagonal_matrix<double,boost::numeric::ublas::row_major> A(this->size1(),this->size2());
	double rn;
	ModelVectorType v;
	for (size_t i=0;i<this->size1();i++)
	{
		v = column(*this,i);
		rn = norm_1(v);
		v/=rn;
	}
}

double ModelMatrix::getSumOfColumn(size_t j) const
{
	ModelVectorType v;
	v = column(*this, j);
	double sum = norm_1(v);
	return sum;
}


double ModelMatrix::getMaximumOfSumsOfColumns() const
{
	double summax = 0;
	double sum = 0;
	for (size_t i = 0; i < this->size2(); i++)
	{
		sum = getSumOfColumn(i);
		if (sum > summax)
			summax = sum;
	}
	return summax;
}

void ModelMatrix::normalizeMatrix(double factor)
{
	*this /= factor;
}

} // namespace parsec
