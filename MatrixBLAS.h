/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_MATRIX_BLAS_HEADER
#define SPIN_MATRIX_BLAS_HEADER
#define WANT_STREAM
#include <newmat.h>
#include <newmatio.h>
#include <cassert>
#include <cstdio>
#include <vector>
#include <map>
#include <fstream>
#include "blas_calls.h"
#include "ObjectMatrix.h"

using namespace std;

namespace SpinAdapted{

void svd(Matrix& M, DiagonalMatrix& e, Matrix& U, Matrix& V);
void xsolve_AxeqB(const Matrix& a, const ColumnVector& b, ColumnVector& x);

void Randomise (Matrix& a);
void SymmetricRandomise (Matrix& a);
void MatrixTensorProduct (const Matrix& a, char conjA, Real scaleA, const Matrix& b, char conjB, Real scaleB, Matrix& c, int rostride, int colstride, bool allocate = false);
void MatrixMultiply (const Matrix& a, char conjA, const Matrix& b, char conjB, Matrix& c, Real scale, double cfactor = 1.);
void MatrixMultiply (double d, const Matrix& a, Matrix& b);
void MatrixScaleAdd (double d, const Matrix& a, Matrix& b);
void MatrixScale(double, Matrix& a);
double MatrixDotProduct(const Matrix& a, const Matrix& b);
void MatrixNormalise(Matrix& a);
void MatrixDiagonalScale(double d, const Matrix& a, double* b);




double dotproduct(const ColumnVector& a, const ColumnVector& b);
double dotproduct(const RowVector& a, const RowVector& b);
double rowdoubleproduct(Matrix& a, int rowa, Matrix& b, int rowb);
void diagonalise(Matrix& sym, DiagonalMatrix& d, Matrix& vec);
void diagonalise_tridiagonal(std::vector<double>& diagonal, std::vector<double>& offdiagonal, int numelements, Matrix& vec);

template<class T> void Clear (T& a)
{
#ifdef BLAS
  memset(a.Store(), 0, a.Storage() * sizeof(double));
#else
    a = 0.;
#endif
}

void CatenateProduct (const ObjectMatrix<Matrix*>& a, Matrix& b, bool allocate = false);  
char TransposeOf (char c);
inline char TransposeOf (char c) { assert (c == 'n' || c == 't'); return (c == 'n') ? 't' : 'n'; }
void MatrixRotate (const Matrix& a, const Matrix& b, const Matrix& c, Matrix& d);
void Save (const Matrix& a, std::ofstream &ofs);
void Load (Matrix& a, std::ifstream &ifs);

template<class T> void print(ostream& os, const std::vector<T>& v)
{
  for (int i = 0; i < v.size(); ++i)
    os << v[i] << " ";
  os << endl;
}
double CheckSum (Matrix& a);
void DebugPrint (std::vector<int>& v);
void DebugPrint (std::vector<double>& v);


template<class T> void CyclicPermute (std::vector<T>& v, std::vector<T>& result, int n);

template<class T> void CyclicPermute (std::vector<T>& v, std::vector<T>& result, int n)
{
  assert (n >= 0);
  result.resize (v.size ());
  for (int i = 0; i < v.size (); ++i)
    result [i] = v [(i + n) % v.size ()];
}

template<class T> bool find(vector<T>& v, const T& value) { return std::find(v.begin(), v.end(), value) != v.end(); }

template<class T> void get_sorted_indices(const std::vector<T>& data, std::vector<int>& indices)
{
  multimap<T, int> sorted_data;
  for (int i = 0; i < data.size(); ++i)
    sorted_data.insert(make_pair<T, int>(data[i], i));
  typename multimap<T, int>::iterator it = sorted_data.begin();
  indices.reserve(data.size());
  while (it != sorted_data.end())
  {
	  indices.push_back(it->second);
	  ++it;
  }

}
}
#endif
