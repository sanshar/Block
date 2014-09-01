/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include <iostream>
#include <cmath>
#include <include/newmatutils.h>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include "MatrixBLAS.h"
#include "global.h"
#ifdef BLAS
#include "blas_calls.h"
#endif
#ifdef USING_DOUBLE
#define GAXPY DAXPY
#define GEMM DGEMM
#endif
#ifdef USING_FLOAT
#define GAXPY SAXPY
#define GEMM SGEMM
#endif

void SpinAdapted::MatrixScale(double d, Matrix& a)
{
#ifdef BLAS
  DSCAL(a.Storage(), d, a.Store(), 1);
#else
  a *= d;
#endif  
}
double SpinAdapted::MatrixDotProduct(const Matrix& a, const Matrix& b)
{
  assert((a.Nrows() == b.Nrows()) && (a.Ncols() == b.Ncols()));
#ifdef BLAS
  return DDOT(a.Storage(), a.Store(), 1, b.Store(), 1);
#else
  abort();
#endif
}
void SpinAdapted::MatrixNormalise(Matrix& a)
{
  double norm = MatrixDotProduct(a, a);
  MatrixScale(1./sqrt(norm), a);
}

void SpinAdapted::Randomise (Matrix& a)
{
  Real* val = a.Store ();
  for (int i = 0; i < a.Storage (); ++i)
    {
      *val = double(rand ()) / RAND_MAX;
      ++val;
    }
}

void SpinAdapted::SymmetricRandomise (Matrix& a)
{
  assert(a.Nrows() == a.Ncols());

  for (int i=0; i<a.Nrows(); i++)
  for (int j=0; j<i+1; j++) {
    a(i+1, j+1) = double(rand())/RAND_MAX;
    a(j+1, i+1) = a(i+1, j+1);
  }

}

double SpinAdapted::dotproduct(const ColumnVector& a, const ColumnVector& b)
{
  assert(a.Nrows() == b.Nrows());
#ifdef BLAS
  return DDOT(a.Storage(), a.Store(), 1, b.Store(), 1);
#else
  return a.t() * b;
#endif
}

double SpinAdapted::dotproduct(const RowVector& a, const RowVector& b)
{
  assert(a.Ncols() == b.Ncols());
#ifdef BLAS
  return DDOT(a.Storage(), a.Store(), 1, b.Store(), 1);
#else
  return a * b.t();
#endif
}

double SpinAdapted::rowdoubleproduct(Matrix& a, int rowa, Matrix& b, int rowb)
{
  assert(a.Ncols() == b.Ncols());
  double* aptr = a.Store() + a.Ncols() * rowa;
  double* bptr = b.Store() + b.Ncols() * rowb;
  return DDOT(a.Ncols(), aptr, 1, bptr, 1);
}

void SpinAdapted::MatrixScaleAdd (double d, const Matrix& a, Matrix& b)
{
  assert (a.Nrows () == b.Nrows () && a.Ncols () == b.Ncols ());
#ifdef BLAS
  int n = a.Nrows () * a.Ncols ();
  assert (n == (b.Nrows () * b.Ncols ()));
  DAXPY (n, d, a.Store (), 1, b.Store (), 1);
#else
  b += d * a;
#endif
}

void SpinAdapted::MatrixDiagonalScale(double d, const Matrix& a, double* b)
{
  //assert (a.Nrows () == a.Ncols () && a.Nrows () == b.Ncols ());
#ifdef BLAS
  int n = a.Nrows ();
  DAXPY (n, d, a.Store (), n+1, b, 1);
#else
  //b += d * a; Should add the non-blas analogue
#endif
}

void SpinAdapted::MatrixTensorProduct (const Matrix& a_ref, char conjA, Real scaleA, const Matrix& b_ref, char conjB, Real scaleB, Matrix& c, int rowstride, int colstride, bool allocate)
{
#ifndef BLAS
  Matrix A;
  Matrix B;
#endif
  Matrix& a = const_cast<Matrix&>(a_ref); // for BLAS calls
  Matrix& b = const_cast<Matrix&>(b_ref);

  int arows = a.Nrows();
  int acols = a.Ncols();
  
  // some specialisations
#ifdef FAST_MTP
  //  if ((brows == 1) && (bcols == 1))
    {
      double b00 = *b.Store();
      if (conjA == 'n')
	{
	  double* cptr = c.Store()+ rowstride*c.Ncols() + colstride;
	  for (int i=0; i< a.Nrows();i++) 
	    DAXPY(a.Ncols(), scaleA * scaleB * b00, a.Store()+i*a.Ncols(), 1, cptr + i*c.Ncols(), 1);
	  return;
	}
      else 	
	{
	  double* aptr = a.Store();
	  double* cptr = c.Store() + rowstride*c.Ncols() + colstride;
	  for (int col = 0; col < acols; ++col)
	    {
	      DAXPY(arows, scaleA * scaleB * b00, aptr, acols, cptr, 1);
	      ++aptr;
	      cptr += c.Ncols();//arows;
	    }

	  return;
	}	
    }
    //  else
    //    abort();
#else 
      try
	{
	  if (conjA == 'n' && conjB == 'n')
	    {
	      if (allocate)
		{
		  c.ReSize (a.Nrows () * b.Nrows (), a.Ncols () * b.Ncols ());
		  Clear (c);
		}
	      //assert ((c.Nrows () == (a.Nrows () * b.Nrows ())) && (c.Ncols () == (a.Ncols () * b.Ncols ())));
#ifdef BLAS
	      int aRows = a.Nrows ();
	      int aCols = a.Ncols ();
	      int bRows = b.Nrows ();
	      int bCols = b.Ncols ();

	      for (int i = 0; i < aRows; ++i)
		for (int j = 0; j < aCols; ++j)
		  {
		    Real scale = scaleA * scaleB * a (i+1,j+1);
		    for (int k = 0; k < bRows; ++k)
		      GAXPY (bCols, scale, &b (k+1,1), 1, &c (i * bRows + k+1 +rowstride,j * bCols+1+colstride), 1);
		  }
	      return;
#else
	      A = a;
	      B = b;
#endif
	    }
	  else if (conjA == 't' && conjB == 'n')
	    {
	      if (allocate)
		{
		  c.ReSize (a.Ncols () * b.Nrows (), a.Nrows () * b.Ncols ());
		  Clear (c);
		}
	      //assert ((c.Nrows () == (a.Ncols () * b.Nrows ())) && (c.Ncols () == (a.Nrows () * b.Ncols ())));
#ifdef BLAS
	      int aRows = a.Ncols ();
	      int aCols = a.Nrows ();
	      int bRows = b.Nrows ();
	      int bCols = b.Ncols ();
	      
	      for (int i = 0; i < aRows; ++i)
		for (int j = 0; j < aCols; ++j)
		  {
		    Real scale = scaleA * scaleB * a (j+1,i+1);
		    for (int k = 0; k < bRows; ++k)
		      GAXPY (bCols, scale, &b (k+1,1), 1, &c (i * bRows + k+1+rowstride,j * bCols+1+colstride), 1);
		  }
	      return;
#else	  
	      A = a.t ();
	      B = b;
#endif
	    }
	  else if (conjA == 'n' && conjB == 't')
	    {
	      if (allocate)
		{
		  c.ReSize (a.Nrows () * b.Ncols (), a.Ncols () * b.Nrows ());
		  Clear (c);
		}
	      //assert ((c.Nrows () == (a.Nrows () * b.Ncols ())) && (c.Ncols () == (a.Ncols () * b.Nrows ())));
#ifdef BLAS
	      int aRows = a.Nrows ();
	      int aCols = a.Ncols ();
	      int bRows = b.Ncols ();
	      int bCols = b.Nrows ();
	      
	      for (int i = 0; i < aRows; ++i)
		for (int j = 0; j < aCols; ++j)
		  {
		    Real scale = scaleA * scaleB * a (i+1,j+1);
		    for (int k = 0; k < bRows; ++k)
		      GAXPY (bCols, scale, &b (1,k+1), bRows, &c (i * bRows + k+1+rowstride,j * bCols+1+colstride), 1);
		  }
	      return;
#else
	      A = a;
	      B = b.t ();
#endif
	    }
	  else if (conjA == 't' && conjB == 't')
	    {
	      if (allocate)
		{
		  c.ReSize (a.Ncols () * b.Ncols (), a.Nrows () * b.Nrows ());
		  Clear (c);
		}
	      //assert ((c.Nrows () == (a.Ncols () * b.Ncols ())) && (c.Ncols () == (a.Nrows () * b.Nrows ())));
#ifdef BLAS
	      int aRows = a.Ncols ();
	      int aCols = a.Nrows ();
	      int bRows = b.Ncols ();
	      int bCols = b.Nrows ();
	      
	      for (int i = 0; i < aRows; ++i)
		for (int j = 0; j < aCols; ++j)
		  {
		    Real scale = scaleA * scaleB * a (j+1,i+1);
		    for (int k = 0; k < bRows; ++k)
		      GAXPY (bCols, scaleA * scaleB * a (j+1,i+1), &b (1,k+1), bRows, &c (i * bRows + k+1+rowstride,j * bCols+1+colstride), 1);
		  }
	      return;
#else
	      A = a.t ();
	      B = b.t ();
#endif
	    }
	  else
	    abort ();
#ifndef BLAS
	  for (int i = 1; i <= A.Nrows (); ++i)
	    for (int j = 1; j <= A.Ncols (); ++j)
	      c.SubMatrix ((i - 1) * B.Nrows () + 1, i * B.Nrows (), (j - 1) * B.Ncols () + 1, j * B.Ncols ()) += (scaleA * scaleB) * A (i,j) * B; 
#endif
	  
	}
      catch (Exception)
	{
	  pout << Exception::what () << endl;
	  abort ();
	}   
#endif
}

void SpinAdapted::xsolve_AxeqB(const Matrix& a, const ColumnVector& b, ColumnVector& x)
{
  FORTINT ar = a.Nrows();
  int bc = 1;
  int info=0;
  FORTINT* ipiv = new FORTINT[ar];
  double* bwork = new double[ar];
  for(int i = 0;i<ar;++i)
    bwork[i] = b.element(i);
  double* workmat = new double[ar*ar];
  for(int i = 0;i<ar;++i)
    for(int j = 0;j<ar;++j)
      workmat[i*ar+j] = a.element(j,i);

  GESV(ar, bc, workmat, ar, ipiv, bwork, ar, info);
  delete[] ipiv;
  delete[] workmat;

  for(int i = 0;i<ar;++i)
    x.element(i) = bwork[i];

  delete[] bwork;

  if(info != 0)
  {
     pout << "Xsolve failed with info error " << info << endl;
     abort();
  }
}

void SpinAdapted::svd(Matrix& M, DiagonalMatrix& d, Matrix& U, Matrix& V)
{
  int nrows = M.Nrows();
  int ncols = M.Ncols();

  assert(nrows >= ncols);

  int minmn = min(nrows, ncols);
  int maxmn = max(nrows, ncols);
  int eigenrows = min(minmn, minmn);
  d.ReSize(minmn);
  Matrix Ut;
  Ut.ReSize(nrows, nrows);
  V.ReSize(ncols, ncols);

  int lwork = maxmn * maxmn + 100;
  double* workspace = new double[lwork];

  // first transpose matrix
  Matrix Mt;
  Mt = M.t();
  int info = 0;
  DGESVD('A', 'A', nrows, ncols, Mt.Store(), nrows, d.Store(), 
	 Ut.Store(), nrows, V.Store(), ncols, workspace, lwork, info);

  U.ReSize(nrows, ncols);
  SpinAdapted::Clear(U);
  for (int i = 0; i < nrows; ++i)
    for (int j = 0; j < ncols; ++j)
      U(i+1,j+1) = Ut(j+1,i+1);
  delete[] workspace;
}

void SpinAdapted::diagonalise_tridiagonal(std::vector<double>& diagonal, std::vector<double>& offdiagonal, int numelements, Matrix& vec)
{
  int nrows = numelements;
  int ncols = numelements;

  vec.ReSize(nrows, nrows);

  Matrix vec_transpose; vec_transpose = vec;
  vector<double> workarray(4*nrows-2,0);
  int info = 0;

  DSTEV('V', nrows, &(diagonal[0]), &(offdiagonal[0]), vec_transpose.Store(), nrows,  &(workarray[0]), info);

  if (info != 0)
    {
      pout << "failed to converge :: " <<info<< endl;
      abort();
    }

  for (int i = 0; i < nrows; ++i)
    for (int j = 0; j < ncols; ++j)
      vec(j+1,i+1) = vec_transpose(i+1,j+1);
}


void SpinAdapted::diagonalise(Matrix& sym, DiagonalMatrix& d, Matrix& vec)
{
  int nrows = sym.Nrows();
  int ncols = sym.Ncols();
  assert(nrows == ncols);
  d.ReSize(nrows);
  vec.ReSize(nrows, nrows);

  Matrix workmat;
  workmat = sym;
  vector<double> workquery(1);
  int info = 0;
  double* dptr = d.Store();

  int query = -1;
  DSYEV('V', 'L', nrows, workmat.Store(), nrows, dptr, &(workquery[0]), query, info); // do query to find best size
  
  int optlength = static_cast<int>(workquery[0]);
  vector<double> workspace(optlength);

  DSYEV('V', 'U', nrows, workmat.Store(), nrows, dptr, &(workspace[0]), optlength, info); // do query to find best size


  
  if (info > 0) 
    {
      pout << "failed to converge " << endl;
      abort(); 
    }
  
  for (int i = 0; i < nrows; ++i)
    for (int j = 0; j < ncols; ++j)
      vec(j+1,i+1) = workmat(i+1,j+1);
}



void SpinAdapted::MatrixMultiply (double d, const Matrix& a, Matrix& b)
{
  //  b += d * a;
#ifdef BLAS 
  assert ((a.Nrows () == b.Nrows ()) && (a.Ncols () == b.Ncols ()));
  int n = a.Nrows () * a.Ncols ();
  GAXPY (n, d, a.Store (), 1, b.Store (), 1);
#else
  b += d * a;
#endif
}

double SpinAdapted::CheckSum (Matrix& a)
{
  double val = 0.;
  for (int i = 0; i < a.Nrows (); ++i)
    for (int j = 0; j < a.Ncols (); ++j)
      val += a.element (i, j);
  return val;
}
void SpinAdapted::MatrixMultiply (const Matrix& a, char conjA, const Matrix& b, char conjB, Matrix& c, Real scale, double cfactor)
{
  //dmrginp.justmultiply.start();
  //dmrginp.justmultiply -> start(); //ROA
  Matrix& a_ref = const_cast<Matrix&>(a); // for BLAS calls
  Matrix& b_ref = const_cast<Matrix&>(b);
  try
    {
      int aRows = a_ref.Nrows ();
      int aCols = a_ref.Ncols ();
      int bRows = b_ref.Nrows ();
      int bCols = b_ref.Ncols ();
      int cRows = c.Nrows ();
      int cCols = c.Ncols ();
      if (conjA == 'n' && conjB == 'n')
	{	  
	  assert ((aCols == bRows) && (cRows == aRows) && (cCols == bCols));
#ifdef BLAS
	  GEMM ('n', 'n', bCols, aRows, bRows, scale, b.Store (), bCols, a.Store (), aCols, cfactor, c.Store (), bCols);
#else
	  c += (scale * a) * b;
#endif
	}
      else if (conjA == 'n' && conjB == 't')
	{
	  assert ((aCols == bCols) && (cRows == aRows) && (cCols == bRows));
#ifdef BLAS
	  GEMM ('t', 'n', bRows, aRows, bCols, scale, b.Store (), bCols, a.Store (), aCols, cfactor, c.Store (), bRows);
#else
	  c += (scale * a) * b.t ();
#endif
	} 
      else if (conjA == 't' && conjB == 'n')
	{
	  assert ((aRows == bRows) && (cRows == aCols) && (cCols == bCols));
#ifdef BLAS
	  GEMM ('n', 't', bCols, aCols, bRows, scale, b.Store (), bCols, a.Store (), aCols, cfactor, c.Store (), bCols);
#else
	c += (scale * a.t ()) * b;
#endif
	}
      else if (conjA == 't' && conjB == 't')
	{
	  assert ((aRows == bCols) && (cRows == aCols) && (cCols == bRows));
#ifdef BLAS
	  GEMM ('t', 't', bRows, aCols, bCols, scale, b.Store (), bCols, a.Store (), aCols, cfactor, c.Store (), bRows);
#else
	  c += (scale * a.t ()) * b.t ();
#endif
	}
      else
	abort ();
    }
  catch (Exception)
    {
      pout << Exception::what () << endl;
      abort ();
   }
  //dmrginp.justmultiply.stop();
  //dmrginp.justmultiply -> stop(); //ROA
}

void SpinAdapted::CatenateProduct (const ObjectMatrix<Matrix*>& a, Matrix& b, bool allocate)
{
  try
    {
      std::vector<int> indexRows (a.Nrows ());
      std::vector<int> indexCols (a.Ncols ());
      int rowLength = 0;
      int colLength = 0;
      for (int i = 0; i < indexRows.size (); ++i)
	{
	  indexRows [i] = (i > 0) ? a (i - 1,0)->Nrows () + indexRows [i - 1] : 1;
	  rowLength += a (i,0)->Nrows ();
	}
      for (int i = 0; i < indexCols.size (); ++i)
	{
	  indexCols [i] = (i > 0) ? a (0,i - 1)->Ncols () + indexCols [i - 1] : 1;
	  colLength += a (0,i)->Ncols ();
	}
      
      if (!allocate) 
	assert (b.Nrows () == rowLength && b.Ncols () == colLength); // precondition
      else
	b.ReSize (rowLength, colLength);

      for (int i = 0; i < a.Nrows (); ++i)
	for (int j = 0; j < a.Ncols (); ++j)
	  {
#ifdef BLAS
	    int bcols = b.Ncols();
	    double* bptr = b.Store() + bcols * (indexRows[i] - 1) + (indexCols[j] - 1);
	    Matrix* aij = a(i, j);
	    double* aptr = aij->Store();
	    int nrows = aij->Nrows();
	    int ncols = aij->Ncols();
	    for (int r = 0; r < nrows; ++r)
	      {
		DCOPY(ncols, aptr, 1, bptr, 1);
		aptr += ncols;
		bptr += bcols;
	      }
#else
	    b.SubMatrix (indexRows [i], indexRows [i] + a (i,j)->Nrows () - 1, indexCols [j], indexCols [j] + a (i,j)->Ncols () - 1) = *(a (i,j));
#endif
	  }
    }
  catch (Exception)
    {
      pout << Exception::what () << endl;
      abort ();
    }
}

  
void SpinAdapted::MatrixRotate (const Matrix& a, const Matrix& b, const Matrix& c, Matrix& d)
{
  try
    {
      assert (d.Nrows () == a.Ncols () && d.Ncols () == c.Ncols ());
#ifdef BLAS
      Matrix work (b.Nrows (), c.Ncols ());
      Clear (work);
      MatrixMultiply (b, 'n', c, 'n', work, 1.);
      MatrixMultiply (a, 't', work, 'n', d, 1.);
#else
      d = a.t () * b * c;
#endif
    }
  catch (Exception)
    {
      pout << Exception::what () << endl;
      abort ();
    }
}

void SpinAdapted::Save(const Matrix& a, std::ofstream &ofs)
{
  boost::archive::binary_oarchive save_mat(ofs);
  save_mat << a;
}

void SpinAdapted::Load(Matrix& a, std::ifstream &ifs)
{
  boost::archive::binary_iarchive load_mat(ifs);
  load_mat >> a;
}

void SpinAdapted::DebugPrint (vector<int>& v)
{
  for (int i = 0; i < v.size(); ++i)
    pout << v[i] << endl;
}
void SpinAdapted::DebugPrint (vector<double>& v)
{
  for (int i = 0; i < v.size(); ++i)
    pout << v[i] << endl;
}

