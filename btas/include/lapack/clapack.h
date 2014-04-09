//
//! \file clapack.h
//! \brief C/C++ interface to lapack calls like clapack_xxx
/*!
 *  This is a part of BTAS library to provide C/C++ interface of LAPACK calls.
 *  Note that all functions work with Column-Major array with respect to FORTRAN-style,
 *  since C-style Row-Major array is not supported in LAPACK library.
 *
 *  Wrapper functions for Row-Major array are provided in BTAS dense array layer.
 */

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// LAPACK collection for double precision real
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// clapack_dgetrf : LU decomposition
// clapack_dpotrf : Cholesky decomposition
// clapack_dsytrf : Bunch-Kaufman decomposition

// clapack_dgetrs : solving linear equation from triangular fragmented matrix
// clapack_dpotrs : (symmetric positive definite)
// clapack_dsytrs : (symmetric)

// clapack_dgesv  : solving linear equation (involving dgetrf)
// clapack_dposv  : (symmetric positive definite)
// clapack_dsysv  : (symmetric)

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// LAPACK collection for double precision complex
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// clapack_dsyev  : eigenvalue decomposition for real-symmetric matrix
// clapack_dgesvd : singularvalue decomposition for real general matrix

// clapack_zheev  : eigenvalue decomposition for complex-hermitian matrix
// clapack_zgesvd : singularvalue decomposition for complex general matrix

#ifndef _CLAPACK_H
#define _CLAPACK_H 1

//! Transposition type used in CLAPACK
enum CLAPACK_TRANSPOSE  { ClapackNoTrans = 111, ClapackTrans = 112, ClapackConjTrans = 113 };
//! Consider upper or lower triangular elements in symmetric matrix
enum CLAPACK_UPLO       { ClapackUpper = 121, ClapackLower = 122 };
//! Calculate eigenvectors
enum CLAPACK_CALCVECTOR { ClapackNoCalcVector = 201, ClapackCalcVector = 202, ClapackCalcThinVector = 203 };

//#ifdef _HAS_INTEL_MKL

//#include <mkl_lapack.h>

//typedef MKL_Complex8  FC_COMPLEX_08;
//typedef MKL_Complex16 FC_COMPLEX_16;

//#else // _HAS_INTEL_MKL

struct FC_COMPLEX_08 {
  float real;
  float imag;
};

struct FC_COMPLEX_16 {
  double real;
  double imag;
};

//#endif // _HAS_INTEL_MKL

#ifdef __cplusplus
extern "C" {
#endif

//####################################################################################################
// LAPACK collection for double precision real number
//####################################################################################################

//! LU decomposition
int clapack_dgetrf(int m, int n, double* a, int lda, int* ipiv);
//! Cholesky decomposition
int clapack_dpotrf(CLAPACK_UPLO Uplo, int n, double* a, int lda, int* ipiv);
//! Bunch-Kaufman decomposition
int clapack_dsytrf(CLAPACK_UPLO Uplo, int n, double* a, int lda, int* ipiv);

//! Solving linear equation from triangular fragmented matrix
int clapack_dgetrs(CLAPACK_TRANSPOSE TransA, int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb);
//! Solving linear equation from triangular fragmented matrix (symmetric positive definite)
int clapack_dpotrs(CLAPACK_UPLO Uplo, int n, int nrhs, double* a, int lda, double* b, int ldb);
//! Solving linear equation from triangular fragmented matrix (symmetric)
int clapack_dsytrs(CLAPACK_UPLO Uplo, int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb);

//! Solving linear equation (involving dgetrf)
int clapack_dgesv (int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb);
//! Solving linear equation (symmetric positive definite)
int clapack_dposv (CLAPACK_UPLO Uplo, int n, int nrhs, double* a, int lda, double* b, int ldb);
//! Solving linear equation (symmetric)
int clapack_dsysv (CLAPACK_UPLO Uplo, int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb);

//! eigenvalue decomposition for real-symmetric matrix
int clapack_dsyev (CLAPACK_CALCVECTOR Jobz, CLAPACK_UPLO Uplo, int n, double* a, int lda, double* w);
//! eigenvalue decomposition for non-hermitian matrix
int clapack_dgeev (CLAPACK_CALCVECTOR Jobl, CLAPACK_CALCVECTOR Jobr, int n, double* a, int lda, double* wr, double* wi, double* vl, int ldvl, double* vr, int ldvr);
//! generalized eigenvalue decomposition for real-symmetric matrix
int clapack_dsygv (int itype, CLAPACK_CALCVECTOR Jobz, CLAPACK_UPLO Uplo, int n, double* a, int lda, double* b, int ldb, double* w);
//! generalized eigenvalue decomposition for non-hermitian matrix pair
int clapack_dggev (CLAPACK_CALCVECTOR Jobl, CLAPACK_CALCVECTOR Jobr, int n, double* a, int lda, double* b, int ldb, double* alphar, double* alphai, double* beta, double* vl, int ldvl, double* vr, int ldvr);
//! singularvalue decomposition for real general matrix
int clapack_dgesvd(CLAPACK_CALCVECTOR Jobu, CLAPACK_CALCVECTOR Jobvt, int m, int n, double* a, int lda, double* s, double* u, int ldu, double* vt, int ldvt);

//####################################################################################################
// LAPACK collection for double precision complex number
//####################################################################################################

//! eigenvalue decomposition for complex-hermitian matrix
int clapack_zheev (CLAPACK_CALCVECTOR Jobz, CLAPACK_UPLO Uplo, int n, void* a, int lda, double* w);
//! singularvalue decomposition for complex general matrix
int clapack_zgesvd(CLAPACK_CALCVECTOR Jobu, CLAPACK_CALCVECTOR Jobvt, int m, int n, void* a, int lda, double* s, void* u, int ldu, void* vt, int ldvt);

#ifdef __cplusplus
}
#endif

#endif // _CLAPACK_H
