//
// Implementation of CLAPACK:
//
// This is a part of BTAS library to provide C/C++ interface of LAPACK calls.
// Note that all functions work with Column-Major array with respect to FORTRAN-style,
// since C-style Row-Major array is not supported in LAPACK library.
//
// Wrapper functions for Row-Major array are provided in BTAS dense array layer.
//

#include <stdlib.h>
#include <lapack/clapack.h>
#include <lapack/lapack_fort_predecl.h>

// Triangular fragmentation

//! LU decomposition
int clapack_dgetrf
(int m, int n, double* a, int lda, int* ipiv) {
  int info;
  dgetrf_(&m, &n, a, &lda, ipiv, &info);
  return info;
}
//! Cholesky decomposition
int clapack_dpotrf
(CLAPACK_UPLO Uplo, int n, double* a, int lda) {
  char FC_UPLO = (Uplo == ClapackUpper) ? 'U' : 'L';
  int info;
  dpotrf_(&FC_UPLO, &n, a, &lda, &info);
  return info;
}
//! Bunch-Kaufman decomposition
int clapack_dsytrf
(CLAPACK_UPLO Uplo, int n, double* a, int lda, int* ipiv) {
  char FC_UPLO = (Uplo == ClapackUpper) ? 'U' : 'L';
  int info, lwork = -1;
  double lwork_opt[1];
  // calc. opt. length of workspace
  dsytrf_(&FC_UPLO, &n, a, &lda, ipiv, lwork_opt, &lwork, &info);
  lwork = (int) lwork_opt[0];
  double* work = (double*) malloc(lwork*sizeof(double));
  // call dsytrf
  dsytrf_(&FC_UPLO, &n, a, &lda, ipiv, work,      &lwork, &info);
  free(work);
  return info;
}

//! Solving linear equation from triangular fragmented matrix
int clapack_dgetrs
(CLAPACK_TRANSPOSE TransA, int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb) {
  char FC_TRANSA = (TransA == ClapackNoTrans) ? 'N' : 'T';
  int info;
  dgetrs_(&FC_TRANSA, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
  return info;
}
//! Solving linear equation from triangular fragmented matrix (symmetric positive definite)
int clapack_dpotrs
(CLAPACK_UPLO Uplo, int n, int nrhs, double* a, int lda, double* b, int ldb) {
  char FC_UPLO = (Uplo == ClapackUpper) ? 'U' : 'L';
  int info;
  dpotrs_(&FC_UPLO, &n, &nrhs, a, &lda, b, &ldb, &info);
  return info;
}
//! Solving linear equation from triangular fragmented matrix (symmetric)
int clapack_dsytrs
(CLAPACK_UPLO Uplo, int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb) {
  char FC_UPLO = (Uplo == ClapackUpper) ? 'U' : 'L';
  int info;
  dsytrs_(&FC_UPLO, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
  return info;
}

//! Solving linear equation (involving dgetrf)
int clapack_dgesv
(int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb) {
  int info = clapack_dgetrf(n, n, a, lda, ipiv);
  if(info == 0)
      info = clapack_dgetrs(ClapackNoTrans, n, nrhs, a, lda, ipiv, b, ldb);
  return info;
}
//! Solving linear equation (involving dpotrf)
int clapack_dposv
(CLAPACK_UPLO Uplo, int n, int nrhs, double* a, int lda, double* b, int ldb) {
  int info = clapack_dpotrf(Uplo, n, a, lda);
  if(info == 0)
      info = clapack_dpotrs(Uplo, n, nrhs, a, lda, b, ldb);
  return info;
}
//! Solving linear equation (involving dsytrf)
int clapack_dsysv
(CLAPACK_UPLO Uplo, int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb) {
  int info = clapack_dsytrf(Uplo, n, a, lda, ipiv);
  if(info == 0)
      info = clapack_dsytrs(Uplo, n, nrhs, a, lda, ipiv, b, ldb);
  return info;
}

//! Solve symmetric eigenvalue problem for double precision real matrix
int clapack_dsyev
(CLAPACK_CALCVECTOR Jobz, CLAPACK_UPLO Uplo, int n, double* a, int lda, double* w) {
  // parse job types
  char FC_JOBZ = (Jobz == ClapackNoCalcVector) ? 'N' : 'V';
  char FC_UPLO = (Uplo == ClapackUpper) ? 'U' : 'L';
  int info, lwork = -1;
  double lwork_opt[1];
  // calc. opt. length of workspace
  dsyev_(&FC_JOBZ, &FC_UPLO, &n, a, &lda, w, lwork_opt, &lwork, &info);
  lwork = (int) lwork_opt[0];
  double* work = (double*) malloc(lwork*sizeof(double));
  // call dsyev
  dsyev_(&FC_JOBZ, &FC_UPLO, &n, a, &lda, w, work,      &lwork, &info);
  // free workspace and return
  free(work);
  return info;
}

//! Solve non-hermitian eigenvalue problem for double precision real matrix
int clapack_dgeev
(CLAPACK_CALCVECTOR Jobl, CLAPACK_CALCVECTOR Jobr,
 int n, double* a, int lda, double* wr, double* wi, double* vl, int ldvl, double* vr, int ldvr) {
  char FC_JOBL = (Jobl == ClapackNoCalcVector) ? 'N' : 'V';
  char FC_JOBR = (Jobr == ClapackNoCalcVector) ? 'N' : 'V';
  int info, lwork = -1;
  double lwork_opt[1];
  // calc. opt. length of workspace
  dgeev_(&FC_JOBL, &FC_JOBR, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, lwork_opt, &lwork, &info);
  lwork = (int) lwork_opt[0];
  double* work = (double*) malloc(lwork*sizeof(double));
  // call dgeev
  dgeev_(&FC_JOBL, &FC_JOBR, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work,      &lwork, &info);
  // free workspace and return
  free(work);
  return info;
}

//! Solve generalized eigenvalue problem for double precision real symmetric matrix pair
/*! itype = 1: A x = w B x
 *          2: A B x = w x
 *          3: B A x = w x
 */
int clapack_dsygv
(int itype, CLAPACK_CALCVECTOR Jobz, CLAPACK_UPLO Uplo, int n, double* a, int lda, double* b, int ldb, double* w) {
  // parse job types
  char FC_JOBZ = (Jobz == ClapackNoCalcVector) ? 'N' : 'V';
  char FC_UPLO = (Uplo == ClapackUpper) ? 'U' : 'L';
  int info, lwork = -1;
  double lwork_opt[1];
  // calc. opt. length of workspace
  dsygv_(&itype, &FC_JOBZ, &FC_UPLO, &n, a, &lda, b, &ldb, w, lwork_opt, &lwork, &info);
  lwork = (int) lwork_opt[0];
  double* work = (double*) malloc(lwork*sizeof(double));
  // call dsygv
  dsygv_(&itype, &FC_JOBZ, &FC_UPLO, &n, a, &lda, b, &ldb, w, work,      &lwork, &info);
  // free workspace and return
  free(work);
  return info;
}

//! Solve generalized eigenvalue problem for double precision real non-symmetric matrix pair
int clapack_dggev
(CLAPACK_CALCVECTOR Jobl, CLAPACK_CALCVECTOR Jobr, int n, double* a, int lda, double* b, int ldb,
 double* alphar, double* alphai, double* beta, double* vl, int ldvl, double* vr, int ldvr) {
  // parse job types
  char FC_JOBL = (Jobl == ClapackNoCalcVector) ? 'N' : 'V';
  char FC_JOBR = (Jobr == ClapackNoCalcVector) ? 'N' : 'V';
  int info, lwork = -1;
  double lwork_opt[1];
  // calc. opt. length of workspace
  dggev_(&FC_JOBL, &FC_JOBR, &n, a, &lda, b, &ldb, alphar, alphai, beta, vl, &ldvl, vr, &ldvr, lwork_opt, &lwork, &info);
  lwork = (int) lwork_opt[0];
  double *work = (double*) malloc(lwork*sizeof(double));
  // call dggev
  dggev_(&FC_JOBL, &FC_JOBR, &n, a, &lda, b, &ldb, alphar, alphai, beta, vl, &ldvl, vr, &ldvr, work,      &lwork, &info);
  // free workspace and return
  free(work);
  return info;
}

//! Singularvalue decomposition for double precision real matrix
int clapack_dgesvd
(CLAPACK_CALCVECTOR Jobu, CLAPACK_CALCVECTOR Jobvt,
 int m, int n, double* a, int lda, double* s, double* u, int ldu, double* vt, int ldvt) {
  // parse job types
  char FC_JOBU;
  if(Jobu  == ClapackCalcVector)     FC_JOBU  = 'A';
  if(Jobu  == ClapackCalcThinVector) FC_JOBU  = 'S';
  if(Jobu  == ClapackNoCalcVector)   FC_JOBU  = 'N';
  char FC_JOBVT;
  if(Jobvt == ClapackCalcVector)     FC_JOBVT = 'A';
  if(Jobvt == ClapackCalcThinVector) FC_JOBVT = 'S';
  if(Jobvt == ClapackNoCalcVector)   FC_JOBVT = 'N';
  int info, lwork = -1;
  double lwork_opt[1];
  // calc. opt. length of workspace
  dgesvd_(&FC_JOBU, &FC_JOBVT, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, lwork_opt, &lwork, &info);
  lwork = (int) lwork_opt[0];
  double* work = (double*) malloc(lwork*sizeof(double));
  // call dgesvd
  dgesvd_(&FC_JOBU, &FC_JOBVT, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work,      &lwork, &info);
  // free workspace and return
  free(work);
  return info;
}

//! Solve hermitian eigenvalue problem (HEP) for double precision complex matrix
int clapack_zheev
(CLAPACK_CALCVECTOR Jobz, CLAPACK_UPLO Uplo, int n, void* a, int lda, double* w) {
  // parse job types
  char FC_JOBZ = (Jobz == ClapackNoCalcVector) ? 'N' : 'V';
  char FC_UPLO = (Uplo == ClapackUpper) ? 'U' : 'L';
  int info, lwork = -1;
  double lwork_opt[2];
  // calc. opt. length of workspace
  double* rwork = (double*) malloc(2*n*sizeof(double));
  zheev_(&FC_JOBZ, &FC_UPLO, &n, (FC_COMPLEX_16*) a, &lda, w, (FC_COMPLEX_16*) lwork_opt, &lwork, rwork, &info);
  lwork = (int) lwork_opt[0];
  FC_COMPLEX_16* work = (FC_COMPLEX_16*) malloc(lwork*sizeof(FC_COMPLEX_16)); // this assumes 16 bytes complex type
  // call zgeev
  zheev_(&FC_JOBZ, &FC_UPLO, &n, (FC_COMPLEX_16*) a, &lda, w, work,      &lwork, rwork, &info);
  // free workspace and return
  free( work);
  free(rwork);
  return info;
}

//! Singularvalue decomposition for double precision complex matrix
int clapack_zgesvd
(CLAPACK_CALCVECTOR Jobu, CLAPACK_CALCVECTOR Jobvt, int m, int n, void* a, int lda, double* s, void* u, int ldu, void* vt, int ldvt) {
  // parse job types
  char FC_JOBU;
  if(Jobu  == ClapackCalcVector)     FC_JOBU  = 'A';
  if(Jobu  == ClapackCalcThinVector) FC_JOBU  = 'S';
  if(Jobu  == ClapackNoCalcVector)   FC_JOBU  = 'N';
  char FC_JOBVT;
  if(Jobvt == ClapackCalcVector)     FC_JOBVT = 'A';
  if(Jobvt == ClapackCalcThinVector) FC_JOBVT = 'S';
  if(Jobvt == ClapackNoCalcVector)   FC_JOBVT = 'N';
  int info, lwork = -1;
  double lwork_opt[2];
  // calc. opt. length of workspace
  int rsize = (m < n) ? m : n;
  double *rwork = (double*) malloc(5*rsize*sizeof(double));
  zgesvd_(&FC_JOBU, &FC_JOBVT, &m, &n, (FC_COMPLEX_16*) a, &lda, s, (FC_COMPLEX_16*) u, &ldu, (FC_COMPLEX_16*) vt, &ldvt, (FC_COMPLEX_16*) lwork_opt, &lwork, rwork, &info);
  lwork = (int) lwork_opt[0];
  FC_COMPLEX_16* work = (FC_COMPLEX_16*) malloc(lwork*sizeof(FC_COMPLEX_16)); // this assumes 16 bytes complex type
  // call zgesvd
  zgesvd_(&FC_JOBU, &FC_JOBVT, &m, &n, (FC_COMPLEX_16*) a, &lda, s, (FC_COMPLEX_16*) u, &ldu, (FC_COMPLEX_16*) vt, &ldvt, work,      &lwork, rwork, &info);
  // free workspace and return
  free( work);
  free(rwork);
  return info;
}
