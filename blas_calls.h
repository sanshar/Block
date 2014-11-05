/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef BLAS_CALLS_HEADER
#define BLAS_CALLS_HEADER
// This may look like C code, but it is really -*- C++ -*-
//AIK, 1998

//This file contains C-wrapped BLAS/Lapack functions calls. No changes in
//function cals were done: user has to care himself to match C and Fortran
//conventions. 

//For crossplatforming out if_def's HERE, in this file only.
//Currently sintaxis of Fortran calls corresponds to AIX/IBM compilers.
//Linking blas/lapack libraries: IBM: blas and /usr/local/lib/lapack.a 
//                               DEC: dxml
//                               SGI: sgimath

#ifndef _FORTINT_64
#define FORTINT int
#else
#define FORTINT long
#endif
 
#ifdef AIX
extern "C"
{
  void dgesvd(char* JOBU, char* JOBVT, int* M, int* N, double* A, int* LDA, double* S, double* U, int* LDU, double* VT, int* LDVT, double* WORK, int* LWORK,
	 int* INFO);
  void daxpy(int *ntot, double *coeff, double *copy_from, int *inc1,
	     double *copy_to, int *inc2);
  void dcopy(int *ntot,double *copy_from,int *inc1, double *copy_to,int *inc2);
  double ddot(int *ntot, double *x, int *incx, double *y, int *incy);
  void dgemm(char *transa, char *transb, int *m, int *n, int *k, double *alpha,
	     double *A, int *lda, double *B, int *ldb, double *beta, double *C,
	     int *ldc);
  void dgemv(char *trans, int *m, int *n, double *alpha,
             double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);
  void dscal(int *size,double *coeff,double *matrix,int *inc);

  void saxpy(int *ntot, float *coeff, float *copy_from, int *inc1,
	     float *copy_to, int *inc2);
  void scopy(int *ntot, float *copy_from, int *inc1, float *copy_to, int *inc2);
  float sdot(int *ntot, float *x, int *incx, float *y, int *incy);
  void sgemm(char *transa, char *transb, int *m, int *n, int *k, float *alpha,
	     float *A, int *lda, float *B, int *ldb, float *beta, float *C,
	     int *ldc);
  void sscal(int *size, float *coeff, float *matrix,int *inc);
  void dstev(char* JOBZ,FORTINT* N,double* A,double* E, double* W, FORTINT* Wlen, double*WORK, FORTINT* INFO);
  void dsyev(char* JOBZ,char* UPLO,int* N,double* A,int* LDA,double* W,double*WORK,int*LWORK,int* INFO);
  void dgesv(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
}
#else //SGI, Linux
extern "C"
{
  void dgesvd_(char* JOBU, char* JOBVT, FORTINT* M, FORTINT* N, double* A, FORTINT* LDA, double* S, double* U, FORTINT* LDU, double* VT, FORTINT* LDVT, double* WORK, FORTINT* LWORK,
	 FORTINT* INFO);
  void daxpy_(FORTINT *ntot, double *coeff, double *copy_from, FORTINT *inc1,
	     double *copy_to, FORTINT *inc2);
  void dcopy_(FORTINT *ntot,double *copy_from,FORTINT *inc1, double *copy_to,FORTINT *inc2);
  double ddot_(FORTINT *ntot, double *x, FORTINT *incx, double *y, FORTINT *incy);
  void dgemm_(char *transa, char *transb, FORTINT *m, FORTINT *n, FORTINT *k, double *alpha,
	     double *A, FORTINT *lda, double *B, FORTINT *ldb, double *beta, double *C,
	     FORTINT *ldc);
  void dgemv_(char *trans, FORTINT *m, FORTINT *n, double *alpha,
              double *a, FORTINT *lda, double *x, FORTINT *incx, double *beta, double *y, FORTINT *incy);
  void dscal_(FORTINT *size, double *coeff, double *matrix,FORTINT *inc);

  void saxpy_(FORTINT *ntot, float *coeff, float *copy_from, FORTINT *inc1,
	     float *copy_to, FORTINT *inc2);
  void scopy_(FORTINT *ntot, float *copy_from, FORTINT *inc1, float *copy_to, FORTINT *inc2);
  float sdot_(FORTINT *ntot, float *x, FORTINT *incx, float *y, FORTINT *incy);
  void sgemm_(char *transa, char *transb, FORTINT *m, FORTINT *n, FORTINT *k, float *alpha,
	     float *A, FORTINT *lda, float *B, FORTINT *ldb, float *beta, float *C,
	     FORTINT *ldc);
  void sscal_(FORTINT *size, float *coeff, float *matrix,FORTINT *inc);
  void dstev_(char* JOBZ,FORTINT* N,double* A,double* E, double* W, FORTINT* Wlen, double*WORK, FORTINT* INFO);
  void dsyev_(char* JOBZ,char* UPLO,FORTINT* N,double* A,FORTINT* LDA,double* W,double*WORK,FORTINT*LWORK,FORTINT* INFO);
  void dgesv_(FORTINT *n, FORTINT *nrhs, double *a, FORTINT *lda, FORTINT *ipiv, double *b, FORTINT *ldb, FORTINT *info);
  int idamax_(FORTINT &n, double* d, FORTINT &indx);
  //int idamax_(int &n, double* d, int &indx);
}
#endif


/*
** 
** void DAXPY(int ntot, double coeff, double *copy_from, int inc1,
**            double *copy_to, int inc2);
** This function adds data from copy_from to copy_to
** 
** int ntot: length of data
** int inc1,inc2: increments for copy_from, copy_to
*/
inline void DAXPY(FORTINT ntot, double coeff, double *copy_from, FORTINT inc1,
	   double *copy_to, FORTINT inc2)
{
#ifdef AIX
  daxpy(&ntot,&coeff,copy_from,&inc1,copy_to,&inc2);
#else
  daxpy_(&ntot,&coeff,copy_from,&inc1,copy_to,&inc2);
#endif
}
inline void SAXPY(FORTINT ntot, float coeff, float *copy_from, FORTINT inc1,
	   float *copy_to, FORTINT inc2)
{
#ifdef AIX
  saxpy(&ntot,&coeff,copy_from,&inc1,copy_to,&inc2);
#else
  saxpy_(&ntot,&coeff,copy_from,&inc1,copy_to,&inc2);
#endif
}

/*
** 
** void DCOPY(int ntot, double *x, int incx, double *y, int *incy);
** This function copies x to y
** 
** int ntot: length of x,y
** int incx,incy: increments for x,y 
*/
inline void DCOPY(FORTINT ntot,double *copy_from,FORTINT inc1,double *copy_to,FORTINT inc2)
{
#ifdef AIX
  dcopy(&ntot,copy_from,&inc1,copy_to,&inc2);
#else
  dcopy_(&ntot,copy_from,&inc1,copy_to,&inc2);
#endif
}
inline void SCOPY(FORTINT ntot,float *copy_from,FORTINT inc1,float *copy_to,FORTINT inc2)
{
#ifdef AIX
  scopy(&ntot,copy_from,&inc1,copy_to,&inc2);
#else
  scopy_(&ntot,copy_from,&inc1,copy_to,&inc2);
#endif
}



/*
** 
** void DDOT(int ntot, double *x, int incx, double *y, int *incy);
**
** This function calculates dot product (x,y)
** 
** int ntot: length of x,y
** int incx,incy: increments for x,y 
*/
inline double DDOT(FORTINT ntot, double *x, FORTINT incx, double *y, FORTINT incy)
{
#ifdef AIX
  return ddot(&ntot,x,&incx,y,&incy);
#else
  return ddot_(&ntot,x,&incx,y,&incy);
#endif
}
inline float SDOT(FORTINT ntot, float *x, FORTINT incx, float *y, FORTINT incy)
{
#ifdef AIX
  return sdot(&ntot,x,&incx,y,&incy);
#else
  return sdot_(&ntot,x,&incx,y,&incy);
#endif
}

/*
** 
** void DGEMM(char transa, char transb, int m, int n, int k, double alpha,
**	      double *A, int lda, double *B, int ldb, double beta, double *C,
**            int ldc);
** This function calculates C(m,n)=alpha*(opT)A(m,k)*(opT)B(k,n)+ beta*C(m,n)
** 
** char transa:       On entry, specifies the form of (op)A used in the
**                    matrix multiplication:
**                    If transa = 'N' or 'n', (op)A = A
**                    If transa = 'T' or 't', (op)A = transp(A)
**                    If transa = 'R' or 'r', (op)A = conjugate(A)
**                    If transa = 'C' or 'c', (op)A = conjug_transp(A)
**                    On exit, transa is unchanged.
** char transb:       On entry, specifies the form of (op)B used in the
**                    matrix multiplication:
**                    If transb = 'N' or 'n', (op)B = B
**                    If transb = 'T' or 't', (op)B = transp(B)
**                    If transb = 'R' or 'r', (op)B = conjugate(B)
** int m:             On entry, the number of rows of the matrix (op)A and of
**                    the matrix C; m >= 0. On exit, m is unchanged.
** int n:             On entry, the number of columns of the matrix (op)B and
**                    of the matrix C; n >= 0. On exit, n is unchanged.
** int k:             On entry, the number of columns of the matrix (op)A and
**                    the number of rows of the matrix (op)B; k >= 0. On exit,
**		      k is unchanged.
** double alpha:      On entry, specifies the scalar alpha. On exit, alpha is
**                    unchanged.
** double *A:         On entry, a two-dimensional array A with dimensions lda
**                    by ka. For (op)A = A  or  conjugate(A), ka >= k and the 
**                    leading m by k portion of the array A contains the matrix
**                    A. For (op)A = transp(A) or conjug_transp(A), ka >= m
**                    and the leading k by m part of the array A contains the
**                    matrix A. On exit, a is unchanged.
** int lda:           On entry, the first dimension of array A.
**                    For (op)A = A  or conjugate(A), lda >= MAX(1,m).
**                    For (op)A=transp(A) or conjug_transp(A), lda >= MAX(1,k).
**                    On exit, lda is unchanged.
** double *B:         On entry, a two-dimensional array B with dimensions ldb
**                    by kb. For (op)B = B or conjugate(B), kb >= n and the
**                    leading k by n portion of the array contains the matrix
**                    B. For (op)B = transp(B) or conjug_transp(B), kb >= k and
**                    the leading n by k part of the array contains the matrix
**		      B. On exit, B is unchanged.
** int ldb:           On entry, the first dimension of array B.
**                    For (op)B = B or <conjugate(B), ldb >= MAX(1,k).
**                    For (op)B = transp(B) or conjug_transp(B), ldb >=
**                    MAX(1,n). On exit, ldb is unchanged.
** double *beta:      On entry, specifies the scalar beta. On exit, beta is
**                    unchanged.
** double C:          On entry, a two-dimensional array with the dimension
**                    ldc by at least n. On exit,  the leading  m by n part of
**                    array C is overwritten by the matrix alpha*(op)A*(op)B +
**                    beta*C.
** int ldc:           On entry, the first dimension  of array C; ldc >=MAX(1,m)
**                    On exit, ldc is unchanged.
**
*/
inline void DGEMM(char transa, char transb, FORTINT m, FORTINT n, FORTINT k, double alpha,
	   double *A, FORTINT lda, double *B, FORTINT ldb, double beta, double *C,
	   FORTINT ldc)
{
#ifdef AIX
  dgemm(&transa,&transb,&m,&n,&k,&alpha,A,&lda,B,&ldb,&beta,C,&ldc);
#else
  dgemm_(&transa,&transb,&m,&n,&k,&alpha,A,&lda,B,&ldb,&beta,C,&ldc);
#endif
}

inline void DGEMV(char trans, FORTINT m, FORTINT n, double alpha,
	   double *A, FORTINT lda, double *X, FORTINT incx, double beta, double *Y, FORTINT incy)
{
#ifdef AIX
  dgemv(&trans,&m,&n,&alpha,A,&lda,X,&incx,&beta,Y,&incy);
#else
  dgemv_(&trans,&m,&n,&alpha,A,&lda,X,&incx,&beta,Y,&incy);
#endif
}

inline void DSYEV(char JOBZ, char UPLO, FORTINT N, double* A, FORTINT LDA, double* W, double* WORK, FORTINT LWORK, FORTINT INFO )
{
#ifdef AIX
  dsyev(&JOBZ,&UPLO,&N,A,&LDA,W,WORK,&LWORK,&INFO);
#else
  dsyev_(&JOBZ,&UPLO,&N,A,&LDA,W,WORK,&LWORK,&INFO);
#endif
}

inline void DSTEV(char JOBZ, FORTINT N, double* D, double* E, double* vec, FORTINT LDA, double* W, FORTINT INFO )
{
#ifdef AIX
  dstev(&JOBZ,&N,D,E,vec,&LDA,W,&INFO);
#else
  dstev_(&JOBZ,&N,D,E,vec,&LDA,W,&INFO);
#endif
}

inline void GESV(FORTINT n, FORTINT nrhs, double* a, FORTINT lda, FORTINT* ipiv, double* b, FORTINT ldb, FORTINT info)
{
#ifdef AIX
  dgesv(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
#else
  dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
#endif
}


inline void SGEMM(char transa, char transb, FORTINT m, FORTINT n, FORTINT k, float alpha,
	   float *A, FORTINT lda, float *B, FORTINT ldb, float beta, float *C,
	   FORTINT ldc)
{
#ifdef AIX
  sgemm(&transa,&transb,&m,&n,&k,&alpha,A,&lda,B,&ldb,&beta,C,&ldc);
#else
  sgemm_(&transa,&transb,&m,&n,&k,&alpha,A,&lda,B,&ldb,&beta,C,&ldc);
#endif
}

// singular value decomposition
inline void DGESVD(char JOBU, char JOBVT, FORTINT M, FORTINT N, double* A, 
		   FORTINT LDA, double* S, double* U, FORTINT LDU, double* VT,
		   FORTINT LDVT,double* WORK, FORTINT LWORK, FORTINT INFO)
{
#ifdef AIX
  dgesvd(&JOBU, &JOBVT, &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK,
	 &INFO);
#else
  dgesvd_(&JOBU, &JOBVT, &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK,
	  &INFO);
#endif
}
/*
** 
** void DSCAL(int ntot, double coeff, double *data, int inc);
**
** This function scales ntot elements from array data with coeff and increment
** inc : data*=coeff
** 
** int ntot:      length of data to scale
** double *data:  data to scale
** int inc:       increments for data
*/
inline void DSCAL(FORTINT size, double coeff, double *data, FORTINT inc)
{
  if( 1.!=coeff )
#ifdef AIX
    dscal(&size,&coeff,data,&inc);
#else
  dscal_(&size,&coeff,data,&inc);
#endif
}
inline void SSCAL(FORTINT size, float coeff, float *data, FORTINT inc)
{
  if( 1.!=coeff )
#ifdef AIX
    sscal(&size,&coeff,data,&inc);
#else
  sscal_(&size,&coeff,data,&inc);
#endif
}


#endif
