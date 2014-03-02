#ifndef _LAPACK_FORT_PREDECL_H
#define _LAPACK_FORT_PREDECL_H 1


#ifdef __cplusplus
extern "C" {
#endif

void dgetrf_(int* m, int* n, double* a, int* lda, int* ipiv, int* info);
void dpotrf_(char* uplo, int* n, double* a, int* lda, int* info);
void dsytrf_(char* uplo, int* n, double* a, int* lda, int* ipiv, double* work, int* lwork, int* info);
void dgetrs_(char* transa, int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);
void dpotrs_(char* uplo, int* n, int* nrhs, double* a, int* lda, double* b, int* ldb, int* info);
void dsytrs_(char* uplo, int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);
void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info);
void dgeev_(char* jobl, char* jobr, int* n, double* a, int* lda, double* wr, double* wi, double* vl, int* ldvl, double* vr, int* ldvr, double* work, int* lwork, int* info);
void dsygv_(int* itype, char* jobz, char* uplo, int* n, double* a, int* lda, double* b, int* ldb, double* w, double* work, int* lwork, int* info);
void dggev_(char* jobl, char* jobr, int* n, double* a, int* lda, double* b, int* ldb, double* alphar, double* alphai, double* beta, double* vl, int* ldvl, double* vr, int* ldvr, double* work, int* lwork, int* info);
void dgesvd_(char* jobu, char* jobvt, int* m, int* n, double* a, int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt, double* work, int* lwork, int* info);

void zheev_(char* jobz, char* uplo, int* n, FC_COMPLEX_16* a, int* lda, double* w, FC_COMPLEX_16* work, int* lwork, double* rwork, int* info);
void zgesvd_(char* jobu, char* jobvt, int* m, int* n, FC_COMPLEX_16* a, int* lda, double* s, FC_COMPLEX_16* u, int* ldu, FC_COMPLEX_16* vt, int* ldvt, FC_COMPLEX_16* work, int* lwork, double* rwork, int* info);

#ifdef __cplusplus
}
#endif



#endif // _LAPACK_FORT_PREDECL_H
