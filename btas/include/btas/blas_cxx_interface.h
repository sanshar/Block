#ifndef _BLAS_CXX_INTERFACE_H
#define _BLAS_CXX_INTERFACE_H 1

#ifdef _HAS_CBLAS

#ifdef _HAS_INTEL_MKL
#include <mkl_cblas.h>
#else  // _HAS_INTEL_MKL
#include <cblas.h>


#endif // _HAS_INTEL_MKL


namespace btas {

typedef CBLAS_TRANSPOSE BTAS_TRANSPOSE;
typedef CBLAS_ORDER     BTAS_ORDER;
typedef CBLAS_UPLO      BTAS_UPLO;
typedef CBLAS_DIAG      BTAS_DIAG;
typedef CBLAS_SIDE      BTAS_SIDE;

extern const BTAS_TRANSPOSE Trans    ; // = CblasTrans;
extern const BTAS_TRANSPOSE NoTrans  ; // = CblasNoTrans;
extern const BTAS_TRANSPOSE ConjTrans; // = CblasConjTrans;

extern const BTAS_ORDER     RowMajor ; // = CblasRowMajor;
extern const BTAS_ORDER     ColMajor ; // = CblasColMajor;

extern const BTAS_UPLO      Upper    ; // = CblasUpper;
extern const BTAS_UPLO      Lower    ; // = CblasLower;

extern const BTAS_DIAG      NonUnit  ; // = CblasNonUnit;
extern const BTAS_DIAG      Unit     ; // = CblasUnit;

extern const BTAS_SIDE      Left     ; // = CblasLeft;
extern const BTAS_SIDE      Right    ; // = CblasRight;

};

#else  // _HAS_CBLAS

namespace btas {

enum BTAS_ORDER     { RowMajor =101, ColMajor =102 };
enum BTAS_TRANSPOSE { NoTrans  =111, Trans    =112, ConjTrans=113 };
enum BTAS_UPLO      { Upper    =121, Lower    =122 };
enum BTAS_DIAG      { NonUnit  =131, Unit     =132 };
enum BTAS_SIDE      { Left     =141, Right    =142 };

// Write BLAS wrappers below

};

#endif // _HAS_CBLAS

#endif // _BLAS_CXX_INTERFACE_H
