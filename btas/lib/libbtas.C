#include <btas/blas_cxx_interface.h>

namespace btas {

const BTAS_TRANSPOSE Trans     = CblasTrans;
const BTAS_TRANSPOSE NoTrans   = CblasNoTrans;
const BTAS_TRANSPOSE ConjTrans = CblasConjTrans;

const BTAS_ORDER     RowMajor  = CblasRowMajor;
const BTAS_ORDER     ColMajor  = CblasColMajor;

const BTAS_UPLO      Upper     = CblasUpper;
const BTAS_UPLO      Lower     = CblasLower;

const BTAS_DIAG      NonUnit   = CblasNonUnit;
const BTAS_DIAG      Unit      = CblasUnit;

const BTAS_SIDE      Left      = CblasLeft;
const BTAS_SIDE      Right     = CblasRight;

};
