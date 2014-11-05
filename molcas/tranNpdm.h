#ifndef __TRANSFORM_NPDM_ARRAY_H
#define __TRANSFORM_NPDM_ARRAY_H

#include "molcas_types.h"

extern "C" {

void block_tran1pdm_ (const FORTINT* N, double* X, const FORTINT* iRoot, const FORTINT* jRoot);
void block_tran2pdm_ (const FORTINT* N, double* X, const FORTINT* iRoot, const FORTINT* jRoot);
void block_tran3pdm_ (const FORTINT* N, double* X, const FORTINT* iRoot, const FORTINT* jRoot);

} // extern "C"

void tran1pdm (FORTINT N, double* X, FORTINT iRoot, FORTINT jRoot);
void tran2pdm (FORTINT N, double* X, FORTINT iRoot, FORTINT jRoot);
void tran3pdm (FORTINT N, double* X, FORTINT iRoot, FORTINT jRoot);

void SF_tran2pdm (FORTINT N, double* X, FORTINT iRoot, FORTINT jRoot);
void SF_tran3pdm (FORTINT N, double* X, FORTINT iRoot, FORTINT jRoot);

#endif // __TRANSFORM_NPDM_ARRAY_H
