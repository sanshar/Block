#ifndef __LOAD_NPDM_ARRAY_H
#define __LOAD_NPDM_ARRAY_H

#include "molcas_types.h"

extern "C" {

void block_load1pdm_   (const FORTINT* N, double* V, const FORTINT* iRoot, const FORTINT* jRoot);

void block_load2pdm_   (const FORTINT* N, double* V, const FORTINT* iRoot, const FORTINT* jRoot);
void block_load2pdm2f_ (const FORTINT* N, double* V, const FORTINT* iRoot, const FORTINT* jRoot, const FORTINT* iP, const FORTINT* iQ);

void block_load3pdm_   (const FORTINT* N, double* V, const FORTINT* iRoot, const FORTINT* jRoot);
void block_load3pdm2f_ (const FORTINT* N, double* V, const FORTINT* iRoot, const FORTINT* jRoot, const FORTINT* iP, const FORTINT* iQ);
void block_load3pdm4f_ (const FORTINT* N, double* V, const FORTINT* iRoot, const FORTINT* jRoot, const FORTINT* iP, const FORTINT* iQ, const FORTINT* jP, const FORTINT* jQ);

} // extern "C"

/// load V(N,N) from disk
void load1pdm   (FORTINT N, double* V, FORTINT iRoot, FORTINT jRoot);

/// load V(N,N,N,N) from disk
void load2pdm   (FORTINT N, double* V, FORTINT iRoot, FORTINT jRoot);
/// load V(N,N,iP,iQ) from disk
void load2pdm2f (FORTINT N, double* V, FORTINT iRoot, FORTINT jRoot, FORTINT iP, FORTINT iQ);

/// load V(N,N,N,N,N,N) from disk
void load3pdm   (FORTINT N, double* V, FORTINT iRoot, FORTINT jRoot);
/// load V(N,N,N,N,iP,iQ) from disk
void load3pdm2f (FORTINT N, double* V, FORTINT iRoot, FORTINT jRoot, FORTINT iP, FORTINT iQ);
/// load V(N,N,iP,iQ,jP,jQ) from disk
void load3pdm4f (FORTINT N, double* V, FORTINT iRoot, FORTINT jRoot, FORTINT iP, FORTINT iQ, FORTINT jP, FORTINT jQ);

#endif // __LOAD_NPDM_ARRAY_H
