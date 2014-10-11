#ifndef __LOAD_NPDM_ARRAY_H
#define __LOAD_NPDM_ARRAY_H

extern "C" {

void block_load1pdm_   (const int* N, double* V, const int* iRoot, const int* jRoot);

void block_load2pdm_   (const int* N, double* V, const int* iRoot, const int* jRoot);
void block_load2pdm2f_ (const int* N, double* V, const int* iRoot, const int* jRoot, const int* iP, const int* iQ);

void block_load3pdm_   (const int* N, double* V, const int* iRoot, const int* jRoot);
void block_load3pdm2f_ (const int* N, double* V, const int* iRoot, const int* jRoot, const int* iP, const int* iQ);
void block_load3pdm4f_ (const int* N, double* V, const int* iRoot, const int* jRoot, const int* iP, const int* iQ, const int* jP, const int* jQ);

} // extern "C"

/// load V(N,N) from disk
void load1pdm   (int N, double* V, int iRoot, int jRoot);

/// load V(N,N,N,N) from disk
void load2pdm   (int N, double* V, int iRoot, int jRoot);
/// load V(N,N,iP,iQ) from disk
void load2pdm2f (int N, double* V, int iRoot, int jRoot, int iP, int iQ);

/// load V(N,N,N,N,N,N) from disk
void load3pdm   (int N, double* V, int iRoot, int jRoot);
/// load V(N,N,N,N,iP,iQ) from disk
void load3pdm2f (int N, double* V, int iRoot, int jRoot, int iP, int iQ);
/// load V(N,N,iP,iQ,jP,jQ) from disk
void load3pdm4f (int N, double* V, int iRoot, int jRoot, int iP, int iQ, int jP, int jQ);

#endif // __LOAD_NPDM_ARRAY_H
