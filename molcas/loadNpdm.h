#ifndef __LOAD_NPDM_ARRAY_H
#define __LOAD_NPDM_ARRAY_H

extern "C" {

void load1pdm_   (const int* N, double* V);

void load2pdm_   (const int* N, double* V);
void load2pdm2f_ (const int* N, double* V, const int* iP, const int* iQ);

void load3pdm_   (const int* N, double* V);
void load3pdm2f_ (const int* N, double* V, const int* iP, const int* iQ);
void load3pdm4f_ (const int* N, double* V, const int* iP, const int* iQ, const int* jP, const int* jQ);

} // extern "C"

/// load V(N,N) from disk
void load1pdm   (int N, double* V);

/// load V(N,N,N,N) from disk
void load2pdm   (int N, double* V);
/// load V(N,N,iP,iQ) from disk
void load2pdm2f (int N, double* V, int iP, int iQ);

/// load V(N,N,N,N,N,N) from disk
void load3pdm   (int N, double* V);
/// load V(N,N,N,N,iP,iQ) from disk
void load3pdm2f (int N, double* V, int iP, int iQ);
/// load V(N,N,iP,iQ,jP,jQ) from disk
void load3pdm4f (int N, double* V, int iP, int iQ, int jP, int jQ);

#endif // __LOAD_NPDM_ARRAY_H
