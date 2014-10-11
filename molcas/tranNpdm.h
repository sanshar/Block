#ifndef __TRANSFORM_NPDM_ARRAY_H
#define __TRANSFORM_NPDM_ARRAY_H

extern "C" {

void block_tran1pdm_ (const int* N, double* X, const int* iRoot, const int* jRoot);
void block_tran2pdm_ (const int* N, double* X, const int* iRoot, const int* jRoot);
void block_tran3pdm_ (const int* N, double* X, const int* iRoot, const int* jRoot);

} // extern "C"

void tran1pdm (int N, double* X, int iRoot, int jRoot);
void tran2pdm (int N, double* X, int iRoot, int jRoot);
void tran3pdm (int N, double* X, int iRoot, int jRoot);

void SF_tran2pdm (int N, double* X, int iRoot, int jRoot);
void SF_tran3pdm (int N, double* X, int iRoot, int jRoot);

#endif // __TRANSFORM_NPDM_ARRAY_H
