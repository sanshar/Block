#ifndef __TRANSFORM_NPDM_ARRAY_H
#define __TRANSFORM_NPDM_ARRAY_H

extern "C" {

void tran1pdm_ (int* N, double* X);
void tran2pdm_ (int* N, double* X);
void tran3pdm_ (int* N, double* X);

} // extern "C"

void tran1pdm (int N, double* X);
void tran2pdm (int N, double* X);
void tran3pdm (int N, double* X);

void SF_tran2pdm (int N, double* X);
void SF_tran3pdm (int N, double* X);

#endif // __TRANSFORM_NPDM_ARRAY_H
