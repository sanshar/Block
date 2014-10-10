#ifndef __BLOCK_DENSI_RASSCF_H
#define __BLOCK_DENSI_RASSCF_H

extern "C" {

void block_densi_rasscf_ (const int* iRoot, const int* N_act, const int* N_elec, double* D, double* DS, double* P, double* PA);

} // extern "C"

void block_densi_rasscf  (const int& iRoot, const int& N_act, const int& N_elec, double* D, double* DS, double* P, double* PA);

#endif // __BLOCK_DENSI_RASSCF_H
