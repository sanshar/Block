#ifndef __BLOCK_CALL_DMRG_H
#define __BLOCK_CALL_DMRG_H

extern "C" {

void block_calldmrg_ (
      const int* Restart, // 0 to full calc., 1 to restart
//    const char* Sym,
      const int* N_roots,
      const int* N_act,
      const int* N_elec,
      const int* M_s,
      const double* E_core,
      const double* h0,
      const double* tuvx,
      const int* M_state,
      const int* N_pdm,
      const double* T_sweep,
      const double* T_noise,
            double* E_sweep);

} // extern "C"

void block_calldmrg (
      const int& Restart,
//    const char* Sym,
      const int& N_roots,
      const int& N_act,
      const int& N_elec,
      const int& M_s,
      const double& E_core,
      const double* h0,
      const double* tuvx,
      const int& M_state,
      const int& N_pdm,
      const double& T_sweep,
      const double& T_noise,
            double* E_sweep);

#endif // __BLOCK_CALL_DMRG_H
