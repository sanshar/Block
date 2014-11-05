#ifndef __BLOCK_CALL_DMRG_H
#define __BLOCK_CALL_DMRG_H

#include "molcas_types.h"

extern "C" {

void block_calldmrg_ (
      const FORTINT* Restart, // 0 to full calc., 1 to restart
      const FORTINT* N_roots,
      const FORTINT* N_act,
      const FORTINT* N_elec,
      const FORTINT* M_s,
      const char* Sym,
      const FORTINT* iSym,
      const FORTINT* OrbSym,
      const double* E_core,
      const double* h0,
      const double* tuvx,
      const FORTINT* M_state,
      const FORTINT* N_pdm,
      const double* T_sweep,
      const double* T_noise,
            double* E_sweep);

} // extern "C"

void block_calldmrg (
      const FORTINT& Restart,
      const FORTINT& N_roots,
      const FORTINT& N_act,
      const FORTINT& N_elec,
      const FORTINT& M_s,
      const char* Sym,
      const FORTINT& iSym,
      const FORTINT* OrbSym,
      const double& E_core,
      const double* h0,
      const double* tuvx,
      const FORTINT& M_state,
      const FORTINT& N_pdm,
      const double& T_sweep,
      const double& T_noise,
            double* E_sweep);

#endif // __BLOCK_CALL_DMRG_H
