#ifndef __MOLPRO_FCIDUMP_H
#define __MOLPRO_FCIDUMP_H

#include "molcas_types.h"

extern "C" {

/// Fortran wrapper
void molpro_fcidump_ (
      const FORTINT* N_act,
      const FORTINT* N_elec,
      const FORTINT* M_s,
      const FORTINT* iSym,
      const FORTINT* OrbSym,
      const double* E_core,
      const double* h0,
      const double* tuvx);

} // extern "C"

/// Dump 1- and 2-particle integrals into formatted file
void molpro_fcidump (
      const FORTINT& N_act,
      const FORTINT& N_elec,
      const FORTINT& M_s,
      const FORTINT& iSym,
      const FORTINT* OrbSym,
      const double& E_core,
      const double* h0,
      const double* tuvx);

#endif // __MOLPRO_FCIDUMP_H
