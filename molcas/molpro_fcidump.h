#ifndef __MOLPRO_FCIDUMP_H
#define __MOLPRO_FCIDUMP_H

extern "C" {

/// Fortran wrapper
void molpro_fcidump_ (
      const int* N_act,
      const int* N_elec,
      const int* M_s,
      const double* E_core,
      const double* h0,
      const double* tuvx);

} // extern "C"

/// Dump 1- and 2-particle integrals into formatted file
void molpro_fcidump (
      const int& N_act,
      const int& N_elec,
      const int& M_s,
      const double& E_core,
      const double* h0,
      const double* tuvx);

#endif // __MOLPRO_FCIDUMP_H
