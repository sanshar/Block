#include "molpro_fcidump.h"

//extern "C" {

/// Fortran wrapper
void molpro_fcidump_ (
      const FORTINT* N_act,
      const FORTINT* N_elec,
      const FORTINT* M_s,
      const FORTINT* iSym,
      const FORTINT* OrbSym,
      const double* E_core,
      const double* h0,
      const double* tuvx)
{
  molpro_fcidump(*N_act, *N_elec, *M_s, *iSym, OrbSym, *E_core, h0, tuvx);
}

//} // extern "C"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

/// Dump 2-el FORTINTegrals to formatted file
void molpro_fcidump (
      const FORTINT& N_act,
      const FORTINT& N_elec,
      const FORTINT& M_s,
      const FORTINT& iSym,
      const FORTINT* OrbSym,
      const double& E_core,
      const double* h0,
      const double* tuvx)
{
  using std::endl;
  using std::setw;
  using std::scientific;

//std::cout << "\t Dumping MOERI" << endl;

  std::ofstream fout("FCIDUMP");

  fout << " &FCI NORB=" << setw(3) << N_act << ",NELEC=" << setw(3) << N_elec << ",MS2=" << setw(2) << M_s << "," << endl;
  fout << "  ORBSYM="; for(int i = 0; i < N_act; ++i) fout << OrbSym[i] << ","; fout << endl;
  fout << "  ISYM=" << iSym << "," << endl;
  fout << " /" << endl;

  fout.precision(16);
  for(int i = 0; i < N_act; ++i) {
    for(int j = 0; j <= i; ++j) {
      int ij = i*(i+1)/2+j;

      for(int k = 0; k < N_act; ++k) {
        for(int l = 0; l <= k; ++l) {
          int kl = k*(k+1)/2+l;

          if(kl > ij) continue;

          int ijkl = ij*(ij+1)/2+kl;

          if(fabs(tuvx[ijkl]) < 1.0e-16) continue;

          fout << setw(24) << scientific << tuvx[ijkl]
               << setw( 4) << i+1
               << setw( 4) << j+1
               << setw( 4) << k+1
               << setw( 4) << l+1 << endl;
        }
      }
    }
  }

  for(int i = 0; i < N_act; ++i) {
    for(int j = 0; j <= i; ++j) {
      int ij = i*(i+1)/2+j;

      if(fabs(h0[ij]) < 1.0e-16) continue;

      fout << setw(24) << scientific << h0[ij]
           << setw( 4) << i+1
           << setw( 4) << j+1
           << setw( 4) << 0
           << setw( 4) << 0 << endl;
    }
  }

  fout << setw(24) << scientific << E_core
       << setw( 4) << 0
       << setw( 4) << 0
       << setw( 4) << 0
       << setw( 4) << 0 << endl;

  fout.close();
}
