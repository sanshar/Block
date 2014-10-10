#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <cassert>

#include <boost/archive/binary_iarchive.hpp>

#include <boost/serialization/vector.hpp> // is missing in Block/include

#include <communicate.h> // enable Boost MPI wrappers if defined SERIAL
#include <multiarray.h> // from Block/include

#include "block_densi_rasscf.h"

void block_densi_rasscf_ (const int* iRoot, const int* N_act, const int* N_elec, double* D, double* DS, double* P, double* PA)
{
  block_densi_rasscf (*iRoot, *N_act, *N_elec, D, DS, P, PA);
}

void print_matrix(const int& N, const double* X)
{
  std::cout.precision(6);
  for(int i = 0; i < N; ++i) {
    std::cout << std::setw(2) << i << ":";
    for(int j = 0; j < N; ++j, ++X) {
      std::cout << std::fixed << std::setw(10) << *X;
    }
    std::cout << std::endl;
  }
}

void print_trimatrix(const int& N, const double* X)
{
  std::cout.precision(6);
  for(int i = 0; i < N; ++i, ++X) {
    std::cout << std::setw(2) << i << ":";
    for(int j = 0; j < i; ++j, ++X) {
      std::cout << std::fixed << std::setw(10) << *X;
    }
    std::cout << std::fixed << std::setw(10) << *X << std::endl;
  }
}

void block_load_twopdm_driver (const int& iRoot, const int& N_act, const int& N_elec, double* DA, double* DB, double* PT)
{
  std::ostringstream filestr;
  filestr << "spatial_binary_twopdm." << iRoot-1 << "." << iRoot-1 << ".bin";

  FILE *fp = fopen(filestr.str().c_str(),"rb");
  int NDUM;
  int istat = fread(&NDUM,sizeof(int),1,fp);

  assert(NDUM == N_act);

  int NMO1 = N_act;
  int NMO2 = NMO1*NMO1;
  int NMO3 = NMO2*NMO1;
  int NMO4 = NMO2*NMO2;

  double *pdm = new double[NMO4];
  int result = fread(pdm,sizeof(double),NMO4,fp);

  // load twopdm from disk and store in chemists ordering
  // <0|ci cj' dk' dl|0> => <0|ci dl cj' dk'|0>
  for(int i = 0; i < NMO1; ++i) {

//  int ix = reorder.at(i);

    int i000 = i*NMO3;

    for(int j = 0; j < NMO1; ++j) {

//    int jx = reorder.at(j);

      int ij = i*NMO1+j;

      int ij00 = i000+j*NMO2;
      int i00j = i000+j;

      DA[ij] = 0.0;

      for(int k = 0; k < NMO1; ++k) {

//      int kx = reorder.at(k);

        int ijk0 = ij00+k*NMO1;
        int ik0j = i00j+k*NMO2;

        for(int l = 0; l < NMO1; ++l) {

//        int lx = reorder.at(l);

          int ijkl = ijk0+l;
          int iklj = ik0j+l*NMO1;

          PT[ijkl] = pdm[iklj];
        }

        int ikkj = ik0j+k*NMO1;

        DA[ij] += pdm[ikkj];
      }

      DA[ij] /= N_elec-1;
      DB[ij]  = DA[ij];
    }
  }

  delete [] pdm;
  fclose(fp);
}

void block_densi_rasscf (const int& iRoot, const int& N_act, const int& N_elec, double* D, double* DS, double* P, double* PA)
{
  int NMO1 = N_act;
  int NMO2 = NMO1*NMO1;
  int NMO4 = NMO2*NMO2;

#ifndef SERIAL
  boost::mpi::communicator world;
#endif

  std::vector<double> DA;
  std::vector<double> DB;
  std::vector<double> PT;

  if(mpigetrank() == 0) {

    DA.resize(NMO2, 0.0);
    DB.resize(NMO2, 0.0);
    PT.resize(NMO4, 0.0);

    // load twopdm from disk and store in chemists ordering
    // <0|ci cj' dk' dl|0> => <0|ci dl cj' dk'|0>
    block_load_twopdm_driver(iRoot, N_act, N_elec, DA.data(), DB.data(), PT.data());
  }
#ifndef SERIAL
  boost::mpi::broadcast(world, DA, 0);
  boost::mpi::broadcast(world, DB, 0);
  boost::mpi::broadcast(world, PT, 0);
#endif

  // compute D (1-electron DM) and DS (spin DM) in MO
  for(int i = 0; i < NMO1; ++i) {
    for(int j = 0; j < i; ++j) {
      int ij = i*NMO1+j;
      int ji = j*NMO1+i;
      int ij_pack = i*(i+1)/2+j;
      D [ij_pack] = (DA[ij]+DA[ji]+DB[ij]+DB[ji])*0.5;  // average ij and ji
      DS[ij_pack] = (DA[ij]+DA[ji]-DB[ij]-DB[ji])*0.25; // average ij and ji
    }
    int ii = i*NMO1+i;
    int ii_pack = i*(i+1)/2+i;
    D [ii_pack] = (DA[ii]+DB[ii]);
    DS[ii_pack] = (DA[ii]-DB[ii])*0.5;
  }

//if(mpigetrank() == 0) {
//  std::cout << "1-particle reduced density matrix:" << std::endl;
//  print_trimatrix(NMO1,D);
//  std::cout << std::endl;

//  std::cout << "Spin density matrix:" << std::endl;
//  print_trimatrix(NMO1,DS);
//  std::cout << std::endl;
//}

  int NMO2_pack = NMO1*(NMO1+1)/2;
  int NMO4_pack = NMO2_pack*(NMO2_pack+1)/2;

  for(int i = 0; i < NMO1; ++i) {
    for(int j = 0; j < i; ++j) {

      int ij = i*NMO1+j;
      int ji = j*NMO1+i;
      int ij_pack = i*(i+1)/2+j;

      for(int k = 0; k < i; ++k) {
        for(int l = 0; l < k; ++l) {

          int kl = k*NMO1+l;
          int lk = l*NMO1+k;
          int kl_pack = k*(k+1)/2+l;

          int ijkl_pack = ij_pack*(ij_pack+1)/2+kl_pack;

          P [ijkl_pack] = (PT[ij*NMO2+kl]+PT[ij*NMO2+lk]+PT[ji*NMO2+kl]+PT[ji*NMO2+lk]
                          +PT[kl*NMO2+ij]+PT[kl*NMO2+ji]+PT[lk*NMO2+ij]+PT[lk*NMO2+ji])*0.25;
          PA[ijkl_pack] = (PT[ij*NMO2+kl]-PT[ij*NMO2+lk]-PT[ji*NMO2+kl]+PT[ji*NMO2+lk]
                          +PT[kl*NMO2+ij]-PT[kl*NMO2+ji]-PT[lk*NMO2+ij]+PT[lk*NMO2+ji])*0.25;
        }
        int kk = k*NMO1+k;
        int kk_pack = k*(k+1)/2+k;

        int ijkk_pack = ij_pack*(ij_pack+1)/2+kk_pack;

        P [ijkk_pack] = (PT[ij*NMO2+kk]+PT[ji*NMO2+kk]+PT[kk*NMO2+ij]+PT[kk*NMO2+ji])*0.25; //?
        PA[ijkk_pack] = 0.0;
      }

      for(int k = 0; k < j; ++k) {

        int ik = i*NMO1+k;
        int ki = k*NMO1+i;
        int ik_pack = i*(i+1)/2+k;

        int ijik_pack = ij_pack*(ij_pack+1)/2+ik_pack;

        P [ijik_pack] = (PT[ij*NMO2+ik]+PT[ij*NMO2+ki]+PT[ji*NMO2+ik]+PT[ji*NMO2+ki]
                        +PT[ik*NMO2+ij]+PT[ik*NMO2+ji]+PT[ki*NMO2+ij]+PT[ki*NMO2+ji])*0.25;
        PA[ijik_pack] = (PT[ij*NMO2+ik]-PT[ij*NMO2+ki]-PT[ji*NMO2+ik]+PT[ji*NMO2+ki]
                        +PT[ik*NMO2+ij]-PT[ik*NMO2+ji]-PT[ki*NMO2+ij]+PT[ki*NMO2+ji])*0.25;
      }

      int ijij_pack = ij_pack*(ij_pack+1)/2+ij_pack;

      P [ijij_pack] = (PT[ij*NMO2+ij]+PT[ij*NMO2+ji]+PT[ji*NMO2+ij]+PT[ji*NMO2+ji])*0.5;
      PA[ijij_pack] = (PT[ij*NMO2+ij]-PT[ij*NMO2+ji]-PT[ji*NMO2+ij]+PT[ji*NMO2+ji])*0.5;
    }
    int ii = i*NMO1+i;
    int ii_pack = i*(i+1)/2+i;

    for(int j = 0; j < i; ++j) {
      for(int k = 0; k < j; ++k) {

        int jk = j*NMO1+k;
        int kj = k*NMO1+j;
        int jk_pack = j*(j+1)/2+k;

        int iijk_pack = ii_pack*(ii_pack+1)/2+jk_pack;

        P [iijk_pack] = (PT[ii*NMO2+jk]+PT[ii*NMO2+kj]+PT[jk*NMO2+ii]+PT[kj*NMO2+ii])*0.5; //ok?
        PA[iijk_pack] = 0.0;
      }
      int jj = j*NMO1+j;
      int jj_pack = j*(j+1)/2+j;

      int iijj_pack = ii_pack*(ii_pack+1)/2+jj_pack;

      P [iijj_pack] = (PT[ii*NMO2+jj]+PT[jj*NMO2+ii])*0.5; //ok?
      PA[iijj_pack] = 0.0;
    }

    for(int j = 0; j < i; ++j) {

      int ij = i*NMO1+j;
      int ji = j*NMO1+i;
      int ij_pack = i*(i+1)/2+j;

      int iiij_pack = ii_pack*(ii_pack+1)/2+ij_pack;

      P [iiij_pack] = (PT[ii*NMO2+ij]+PT[ii*NMO2+ji]+PT[ij*NMO2+ii]+PT[ji*NMO2+ii])*0.5; //?
      PA[iiij_pack] = 0.0;
    }

    int iiii_pack = ii_pack*(ii_pack+1)/2+ii_pack;

    P [iiii_pack] = PT[ii*NMO2+ii]; //ok?
    PA[iiii_pack] = 0.0;
  }

//if(mpigetrank() == 0) {
//  std::cout << "Symmetric 2-particle density matrix:" << std::endl;
//  print_trimatrix(NMO2_pack,P);
//  std::cout << std::endl;

//  std::cout << "Anti-symmetric 2-particle density matrix:" << std::endl;
//  print_trimatrix(NMO2_pack,PA);
//  std::cout << std::endl;
//}
}
