#include "block_calldmrg.h"
#include "molpro_fcidump.h"

/// Fortran wrapper
void block_calldmrg_ (
      const int* Restart,
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
            double* E_sweep)
{
  block_calldmrg(*Restart, *N_roots, *N_act, *N_elec, *M_s, *E_core, h0, tuvx, *M_state, *N_pdm, *T_sweep, *T_noise, E_sweep);
}

extern int calldmrg(char*, char*);

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <vector>
using std::vector;
#include <sweep_params.h>

#include <communicate.h> // enable Boost MPI wrappers if defined SERIAL

#include <boost/filesystem.hpp>
#include <boost/serialization/string.hpp>
#include <boost/algorithm/string.hpp>

/// Dump 2-el integrals to formatted file
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
            double* E_sweep)
{
  using std::endl;
  using std::setw;

#ifndef SERIAL
  boost::mpi::communicator world;
#endif

//std::string symlab(Sym,3);
//boost::algorithm::tolower(symlab);

  std::string prefix;
  if(mpigetrank() == 0) prefix = boost::filesystem::current_path().c_str();

#ifndef SERIAL
  boost::mpi::broadcast(world, prefix, 0);
#endif

  char input[64];
  {
    std::ostringstream oss; oss << "dmrg.conf." << mpigetrank();
    strcpy(input, oss.str().c_str());
  }
  char output[64];
  {
    std::ostringstream oss; oss << "dmrg.out." << mpigetrank();
    strcpy(output, oss.str().c_str());
  }

// if(mpigetrank() == 0) // branch by mpi rank == 0
  if(1) // create config file for every procs/nodes
  {
    /// Dump integrals as MOLPRO format
    molpro_fcidump(N_act,N_elec,M_s,E_core,h0,tuvx);

    /// Create config file
    std::ofstream fcon(input);

    fcon << "nelec " << setw(2) << N_elec << endl;
    fcon << "spin  " << setw(2) << M_s << endl;
    fcon << "irrep  1" << endl; // nosymm for the meanwhile

    int N_sweep = 0;
    int M_start = (Restart == 1) ? M_state : 250;
    double T_start = T_noise;

    fcon << "schedule" << endl;
    while(M_start < M_state) {
      fcon << setw(2) << N_sweep << setw(5) << M_start << " " << T_start << " " << T_start << endl;
      N_sweep += 2;
      M_start *= 2;
    }
    while(T_start > 1.0e-6) {
      fcon << setw(2) << N_sweep << setw(5) << M_state << " " << T_start << " " << T_start << endl;
      N_sweep += 4;
      T_start /= 10;
    }
    while(T_start > T_sweep) {
      fcon << setw(2) << N_sweep << setw(5) << M_state << " " << T_start << " " << 0.0 << endl;
      N_sweep += 2;
      T_start /= 10;
    }
    fcon << setw(2) << N_sweep << setw(5) << M_state << " " << T_sweep << " 0.0" << endl;
    fcon << "end" << endl;
    fcon << "maxiter 100" << endl;

    if(N_sweep > 0)
      fcon << "twodot_to_onedot " << N_sweep+4 << endl;
    else
      fcon << "onedot" << endl;

    fcon << "sweep_tol " << T_sweep << endl;
    if(Restart == 1) {
      // FIXME:
      // when using restart for N_roots = 1, energy oscillation occurs somehow...
      // fullrestart with SA-DMRG and onedot fails
      if(N_sweep > 0 || N_roots == 1)
        fcon << "fullrestart" << endl;
      else
        fcon << "restart" << endl;

      fcon << "reset_iter" << endl;
    }

    switch (N_pdm) {
      case 1:
        fcon << "onepdm" << endl;
        break;
      case 2:
        fcon << "twopdm" << endl;
        break;
      case 3:
        fcon << "threepdm" << endl;
        fcon << "disk_dump_pdm" << endl;
        break;
      case 4:
        fcon << "fourpdm" << endl;
        fcon << "disk_dump_pdm" << endl;
        break;
      default:
        exit(1); // Block only supports 1-4RDMs
    }

//  fcon << "store_spinpdm" << endl;
//  fcon << "prefix " << prefix << endl;
    fcon << "orbitals FCIDUMP" << endl;
//  fcon << "symmetry " << symlab << endl;
    fcon << "gaopt default" << endl;
    fcon << "hf_occ integral" << endl;

    if(N_roots > 1) {
      fcon << "nroots " << N_roots << endl;
      fcon << "weights ";
      for(int i = 0; i < N_roots; ++i) fcon << 1.0/N_roots << " ";
      fcon << endl;
    }

//  fcon << "outputlevel 2" << endl;

    fcon.close();
  }

  if(mpigetrank() == 0) {
    boost::filesystem::path path_to_reorder("./RestartReorder.dat");
    if(Restart == 0 && boost::filesystem::exists(path_to_reorder))
      boost::filesystem::remove(path_to_reorder);
  }

  std::streambuf *backup;
  backup = std::cout.rdbuf();
  std::ofstream fout;

  fout.open(output);
  std::cout.rdbuf(fout.rdbuf());

  calldmrg(input,0);
//calldmrg(input,output);

  SpinAdapted::SweepParams param;
  bool forward; int size;
  param.restorestate(forward, size);

//assert(param.get_lowest_energy().size() == N_roots);
  for(int i = 0; i < N_roots; ++i)
    E_sweep[i] = param.get_lowest_energy()[i];

  // Sorting NPDMs as Chemist's order
  // TODO: This should be re-implemented as seperate subroutines
  if(N_pdm == 3 && mpigetrank() == 0) // branch by mpi rank == 0
  {
    boost::filesystem::path path_to_3pdm("./spatial_threepdm.0.0.bin");
    assert(boost::filesystem::exists(path_to_3pdm));

    int N1 = N_act;
    int N2 = N1*N1;
    int N3 = N2*N1;
    int N4 = N3*N1;
    int N5 = N4*N1;

    // sort 3pdm as chemist's notation in column-major order
    FILE *ifp = fopen("spatial_threepdm.0.0.bin","rb");
    FILE *ofp = fopen("SORTED3PDM.0","wb");

    // sort <i,j,k,l,m,n> (in row-major) to G(i,n,j,m,k,l) (in col-major)
    // i.e. G(i,j,k,l,m,n) = <i,k,m,n,l,j>

    // O(N^4) memory is allocated as buffer (lower overheads?)
    double *vbuf = new double[N4];
    for(int n = 0; n < N_act; ++n) {
      size_t Lxxxnxx = N2*n;
      for(int m = 0; m < N_act; ++m) {
        size_t Lxxmnxx = Lxxxnxx + N3*m;
        double *p = vbuf;
        for(int l = 0; l < N_act; ++l) {
          size_t Lxxmnlx = Lxxmnxx + N1*l;
          for(int k = 0; k < N_act; ++k) {
            size_t Lxkmnlx = Lxxmnlx + N4*k;
            for(int j = 0; j < N_act; ++j) {
              size_t Lxkmnlj = Lxkmnlx + j;
              for(int i = 0; i < N_act; ++i, ++p) {
                size_t Likmnlx = Lxkmnlj + N5*i;
                fseek(ifp,sizeof(double)*Likmnlx,SEEK_SET);
                fread(p,sizeof(double),1,ifp);
              }
            }
          }
        }
        fwrite(vbuf,sizeof(double),N4,ofp);
      }
    }
    // O(N^4+N^2) memory is allocated as buffer (least overheads?)
//  double *vbuf = new double[N4];
//  double *xbuf = new double[N2];
//  for(int n = 0; n < N_act; ++n) {
//    size_t Lxxxnxx = N2*n;
//    for(int m = 0; m < N_act; ++m) {
//      size_t Lxxmnxx = Lxxxnxx + N3*m;
//      double *p = vbuf;
//      for(int k = 0; k < N_act; ++k) {
//        size_t Lxkmnxx = Lxxmnxx + N4*k;
//        size_t Gxxkx = N2*k;
//        for(int i = 0; i < N_act; ++i) {
//          size_t Likmnxx = Lxkmnxx + N5*i;
//          size_t Gixkx = Gxxkx + i;
//          fseek(ifp,sizeof(double)*Likmnxx,SEEK_SET);
//          fread(xbuf,sizeof(double),N2,ifp);
//          double *p = xbuf;
//          for(int l = 0; l < N_act; ++l) {
//            size_t Gixkl = Gixkx + N3*l;
//            for(int j = 0; j < N_act; ++j, ++p) {
//              size_t Gijkl = Gixkx + N1*j;
//              vbuf[Gijkl] = *p;
//            }
//          }
//        }
//      }
//      fwrite(vbuf,sizeof(double),N4,ofp);
//    }
//  }
//  delete [] xbuf;
    delete [] vbuf;

    fclose(ifp);
    fclose(ofp);

    boost::filesystem::remove(path_to_3pdm);
  } // if(N_pdm == 3 && mpigetrank() == 0)

  std::cout.rdbuf(backup);
  fout.close();
}
