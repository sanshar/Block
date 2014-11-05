#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>

#include <boost/filesystem.hpp>

#include <communicate.h> // enable Boost MPI wrappers if defined SERIAL

#include "loadNpdm.h"

void block_load1pdm_   (const FORTINT* N, double* V, const FORTINT* iRoot, const FORTINT* jRoot)
{ load1pdm(*N,V,*iRoot-1,*jRoot-1); }

void block_load2pdm_   (const FORTINT* N, double* V, const FORTINT* iRoot, const FORTINT* jRoot)
{ load2pdm(*N,V,*iRoot-1,*jRoot-1); }

void block_load2pdm2f_ (const FORTINT* N, double* V, const FORTINT* iRoot, const FORTINT* jRoot, const FORTINT* iP, const FORTINT* iQ)
{ load2pdm2f(*N,V,*iRoot-1,*jRoot-1,*iP-1,*iQ-1); }

void block_load3pdm_   (const FORTINT* N, double* V, const FORTINT* iRoot, const FORTINT* jRoot)
{ load3pdm(*N,V,*iRoot-1,*jRoot-1); }

void block_load3pdm2f_ (const FORTINT* N, double* V, const FORTINT* iRoot, const FORTINT* jRoot, const FORTINT* iP, const FORTINT* iQ)
{ load3pdm2f(*N,V,*iRoot-1,*jRoot-1,*iP-1,*iQ-1); }

void block_load3pdm4f_ (const FORTINT* N, double* V, const FORTINT* iRoot, const FORTINT* jRoot, const FORTINT* iP, const FORTINT* iQ, const FORTINT* jP, const FORTINT* jQ)
{ load3pdm4f(*N,V,*iRoot-1,*jRoot-1,*iP-1,*iQ-1,*jP-1,*jQ-1); }

void load1pdm (FORTINT N, double* V, FORTINT iRoot, FORTINT jRoot)
{
#ifndef SERIAL
  boost::mpi::communicator world;
#endif
  size_t N2 = N*N;
  if(mpigetrank() == 0) {
    std::ostringstream filestr;
    filestr << "./SORTED1PDM." << iRoot << "." << jRoot << ".0";

//  boost::filesystem::path path_to_1pdm(filestr.str().c_str());
//  assert(boost::filesystem::exists(path_to_1pdm));

    FILE *fp = fopen(filestr.str().c_str(),"rb");
    fread(V,sizeof(double),N2,fp);
    fclose(fp);
  }
#ifndef SERIAL
  boost::mpi::broadcast(world,V,N2,0);
#endif
}

void load2pdm (FORTINT N, double* V, FORTINT iRoot, FORTINT jRoot)
{
#ifndef SERIAL
  boost::mpi::communicator world;
#endif
  size_t N2 = N*N;
  size_t N4 = N2*N2;
  if(mpigetrank() == 0) {
    std::ostringstream filestr;
    filestr << "./SORTED2PDM." << iRoot << "." << jRoot << ".0";

//  boost::filesystem::path path_to_2pdm(filestr.str().c_str());
//  assert(boost::filesystem::exists(path_to_2pdm));

    FILE *fp = fopen(filestr.str().c_str(),"rb");
    fread(V,sizeof(double),N4,fp);
    fclose(fp);
  }
#ifndef SERIAL
  boost::mpi::broadcast(world,V,N4,0);
#endif
}

void load2pdm2f (FORTINT N, double* V, FORTINT iRoot, FORTINT jRoot, FORTINT iP, FORTINT iQ)
{
#ifndef SERIAL
  boost::mpi::communicator world;
#endif
  size_t N2 = N*N;
  if(mpigetrank() == 0) {
    std::ostringstream filestr;
    filestr << "./SORTED2PDM." << iRoot << "." << jRoot << ".0";

//  boost::filesystem::path path_to_2pdm(filestr.str().c_str());
//  assert(boost::filesystem::exists(path_to_2pdm));

    FILE *fp = fopen(filestr.str().c_str(),"rb");
    int iOffSet = iP+iQ*N;
    fseek(fp,sizeof(double)*iOffSet*N2,SEEK_SET);
    fread(V,sizeof(double),N2,fp);
    fclose(fp);
  }
#ifndef SERIAL
  boost::mpi::broadcast(world,V,N2,0);
#endif
}

void load3pdm (FORTINT N, double* V, FORTINT iRoot, FORTINT jRoot)
{
#ifndef SERIAL
  boost::mpi::communicator world;
#endif
  size_t N2 = N*N;
  size_t N4 = N2*N2;
  size_t N6 = N4*N2;
  if(mpigetrank() == 0) {
    std::ostringstream filestr;
    filestr << "./SORTED3PDM." << iRoot << "." << jRoot << ".0";

//  boost::filesystem::path path_to_3pdm(filestr.str().c_str());
//  assert(boost::filesystem::exists(path_to_3pdm));

    FILE *fp = fopen(filestr.str().c_str(),"rb");
    fread(V,sizeof(double),N6,fp);
    fclose(fp);
  }
#ifndef SERIAL
  boost::mpi::broadcast(world,V,N6,0);
#endif
}

void load3pdm2f (FORTINT N, double* V, FORTINT iRoot, FORTINT jRoot, FORTINT iP, FORTINT iQ)
{
#ifndef SERIAL
  boost::mpi::communicator world;
#endif
  size_t N2 = N*N;
  size_t N4 = N2*N2;
  if(mpigetrank() == 0) {
    std::ostringstream filestr;
    filestr << "./SORTED3PDM." << iRoot << "." << jRoot << ".0";

//  boost::filesystem::path path_to_3pdm(filestr.str().c_str());
//  assert(boost::filesystem::exists(path_to_3pdm));

    FILE *fp = fopen(filestr.str().c_str(),"rb");
    int iOffSet = iP+iQ*N;
    fseek(fp,sizeof(double)*iOffSet*N4,SEEK_SET);
    fread(V,sizeof(double),N4,fp);
    fclose(fp);
  }
#ifndef SERIAL
  boost::mpi::broadcast(world,V,N4,0);
#endif
}

void load3pdm4f (FORTINT N, double* V, FORTINT iRoot, FORTINT jRoot, FORTINT iP, FORTINT iQ, FORTINT jP, FORTINT jQ)
{
#ifndef SERIAL
  boost::mpi::communicator world;
#endif
  size_t N2 = N*N;
  size_t N3 = N*N2;
  if(mpigetrank() == 0) {
    std::ostringstream filestr;
    filestr << "./SORTED3PDM." << iRoot << "." << jRoot << ".0";

//  boost::filesystem::path path_to_3pdm(filestr.str().c_str());
//  assert(boost::filesystem::exists(path_to_3pdm));

    FILE *fp = fopen(filestr.str().c_str(),"rb");
    int iOffSet = iP+iQ*N+jP*N2+jQ*N3;
    fseek(fp,sizeof(double)*iOffSet*N2,SEEK_SET);
    fread(V,sizeof(double),N2,fp);
    fclose(fp);
  }
#ifndef SERIAL
  boost::mpi::broadcast(world,V,N2,0);
#endif
}

