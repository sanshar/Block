#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

#include <communicate.h> // enable Boost MPI wrappers if defined SERIAL

#include "loadNpdm.h"

void load1pdm_   (const int* N, double* V)
{ load1pdm(*N,V); }

void load2pdm_   (const int* N, double* V)
{ load2pdm(*N,V); }

void load2pdm2f_ (const int* N, double* V, const int* iP, const int* iQ)
{ load2pdm2f(*N,V,*iP-1,*iQ-1); }

void load3pdm_   (const int* N, double* V)
{ load3pdm(*N,V); }

void load3pdm2f_ (const int* N, double* V, const int* iP, const int* iQ)
{ load3pdm2f(*N,V,*iP-1,*iQ-1); }

void load3pdm4f_ (const int* N, double* V, const int* iP, const int* iQ, const int* jP, const int* jQ)
{ load3pdm4f(*N,V,*iP-1,*iQ-1,*jP-1,*jQ-1); }

void load1pdm (int N, double* V)
{
#ifndef SERIAL
  boost::mpi::communicator world;
#endif
  size_t N2 = N*N;
  if(mpigetrank() == 0) {
    FILE *fp = fopen("SORTED1PDM.0","rb");
    fread(V,sizeof(double),N2,fp);
    fclose(fp);
  }
#ifndef SERIAL
  boost::mpi::broadcast(world,V,N2,0);
#endif
}

void load2pdm (int N, double* V)
{
#ifndef SERIAL
  boost::mpi::communicator world;
#endif
  size_t N2 = N*N;
  size_t N4 = N2*N2;
  if(mpigetrank() == 0) {
    FILE *fp = fopen("SORTED2PDM.0","rb");
    fread(V,sizeof(double),N4,fp);
    fclose(fp);
  }
#ifndef SERIAL
  boost::mpi::broadcast(world,V,N4,0);
#endif
}

void load2pdm2f (int N, double* V, int iP, int iQ)
{
#ifndef SERIAL
  boost::mpi::communicator world;
#endif
  size_t N2 = N*N;
  if(mpigetrank() == 0) {
    FILE *fp = fopen("SORTED2PDM.0","rb");
    int iOffSet = iP+iQ*N;
    fseek(fp,sizeof(double)*iOffSet*N2,SEEK_SET);
    fread(V,sizeof(double),N2,fp);
    fclose(fp);
  }
#ifndef SERIAL
  boost::mpi::broadcast(world,V,N2,0);
#endif
}

void load3pdm (int N, double* V)
{
#ifndef SERIAL
  boost::mpi::communicator world;
#endif
  size_t N2 = N*N;
  size_t N4 = N2*N2;
  size_t N6 = N4*N2;
  if(mpigetrank() == 0) {
    FILE *fp = fopen("SORTED3PDM.0","rb");
    fread(V,sizeof(double),N6,fp);
    fclose(fp);
  }
#ifndef SERIAL
  boost::mpi::broadcast(world,V,N6,0);
#endif
}

void load3pdm2f (int N, double* V, int iP, int iQ)
{
#ifndef SERIAL
  boost::mpi::communicator world;
#endif
  size_t N2 = N*N;
  size_t N4 = N2*N2;
  if(mpigetrank() == 0) {
    FILE *fp = fopen("SORTED3PDM.0","rb");
    int iOffSet = iP+iQ*N;
    fseek(fp,sizeof(double)*iOffSet*N4,SEEK_SET);
    fread(V,sizeof(double),N4,fp);
    fclose(fp);
  }
#ifndef SERIAL
  boost::mpi::broadcast(world,V,N4,0);
#endif
}

void load3pdm4f (int N, double* V, int iP, int iQ, int jP, int jQ)
{
#ifndef SERIAL
  boost::mpi::communicator world;
#endif
  size_t N2 = N*N;
  size_t N3 = N*N2;
  if(mpigetrank() == 0) {
    FILE *fp = fopen("SORTED3PDM.0","rb");
    int iOffSet = iP+iQ*N+jP*N2+jQ*N3;
    fseek(fp,sizeof(double)*iOffSet*N2,SEEK_SET);
    fread(V,sizeof(double),N2,fp);
    fclose(fp);
  }
#ifndef SERIAL
  boost::mpi::broadcast(world,V,N2,0);
#endif
}

