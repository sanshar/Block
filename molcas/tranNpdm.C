#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <iostream>

//#include <mkl_cblas.h>

#include <boost/filesystem.hpp>

#include <communicate.h> // enable Boost MPI wrappers if defined SERIAL
#include <blas_calls.h>

#include "tranNpdm.h"

void block_tran1pdm_ (const FORTINT* N, double* X, const FORTINT* iRoot, const FORTINT* jRoot) { tran1pdm(*N, X, *iRoot-1, *jRoot-1); }
void block_tran2pdm_ (const FORTINT* N, double* X, const FORTINT* iRoot, const FORTINT* jRoot) { tran2pdm(*N, X, *iRoot-1, *jRoot-1); }
void block_tran3pdm_ (const FORTINT* N, double* X, const FORTINT* iRoot, const FORTINT* jRoot) { tran3pdm(*N, X, *iRoot-1, *jRoot-1); }

void tran1pdm (FORTINT N, double* X, FORTINT iRoot, FORTINT jRoot)
{
#ifndef SERIAL
  boost::mpi::communicator world;
#endif
  if(mpigetrank() == 0) {

    std::ostringstream filestr;
    filestr << "./SORTED1PDM." << iRoot << "." << jRoot << ".0";

    size_t N2 = N *N;

    FILE *ifp = fopen(filestr.str().c_str(),"rb");

    double *V1 = new double[N2];
    double *G1 = new double[N2];

    if(fread(G1, sizeof(double), N2, ifp) != N2) exit(1);

    DGEMM('T','N',N,N,N,1.0,X,N,G1,N,0.0,V1,N);
    DGEMM('N','N',N,N,N,1.0,V1,N,X,N,0.0,G1,N);

    fclose(ifp);

    FILE *ofp = fopen(filestr.str().c_str(), "wb");

    fwrite(G1, sizeof(double), N2, ofp);

    fclose(ofp);

    delete [] G1;
    delete [] V1;
  }
}

void tran2pdm (FORTINT N, double* X, FORTINT iRoot, FORTINT jRoot)
{
  SF_tran2pdm(N, X, iRoot, jRoot);
  SF_tran2pdm(N, X, iRoot, jRoot);
  SF_tran2pdm(N, X, iRoot, jRoot);
  SF_tran2pdm(N, X, iRoot, jRoot);
}

void tran3pdm (FORTINT N, double* X, FORTINT iRoot, FORTINT jRoot)
{
  SF_tran3pdm(N, X, iRoot, jRoot);
  SF_tran3pdm(N, X, iRoot, jRoot);
  SF_tran3pdm(N, X, iRoot, jRoot);
  SF_tran3pdm(N, X, iRoot, jRoot);
  SF_tran3pdm(N, X, iRoot, jRoot);
  SF_tran3pdm(N, X, iRoot, jRoot);
}

void SF_tran2pdm (FORTINT N, double* X, FORTINT iRoot, FORTINT jRoot)
{
#ifndef SERIAL
  boost::mpi::communicator world;
#endif 
  if(mpigetrank() == 0) {

    std::ostringstream filestr;
    filestr << "./SORTED2PDM." << iRoot << "." << jRoot << ".0";

    size_t N2 = N *N;
    size_t N3 = N2*N;
    size_t N4 = N3*N;

    size_t i,j,k,l,p;

    FILE *ifp = fopen(filestr.str().c_str(),"rb");

    double *V2 = new double[N4];
    double *G2 = new double[N];

    p = 0;

    for(i = 0; i < N; ++i)
      for(j = 0; j < N; ++j)
        for(k = 0; k < N; ++k, ++p) {
          if(fread(G2,sizeof(double),N,ifp) != N) exit(1);
          DGEMV('T',N,N,1.0,X,N,G2,1,0.0,V2+p,N3);
        }

    fclose(ifp);

    FILE *ofp = fopen(filestr.str().c_str(),"wb");

    fwrite(V2,sizeof(double),N4,ofp);

    fclose(ofp);

    delete [] G2;
    delete [] V2;
  }
}

void SF_tran3pdm (FORTINT N, double* X, FORTINT iRoot, FORTINT jRoot)
{
#ifndef SERIAL
  boost::mpi::communicator world;
#endif

  if(mpigetrank() == 0) {

    int istat;

    std::ostringstream ifname;
    ifname << "./SORTED3PDM." << iRoot << "." << jRoot << "." << mpigetrank();

    FILE *ifp = fopen(ifname.str().c_str(),"rb");

    std::ostringstream ofname;
    ofname << "./SORTED3PDM." << iRoot << "." << jRoot << "." << mpigetrank() << ".scr";

    FILE *ofp = fopen(ofname.str().c_str(),"wb");

    int pOff = N*N*N*N*N-1;

    double *T3 = new double[N];
    double *G3 = new double[N];

    size_t p = 0;

    for(int i = 0; i < N; ++i)
      for(int j = 0; j < N; ++j)
        for(int k = 0; k < N; ++k)
          for(int l = 0; l < N; ++l)
            for(int m = 0; m < N; ++m, ++p) {
              // load partial vec from disk
              istat = fread(G3,sizeof(double),N,ifp);
              // rotate vec by X
              DGEMV('T',N,N,1.0,X,N,G3,1,0.0,T3,1);
              // sort and store vec on disk
              fseek(ofp,sizeof(double)*p,SEEK_SET);
              for(int n = 0; n < N; ++n) {
                fwrite(&T3[n],sizeof(double),1,ofp);
                fseek(ofp,sizeof(double)*pOff,SEEK_CUR);
              }
            }

    delete [] G3;
    delete [] T3;

    fclose(ifp);
    fclose(ofp);

    // move -f ./SORTED3PDM.i.j.x.scr ./SORTED3PDM.i.j.x
    boost::filesystem::path oldPath(ofname.str().c_str());
    boost::filesystem::path newPath(ifname.str().c_str());
    boost::filesystem::rename(oldPath, newPath);
  } // mpigetrank() == 0
}

/*
void para_SF3_tran3pdm (FORTINT N, double* X, FORTINT iRoot, FORTINT jRoot)
{
  boost::mpi::communicator world;

  FILE *ifp, *ofp;
  int istat;

  std::ostringstream ifname;
  ifname << "./SORTED3PDM." << iRoot << "." << jRoot << "." << mpigetrank();

  ifp = fopen(ifname.str().c_str(),"rb");

  std::ostringstream ofname;
  ofname << "./SORTED3PDM." << iRoot << "." << jRoot << "." << mpigetrank() << ".scr";

  ofp = fopen(ofname.str().c_str(),"wb");

  int N3 = N*N*N;

  double *V3 = new double[N3];
  double *G3 = new double[N3];

  int nRank = world.size();
  int myRank = mpigetrank();

  int NG3 = pG3[myRank+1]-pG3[myRank];
  double *T3 = new double[NG3*nRank];

  for(int iG3 = 0; iG3 < MxG3; ++iG3) {

    if(iG3 < NG3) {
      istat = fread(V3,sizeof(double),N3,ifp);
      DGEMM('T','T',N,N2,N,1.0,X,N,V3,N2,0.0,G3);
      DGEMM('T','T',N,N2,N,1.0,X,N,G3,N2,0.0,V3);
      DGEMM('T','T',N,N2,N,1.0,X,N,V3,N2,0.0,G3);

      for(int iRank = 0; iRank < world.size(); ++iRank) {
        if(iRank != myRank) {
          // sending data
          int LG3 = pG3[iRank+1]-pG3[iRank];
          world.send(iRank,myRank,G3+pG3[iRank],LG3);
        }
      }
    }

    for(int iRank = 0; iRank < world.size(); ++iRank) {
      int LG3 = pG3[iRank+1]-pG3[iRank];
      if(iG3 < LG3) {
        if(iRank == myRank) {
          dcopy(LG3,G3+pG3[iRank],1,T3+iRank*LG3,1);
        }
        else {
          // receiving data
          world.recv(iRank,myRank,T3+iRank*NG3,NG3);
        }
        fseek(ofp,sizeof(double)*(pG3[iRank]+iG3),SEEK_SET);
        for(int jG3 = 0; jG3 < NG3; ++jG3) {
          fwrite(&T3[jG3+iRank*NG3],sizeof(double),1,ofp);
          fseek(ofp,sizeof(double)*N3,SEEK_CUR);
        }
      }
    }
  }

  delete [] V3;
  delete [] G3;
  delete [] T3;

  fclose(ifp);
  fclose(ofp);

  // move -f ./SORTED3PDM.i.j.x.scr ./SORTED3PDM.i.j.x
  boost::filesystem::path oldPath(ofname.str().c_str());
  boost::filesystem::path newPath(ifname.str().c_str());
  boost::filesystem::rename(oldPath, newPath);
}
*/
