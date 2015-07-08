#include <sstream>

#include <boost/filesystem.hpp>

#include <communicate.h> // enable Boost MPI wrappers if defined SERIAL

#include "sortNpdm.h"

void sort1pdm (FORTINT N_act, FORTINT iRoot, FORTINT jRoot)
{
  // Sorting 3PDMs as Chemist's order
  if(mpigetrank() == 0) {

    std::ostringstream ifname;
    ifname << "./node0/spatial_binary_onepdm." << iRoot << "." << jRoot << ".bin";

    std::ostringstream ofname;
    ofname << "./SORTED1PDM." << iRoot << "." << jRoot << ".0";

    boost::filesystem::path path_to_1pdm(ifname.str().c_str());
    assert(boost::filesystem::exists(path_to_1pdm));

    int N1 = N_act;
    int N2 = N1*N1;

    // sort 1pdm as chemist's notation in column-major order
    FILE *ifp = fopen(ifname.str().c_str(),"rb");
    FILE *ofp = fopen(ofname.str().c_str(),"wb");

    // sort <i,j> (in row-major) to G(i,j) (in col-major)
    // i.e. G(i,j) = <i,j>

    int Ndum;
    fread(&Ndum,sizeof(int),1,ifp);
    assert(Ndum == N_act);

    double *xbuf = new double[N2];
    fread(xbuf,sizeof(double),N2,ifp);

    double *vbuf = new double[N2];
    double *p = vbuf;
    for(int j = 0; j < N_act; ++j) {
      for(int i = 0; i < N_act; ++i, ++p) {
        *p = xbuf[N1*i+j];
      }
    }
    fwrite(vbuf,sizeof(double),N2,ofp);

    delete [] vbuf;
    delete [] xbuf;

    fclose(ifp);
    fclose(ofp);

    boost::filesystem::remove(path_to_1pdm);
  }
}

void sort2pdm (FORTINT N_act, FORTINT iRoot, FORTINT jRoot)
{
  // Sorting 3PDMs as Chemist's order
  if(mpigetrank() == 0) {

    std::ostringstream ifname;
    ifname << "./node0/spatial_binary_twopdm." << iRoot << "." << jRoot << ".bin";

    std::ostringstream ofname;
    ofname << "./SORTED2PDM." << iRoot << "." << jRoot << ".0";

    boost::filesystem::path path_to_2pdm(ifname.str().c_str());
    assert(boost::filesystem::exists(path_to_2pdm));

    int N1 = N_act;
    int N2 = N1*N1;
    int N3 = N2*N1;
    int N4 = N3*N1;

    // sort 2pdm as chemist's notation in column-major order
    FILE *ifp = fopen(ifname.str().c_str(),"rb");
    FILE *ofp = fopen(ofname.str().c_str(),"wb");

    // sort <i,j,k,l> (in row-major) to G(i,l,j,k) (in col-major)
    // i.e. G(i,j,k,l) = <i,k,l,j>

    int Ndum;
    fread(&Ndum,sizeof(int),1,ifp);
    assert(Ndum == N_act);

    double *xbuf = new double[N4];
    fread(xbuf,sizeof(double),N4,ifp);

    double *vbuf = new double[N2];
    for(int l = 0; l < N_act; ++l) {
      size_t Lxxlx = N1*l;
      for(int k = 0; k < N_act; ++k) {
        size_t Lxklx = Lxxlx + N2*k;
        double *p = vbuf;
        for(int j = 0; j < N_act; ++j) {
          size_t Lxklj = Lxklx + j;
          for(int i = 0; i < N_act; ++i, ++p) {
            size_t Liklj = Lxklj + N3*i;
            // Re-scale by 2, because twopdm module computes (1/2) d_iklj
            *p = 2.0*xbuf[Liklj];
          }
        }
        fwrite(vbuf,sizeof(double),N2,ofp);
      }
    }

    delete [] vbuf;
    delete [] xbuf;

    fclose(ifp);
    fclose(ofp);

    boost::filesystem::remove(path_to_2pdm);
  }
}

void sort3pdm (FORTINT N_act, FORTINT iRoot, FORTINT jRoot)
{
  // Sorting 3PDMs as Chemist's order
  if(mpigetrank() == 0) {

    std::ostringstream ifname;
    ifname << "./node0/spatial_threepdm." << iRoot << "." << jRoot << ".bin";

    std::ostringstream ofname;
    ofname << "./SORTED3PDM." << iRoot << "." << jRoot << ".0";

    boost::filesystem::path path_to_3pdm(ifname.str().c_str());
    assert(boost::filesystem::exists(path_to_3pdm));

    int N1 = N_act;
    int N2 = N1*N1;
    int N3 = N2*N1;
    int N4 = N3*N1;
    int N5 = N4*N1;

    // sort 3pdm as chemist's notation in column-major order
    FILE *ifp = fopen(ifname.str().c_str(),"rb");
    FILE *ofp = fopen(ofname.str().c_str(),"wb");

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

    delete [] vbuf;

    fclose(ifp);
    fclose(ofp);

    boost::filesystem::remove(path_to_3pdm);
  }
}

