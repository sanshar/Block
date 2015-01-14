#ifndef SPIN_OVERLAP_TENSOR_HEADER
#define SPIN_OVERLAP_TENSOR_HEADER 
#include "btas/SPARSE/STArray.h"
#include "pario.h"
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
using namespace boost;
using namespace std;

namespace SpinAdapted {

template<class T, size_t N>
void SaveOverlapTensor (const std::vector<int>& sites, const btas::STArray<T, N>& m1, int state1, int state2)
{
  int rank = mpigetrank();
  if (rank == 0)
    {
      int first = min(sites[0], *sites.rbegin()), last = max(sites[0], *sites.rbegin());
      char file [5000];
      sprintf (file, "%s%s%d%s%d%s%d%s%d%s%d%s", dmrginp.save_prefix().c_str(), "/Overlap-", first, "-", last, ".", mpigetrank(),".state",state1,".",state2, ".tmp");
      p1out << "\t\t\t Saving Overlap Tensor :: " << file << endl;
      std::ofstream ofs(file, std::ios::binary);
      boost::archive::binary_oarchive save_mat(ofs);
      save_mat << m1;
      ofs.close();
    }
}

template<class T, size_t N>
void LoadOverlapTensor (const std::vector<int>& sites, btas::STArray<T, N>& m1, int state1, int state2)
{
  int rank = mpigetrank();
  if (rank == 0)
    {
      int first = min(sites[0], *sites.rbegin()), last = max(sites[0], *sites.rbegin());
      char file [5000];
      sprintf (file, "%s%s%d%s%d%s%d%s%d%s%d%s", dmrginp.save_prefix().c_str(), "/Overlap-", first, "-", last, ".", mpigetrank(),".state",state1,".",state2, ".tmp");
      p1out << "\t\t\t Loading Overlap Tensor :: " << file << endl;

      std::ifstream ifs(file, std::ios::binary);
      boost::archive::binary_iarchive load_mat(ifs);
      load_mat >> m1;
      ifs.close();
    }
}
};

#endif


