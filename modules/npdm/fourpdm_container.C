/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "global.h"
#include <boost/format.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include "fourpdm_container.h"
#include "npdm_permutations.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

Fourpdm_container::Fourpdm_container(int sites)
{
  store_full_spin_array_ = true;
  store_full_spatial_array_ = true;

  if(mpigetrank() == 0) {
    if(store_full_spin_array_) {
      if(dmrginp.spinAdapted())
        fourpdm.resize(2*sites);
      else
        fourpdm.resize(sites);

      fourpdm.fill(0.0);
    }
    if(store_full_spatial_array_) {
//    if(dmrginp.spinAdapted())
//      spatial_fourpdm.resize(sites);
//    else
//      spatial_fourpdm.resize(sites/2);

//    spatial_fourpdm.fill(0.0);
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::save_npdms(const int& i, const int& j)
{
#ifndef SERIAL
  boost::mpi::communicator world;
  world.barrier();
#endif
  Timer timer;

  save_npdm_binary(i,j);
  save_npdm_text(i,j);

//save_spatial_npdm_binary(i,j);
//save_spatial_npdm_text(i,j);

#ifndef SERIAL
  world.barrier();
#endif
  pout << "4PDM save full array time " << timer.elapsedwalltime() << " " << timer.elapsedcputime() << endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::save_npdm_text(const int &i, const int &j)
{
  if(mpigetrank() == 0) {
    char file[5000];
    sprintf(file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/fourpdm.", i, j,".txt");
    ofstream ofs(file);

    size_t ext = fourpdm.extent();

    ofs << ext << endl;
    double trace = 0.0;
    for(int i = 0; i < ext; ++i)
      for(int j = 0; j < ext; ++j)
        for(int k = 0; k < ext; ++k)
          for(int l = 0; l < ext; ++l)
            for(int m = 0; m < ext; ++m)
              for(int n = 0; n < ext; ++n)
                for(int p = 0; p < ext; ++p)
                  for(int q = 0; q < ext; ++q) {
                    double value = fourpdm(i,j,k,l,m,n,p,q);
                    if(i == q && j == p && k == n && l == m) trace += value;

                    // printing non-redundant elements only
                    if(i > j && j > k && k > l && m > n && n > p && p > q) {
                      int ijkl = i*(i-1)*(i-2)*(i-3)/24+j*(j-1)*(j-2)/6+k*(k-1)/2+l;
                      int mnpq = m*(m-1)*(m-2)*(m-3)/24+n*(n-1)*(n-2)/6+p*(p-1)/2+q;
                      if(ijkl > mnpq && abs(value) > NUMERICAL_ZERO) {
                        ofs << boost::format("%d %d %d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % p % q % value;
                      }
                    }
                  }

    ofs.close();
    std::cout << "Spin-orbital 4PDM trace = " << trace << "\n";
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::save_spatial_npdm_text(const int &i, const int &j)
{
  if(mpigetrank() == 0) {
    char file[5000];
    sprintf(file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_fourpdm.", i, j,".txt");
    ofstream ofs(file);

//  ofs << spatial_fourpdm.extent() << endl;

//  for(auto it = spatial_onepdm.begin(); it != spatial_onepdm.end(); ++it) {
//    int i = it.index()[0];
//    int j = it.index()[1];
//    int k = it.index()[2];
//    int l = it.index()[3];
//    int m = it.index()[4];
//    int n = it.index()[5];
//    int p = it.index()[6];
//    int q = it.index()[7];

//    if(abs(*it) > NUMERICAL_ZERO ) {
//      ofs << boost::format("%d %d %d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % p % q % (*it);
//    }
//  }

    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::save_npdm_binary(const int &i, const int &j)
{
  if(mpigetrank() == 0) {
    char file[5000];
    sprintf(file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/fourpdm.", i, j,".bin");
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << fourpdm;
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::save_spatial_npdm_binary(const int &i, const int &j)
{
  if(mpigetrank() == 0) {
    char file[5000];
    sprintf(file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_fourpdm.", i, j,".bin");
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
//  save << spatial_fourpdm;
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::load_npdm_binary(const int &i, const int &j) { abort(); }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

/// this accumulates temporary storage to master process
void Fourpdm_container::update_array_component()
{
#ifndef SERIAL
  boost::mpi::communicator world;
  if(world.rank() == 0) {
#endif
    for(auto it = tmp_store_.begin(); it != tmp_store_.end(); ++it) {
      double value = it->second;
      if(abs(value) < NUMERICAL_ZERO) continue;

      assert(it->first.size() == 8);
      int i = it->first[0];
      int j = it->first[1];
      int k = it->first[2];
      int l = it->first[3];
      int m = it->first[4];
      int n = it->first[5];
      int p = it->first[6];
      int q = it->first[7];

      // DEBUG: whether or not duplication occurs
      // FIXME: Really need? -> No? checked with (8e, 8o)
      assert(abs(fourpdm(i,j,k,l,m,n,p,q)) == 0.0);

      fourpdm(i,j,k,l,m,n,p,q) = value;
    }
#ifndef SERIAL
    for(int iproc = 1; iproc < world.size(); ++iproc) {
      std::vector<std::pair<std::vector<int>, double> > tmp_recv;
      world.recv(iproc, iproc, tmp_recv);
      for(auto it = tmp_recv.begin(); it != tmp_recv.end(); ++it) {
        double value = it->second;
        if(abs(value) < NUMERICAL_ZERO) continue;

        assert(it->first.size() == 8);
        int i = it->first[0];
        int j = it->first[1];
        int k = it->first[2];
        int l = it->first[3];
        int m = it->first[4];
        int n = it->first[5];
        int p = it->first[6];
        int q = it->first[7];

        fourpdm(i,j,k,l,m,n,p,q) = value;
      }
    }
  }
  else {
    world.send(0, mpigetrank(), tmp_store_);
  }
#endif
  tmp_store_.clear();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::store_npdm_elements(const std::vector<std::pair<std::vector<int>, double> >& in)
{
  // dims of spin-adapted -> non-spin-adapted transformation
  if(dmrginp.spinAdapted())
    assert(in.size() == 70);
  else
    assert(in.size() == 1);

  if(store_full_spin_array_ || store_full_spatial_array_) tmp_store_.insert(tmp_store_.end(), in.begin(), in.end());
}

//===========================================================================================================================================================

}
}


