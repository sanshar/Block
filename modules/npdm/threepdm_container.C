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
#include "threepdm_container.h"
#include "npdm_permutations.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

Threepdm_container::Threepdm_container(int sites)
{
  store_full_spin_array_ = true;
  store_full_spatial_array_ = true;

  if(mpigetrank() == 0) {
    if(store_full_spin_array_) {
      if(dmrginp.spinAdapted())
        threepdm.resize(2*sites);
      else
        threepdm.resize(sites);

      threepdm.fill(0.0);
    }
    if(store_full_spatial_array_) {
      if(dmrginp.spinAdapted())
        spatial_threepdm.resize(sites);
      else
        spatial_threepdm.resize(sites/2);

      spatial_threepdm.fill(0.0);
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::save_npdms(const int& i, const int& j)
{
#ifndef SERIAL
  boost::mpi::communicator world;
  world.barrier();
#endif
  Timer timer;

  save_npdm_binary(i,j);
  save_npdm_text(i,j);

  if(store_full_spatial_array_) {
    build_full_spatial_array();
    save_spatial_npdm_binary(i,j);
    save_spatial_npdm_text(i,j);
  }

#ifndef SERIAL
  world.barrier();
#endif
  pout << "3PDM save full array time " << timer.elapsedwalltime() << " " << timer.elapsedcputime() << endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::save_npdm_text(const int &i, const int &j)
{
  if(mpigetrank() == 0) {
    char file[5000];
    sprintf(file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/threepdm.", i, j,".txt");
    ofstream ofs(file);

    size_t ext = threepdm.extent();

    ofs << ext << endl;
    double trace = 0.0;
    for(int i = 0; i < ext; ++i)
      for(int j = 0; j < ext; ++j)
        for(int k = 0; k < ext; ++k)
          for(int l = 0; l < ext; ++l)
            for(int m = 0; m < ext; ++m)
              for(int n = 0; n < ext; ++n) {
                double value = threepdm(i,j,k,l,m,n);
                if(i == n && j == m && k == l) trace += value;

                // printing non-redundant elements only
                if(i > j && j > k && l > m && m > n) {
                  int ijk = i*(i-1)*(i-2)/6+j*(j-1)/2+k;
                  int lmn = l*(l-1)*(l-2)/6+m*(m-1)/2+n;
                  if(ijk > lmn && abs(value) > NUMERICAL_ZERO) {
                    ofs << boost::format("%d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % value;
                  }
                }
                // printing all elements for DEBUG
//              if(abs(value) > NUMERICAL_ZERO)
//                ofs << boost::format("%d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % value;
              }

    ofs.close();
    std::cout << "Spin-orbital 3PDM trace = " << trace << "\n";
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::save_spatial_npdm_text(const int &i, const int &j)
{
  if(mpigetrank() == 0) {
    char file[5000];
    sprintf(file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.", i, j,".txt");
    ofstream ofs(file);

    size_t ext = spatial_threepdm.extent();

    ofs << ext << endl;
    double trace = 0.0;

//  TODO: NN prefers to following, symmetric_spatial_iterator is to be implemented.

//  // printing spatial array (only non-redundant elements)
//  for(auto it = spatial_threepdm.begin(); it != spatial_threepdm.end(); ++it) {
//    int i = it.index()[0];
//    int j = it.index()[1];
//    int k = it.index()[2];
//    int l = it.index()[3];
//    int m = it.index()[4];
//    int n = it.index()[5];
//    if(abs(*it) > NUMERICAL_ZERO ) {
//      ofs << boost::format("%d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % (*it);
//    }
//  }

    // printing full all elements
    for(int i = 0; i < ext; ++i)
      for(int j = 0; j < ext; ++j)
        for(int k = 0; k < ext; ++k)
          for(int l = 0; l < ext; ++l)
            for(int m = 0; m < ext; ++m)
              for(int n = 0; n < ext; ++n) {
                double& value = spatial_threepdm(i,j,k,l,m,n);
                if(i == n  && j == m && k == l) trace += value;

                if(abs(value) > NUMERICAL_ZERO)
                  ofs << boost::format("%d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % value;
              }

    ofs.close();
    std::cout << "Spatial      3PDM trace = " << trace << "\n";
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::save_npdm_binary(const int &i, const int &j)
{
  if(mpigetrank() == 0) {
    char file[5000];
    sprintf(file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/threepdm.", i, j,".bin");
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << threepdm;
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::save_spatial_npdm_binary(const int &i, const int &j)
{
  if(mpigetrank() == 0) {
    char file[5000];
    sprintf(file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.", i, j,".bin");
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << spatial_threepdm;
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::load_npdm_binary(const int &i, const int &j) { abort(); }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::build_full_spatial_array()
{
  if(mpigetrank() == 0) {
    // Take into account orbital reordering
    const std::vector<int>& ro = dmrginp.reorder_vector();

    size_t ext = spatial_threepdm.extent();

    for(int i = 0; i < ext; ++i) {
      int i2 = 2*i;
      for(int j = 0; j < ext; ++j) {
        int j2 = 2*j;
        for(int k = 0; k < ext; ++k) {
          int k2 = 2*k;
          for(int l = 0; l < ext; ++l) {
            int l2 = 2*l;
            for(int m = 0; m < ext; ++m) {
              int m2 = 2*m;
              for(int n = 0; n < ext; ++n) {
                int n2 = 2*n;
                double& value = spatial_threepdm(ro.at(i), ro.at(j), ro.at(k), ro.at(l), ro.at(m), ro.at(n));
                if(abs(value) == 0.0) {
                  value = threepdm(i2  ,j2  ,k2  ,l2  ,m2  ,n2  )
                        + threepdm(i2+1,j2  ,k2  ,l2  ,m2  ,n2+1)
                        + threepdm(i2  ,j2+1,k2  ,l2  ,m2+1,n2  )
                        + threepdm(i2  ,j2  ,k2+1,l2+1,m2  ,n2  )
                        + threepdm(i2+1,j2+1,k2  ,l2  ,m2+1,n2+1)
                        + threepdm(i2+1,j2  ,k2+1,l2+1,m2  ,n2+1)
                        + threepdm(i2  ,j2+1,k2+1,l2+1,m2+1,n2  )
                        + threepdm(i2+1,j2+1,k2+1,l2+1,m2+1,n2+1);
                }
              }
            }
          }
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

/// this accumulates temporary storage to master process
void Threepdm_container::update_array_component()
{
#ifndef SERIAL
  boost::mpi::communicator world;
  if(world.rank() == 0) {
#endif
    for(auto it = tmp_store_.begin(); it != tmp_store_.end(); ++it) {
      double value = it->second;
      if(abs(value) < NUMERICAL_ZERO) continue;

      assert(it->first.size() == 6);
      int i = it->first[0];
      int j = it->first[1];
      int k = it->first[2];
      int l = it->first[3];
      int m = it->first[4];
      int n = it->first[5];

// DEBUG: whether or not duplication occurs
//      : NN checked there's no duplication in case (8e, 8o)
//    assert(abs(threepdm(i,j,k,l,m,n)) == 0.0);

      threepdm(i,j,k,l,m,n) = value;
    }
#ifndef SERIAL
    for(int iproc = 1; iproc < world.size(); ++iproc) {
      std::vector<std::pair<std::vector<int>, double> > tmp_recv;
      world.recv(iproc, iproc, tmp_recv);
      for(auto it = tmp_recv.begin(); it != tmp_recv.end(); ++it) {
        double value = it->second;
        if(abs(value) < NUMERICAL_ZERO) continue;

        assert(it->first.size() == 6);
        int i = it->first[0];
        int j = it->first[1];
        int k = it->first[2];
        int l = it->first[3];
        int m = it->first[4];
        int n = it->first[5];

        threepdm(i,j,k,l,m,n) = value;
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

void Threepdm_container::store_npdm_elements(const std::vector<std::pair<std::vector<int>, double> >& in)
{
  // dims of spin-adapted -> non-spin-adapted transformation
  if(dmrginp.spinAdapted())
    assert(in.size() == 20);
  else
    assert(in.size() == 1);

  if(store_full_spin_array_ || store_full_spatial_array_) tmp_store_.insert(tmp_store_.end(), in.begin(), in.end());
}

//===========================================================================================================================================================

}
}


