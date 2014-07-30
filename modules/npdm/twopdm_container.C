/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "global.h"
#include "timer.h"
#include <boost/format.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include "twopdm_container.h"
#include "npdm_permutations.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

Twopdm_container::Twopdm_container( int sites )
{
  store_full_spin_array_ = true;
  store_full_spatial_array_ = true;

  if(mpigetrank() == 0) {
    if(store_full_spin_array_) {
      if(dmrginp.spinAdapted())
        twopdm.resize(2*sites);
      else
        twopdm.resize(sites);

      twopdm.fill(0.0);
    }
    if(store_full_spatial_array_) {
      if(dmrginp.spinAdapted())
        spatial_twopdm.resize(sites);
      else
        spatial_twopdm.resize(sites/2);

      spatial_twopdm.fill(0.0);
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
  
void Twopdm_container::save_npdms(const int& i, const int& j)
{
#ifndef SERIAL
  boost::mpi::communicator world;
  world.barrier();
#endif
  Timer timer;

  save_npdm_binary(i, j);
  save_npdm_text(i, j);

  if(store_full_spatial_array_) {
    build_full_spatial_array();
    save_spatial_npdm_binary(i, j);
    save_spatial_npdm_text(i, j);
  }

#ifndef SERIAL
  world.barrier();
#endif
  pout << "2PDM save full array time " << timer.elapsedwalltime() << " " << timer.elapsedcputime() << endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_container::save_npdm_text(const int &i, const int &j)
{
  if(mpigetrank() == 0) {
    char file[5000];
    sprintf(file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/twopdm.", i, j,".txt");
    ofstream ofs(file);

    size_t ext = twopdm.extent();

    ofs << ext << endl;
    double trace = 0.0;
    for(int i = 0; i < ext; ++i)
      for(int j = 0; j < ext; ++j)
        for(int k = 0; k < ext; ++k)
          for(int l = 0; l < ext; ++l) {
            double value = twopdm(i,j,k,l);
            if(i == l && j == k) trace += value;

            // printing non-redundant elements only
            if(i > j && k > l) {
              int ij = i*(i-1)/2+j;
              int kl = k*(k-1)/2+l;
              if(ij > kl && abs(value) > NUMERICAL_ZERO) {
                ofs << boost::format("%d %d %d %d %20.14e\n") % i % j % k % l % value;
              }
            }
          }

    ofs.close();
    std::cout << "Spin-orbital 2PDM trace = " << trace << "\n";
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_container::save_spatial_npdm_text(const int &i, const int &j)
{
  if(mpigetrank() == 0) {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_twopdm.", i, j,".txt");
    ofstream ofs(file);

    size_t ext = spatial_twopdm.extent();

    ofs << ext << endl;
    double trace = 0.0;

    ofs << spatial_twopdm.extent() << endl;

//  for(auto it = spatial_twopdm.begin(); it != spatial_twopdm.end(); ++it) {
//    int i = it.index()[0];
//    int j = it.index()[1];
//    int k = it.index()[2];
//    int l = it.index()[3];
//    if(abs(*it) > NUMERICAL_ZERO) {
//      ofs << boost::format("%d %d %d %d %20.14e\n") % i % j % k % l % (*it);
//    }
//  }

    // printing full spatial array (no use of permutation symmetry)
    for(int i = 0; i < ext; ++i)
      for(int j = 0; j < ext; ++j)
        for(int k = 0; k < ext; ++k)
          for(int l = 0; l < ext; ++l) {
            double value = spatial_twopdm(i,j,k,l);
            if(i == l  && j == k) trace += value;

            if(abs(value) > NUMERICAL_ZERO)
              ofs << boost::format("%d %d %d %d %20.14e\n") % i % j % k % l % value;
          }

    ofs.close();
    std::cout << "Spatial      2PDM trace = " << trace << "\n";
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_container::save_npdm_binary(const int &i, const int &j)
{
  if(mpigetrank() == 0) {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/twopdm.", i, j,".bin");
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << twopdm;
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_container::save_spatial_npdm_binary(const int &i, const int &j)
{
  if(mpigetrank() == 0) {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_twopdm.", i, j,".bin");
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << spatial_twopdm;
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_container::build_full_spatial_array()
{
  if(mpigetrank() == 0) {
    // Take into account orbital reordering
    const std::vector<int>& ro = dmrginp.reorder_vector();

    size_t ext = spatial_twopdm.extent();

    for(int i = 0; i < ext; ++i) {
      int i2 = 2*i;
      for(int j = 0; j < ext; ++j) {
        int j2 = 2*j;
        for(int k = 0; k < ext; ++k) {
          int k2 = 2*k;
          for(int l = 0; l < ext; ++l) {
            int l2 = 2*l;
            double& value = spatial_twopdm(ro.at(i), ro.at(j), ro.at(k), ro.at(l));
            if(abs(value) == 0.0) {
              value = twopdm(i2  ,j2  ,k2  ,l2  )
                    + twopdm(i2+1,j2  ,k2  ,l2+1)
                    + twopdm(i2  ,j2+1,k2+1,l2  )
                    + twopdm(i2+1,j2+1,k2+1,l2+1);
            }
          }
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_container::update_array_component()
{
#ifndef SERIAL
  boost::mpi::communicator world;
  if(world.rank() == 0) {
#endif
    for (auto it = tmp_store_.begin(); it != tmp_store_.end(); ++it) {
      double value = it->second;
      if(abs(value) < NUMERICAL_ZERO) continue;

      assert(it->first.size() == 4);
      int i = it->first[0];
      int j = it->first[1];
      int k = it->first[2];
      int l = it->first[3];

      twopdm(i,j,k,l) = value;
    }
#ifndef SERIAL
    for(int iproc = 1; iproc < world.size(); ++iproc) {
      std::vector<std::pair<std::vector<int>, double> > tmp_recv;
      world.recv(iproc, iproc, tmp_recv);
      for(auto it = tmp_recv.begin(); it != tmp_recv.end(); ++it) {
        double value = it->second;
        if(abs(value) < NUMERICAL_ZERO) continue;

        assert(it->first.size() == 4);
        int i = it->first[0];
        int j = it->first[1];
        int k = it->first[2];
        int l = it->first[3];

        twopdm(i,j,k,l) = value;
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

void Twopdm_container::store_npdm_elements(const std::vector<std::pair<std::vector<int>, double> > & in)
{
  // dims of spin-adapted -> non-spin-adapted transformation
  if(dmrginp.spinAdapted())
    assert(in.size() == 6);
  else
    assert(in.size() == 1);

  if(store_full_spin_array_ || store_full_spatial_array_) tmp_store_.insert(tmp_store_.end(), in.begin(), in.end());
}

//===========================================================================================================================================================

}
}

