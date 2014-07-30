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
#include "onepdm_container.h"
#include "npdm_permutations.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

Onepdm_container::Onepdm_container( int sites )
{
  store_full_spin_array_ = true;
  store_full_spatial_array_ = true;

  if(mpigetrank() == 0) {
    if(store_full_spin_array_) {
      if(dmrginp.spinAdapted())
        onepdm.resize(2*sites);
      else
        onepdm.resize(sites);

      onepdm.fill(0.0);
    }
    if(store_full_spatial_array_) {
      if(dmrginp.spinAdapted())
        spatial_onepdm.resize(sites);
      else
        spatial_onepdm.resize(sites/2);

      spatial_onepdm.fill(0.0);
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
  
void Onepdm_container::save_npdms(const int& i, const int& j)
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
  pout << "1PDM save full array time " << timer.elapsedwalltime() << " " << timer.elapsedcputime() << endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Onepdm_container::save_npdm_text(const int &i, const int &j)
{
  if(mpigetrank() == 0) {
    char file[5000];
    sprintf(file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/onepdm.", i, j,".txt");
    ofstream ofs(file);

    size_t ext = onepdm.extent();

    ofs << ext << endl;
    double trace = 0.0;
    for(int i = 0; i < ext; ++i)
      for(int j = 0; j < ext; ++j) {
        double value = onepdm(i,j);
        if(i == j) trace += value;

        // printing non-redundant elements only
        if(i > j && abs(value) > NUMERICAL_ZERO) {
          ofs << boost::format("%d %d %20.14e\n") % i % j % value;
        }
      }

    ofs.close();
    std::cout << "Spin-orbital 1PDM trace = " << trace << "\n";
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Onepdm_container::save_spatial_npdm_text(const int &i, const int &j)
{
  if(mpigetrank() == 0) {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_onepdm.", i, j,".txt");
    ofstream ofs(file);

    size_t ext = spatial_onepdm.extent();

    ofs << ext << endl;
    double trace = 0.0;

//  for(auto it = spatial_onepdm.begin(); it != spatial_onepdm.end(); ++it) {
//    int i = it.index()[0];
//    int j = it.index()[1];
//    if(abs(*it) > NUMERICAL_ZERO) {
//      ofs << boost::format("%d %d %20.14e\n") % i % j % (*it);
//    }
//  }

    // printing full spatial array (no use of permutation symmetry)
    for(int i = 0; i < ext; ++i)
      for(int j = 0; j < ext; ++j) {
        double value = spatial_onepdm(i,j);
        if(i == j) trace += value;

        if(abs(value) > NUMERICAL_ZERO)
          ofs << boost::format("%d %d %20.14e\n") % i % j % value;
      }


    ofs.close();
    std::cout << "Spatial      1PDM trace = " << trace << "\n";
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Onepdm_container::save_npdm_binary(const int &i, const int &j)
{
  if(mpigetrank() == 0) {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/onepdm.", i, j,".bin");
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << onepdm;
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Onepdm_container::save_spatial_npdm_binary(const int &i, const int &j)
{
  if(mpigetrank() == 0) {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_onepdm.", i, j,".bin");
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << spatial_onepdm;
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Onepdm_container::build_full_spatial_array()
{
  if(mpigetrank() == 0) {
    // Take into account orbital reordering
    const std::vector<int>& ro = dmrginp.reorder_vector();

    size_t ext = spatial_onepdm.extent();

    for(int i = 0; i < ext; ++i) {
      int i2 = 2*i;
      for(int j = 0; j < ext; ++j) {
        int j2 = 2*j;
        double& value = spatial_onepdm(ro.at(i), ro.at(j));
        if(abs(value) == 0.0) {
          value = onepdm(i2  ,j2  ) + onepdm(i2+1,j2+1);
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Onepdm_container::update_array_component()
{
#ifndef SERIAL
  boost::mpi::communicator world;
  if(world.rank() == 0) {
#endif
    for (auto it = tmp_store_.begin(); it != tmp_store_.end(); ++it) {
      double value = it->second;
      if(abs(value) < NUMERICAL_ZERO) continue;

      assert(it->first.size() == 2);
      int i = it->first[0];
      int j = it->first[1];

      onepdm(i,j) = value;
    }
#ifndef SERIAL
    for(int iproc = 1; iproc < world.size(); ++iproc) {
      std::vector<std::pair<std::vector<int>, double> > tmp_recv;
      world.recv(iproc, iproc, tmp_recv);
      for(auto it = tmp_recv.begin(); it != tmp_recv.end(); ++it) {
        double value = it->second;
        if(abs(value) < NUMERICAL_ZERO) continue;

        assert(it->first.size() == 2);
        int i = it->first[0];
        int j = it->first[1];

        onepdm(i,j) = value;
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

void Onepdm_container::store_npdm_elements(const std::vector<std::pair<std::vector<int>, double> > & in)
{
  // dims of spin-adapted -> non-spin-adapted transformation
  if(dmrginp.spinAdapted())
    assert(in.size() == 2);
  else
    assert(in.size() == 1);

  if(store_full_spin_array_ || store_full_spatial_array_) tmp_store_.insert(tmp_store_.end(), in.begin(), in.end());
}

//===========================================================================================================================================================

}
}

