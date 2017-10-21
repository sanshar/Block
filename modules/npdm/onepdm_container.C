/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "global.h"
#include "timer.h"
#include <boost/format.hpp>
#include "onepdm_container.h"
#include "npdm_permutations.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

Onepdm_container::Onepdm_container( int sites )
{

  if(dmrginp.spinAdapted()){
    onepdm.resize(2*sites,2*sites);
    spatial_onepdm.resize(sites,sites);

  }
  else{
    onepdm.resize(sites,sites);
    spatial_onepdm.resize(sites/2,sites/2);
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
  accumulate_npdm();
  save_npdm_binary(i, j);
  save_npdm_text(i, j);
  if(dmrginp.spinAdapted())
    accumulate_spatial_npdm();
  else
    calculate_spatial_npdm();
  save_spatial_npdm_binary(i, j);
  save_spatial_npdm_text(i, j);
#ifndef SERIAL
  world.barrier();
#endif
  ecpu = timer.elapsedcputime();ewall=timer.elapsedwalltime();
  p3out << "1PDM save full array time " << ewall << " " << ecpu << endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Onepdm_container::save_npdm_text(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/onepdm.", i, j,".txt");
    ofstream ofs(file);
    ofs << onepdm.dim1() << endl;
    double trace = 0.0;
    for(int k=0;k<onepdm.dim1();++k)
      for(int l=0;l<onepdm.dim2();++l) {
        ofs << boost::format("%d %d %20.14e\n") % k % l % onepdm(k,l);
        if ( k==l ) trace += onepdm(k,l);
      }
    ofs.close();
    pout << "Spin-orbital 1PDM trace = " << trace << "\n";
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Onepdm_container::save_spatial_npdm_text(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_onepdm.", i, j,".txt");
    ofstream ofs(file);
    ofs << spatial_onepdm.dim1() << endl;
    double trace = 0.0;
    for(int k=0;k<spatial_onepdm.dim1();++k)
      for(int l=0;l<spatial_onepdm.dim2();++l) {
        ofs << boost::format("%d %d %20.14e\n") % k % l % spatial_onepdm(k,l);
        if ( k==l ) trace += spatial_onepdm(k,l);
      }
    ofs.close();
    pout << "Spatial      1PDM trace = " << trace << "\n";
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Onepdm_container::save_npdm_binary(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
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
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_onepdm.", i, j,".bin");
#ifndef MOLCAS
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << spatial_onepdm;
    ofs.close();
#else
    FILE* f = fopen(file,"wb");
    fwrite(spatial_onepdm.data(),sizeof(double),spatial_onepdm.size(),f);
    fclose(f);
#endif
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Onepdm_container::accumulate_npdm()
{
#ifndef SERIAL
  mpi::communicator world;
#ifdef NDEBUG
  if( mpigetrank() == 0)
    MPI_Reduce(MPI_IN_PLACE, onepdm.data(),  onepdm.size(), MPI_DOUBLE, MPI_SUM, 0, world);
  else
    MPI_Reduce(onepdm.data(), onepdm.data(),  onepdm.size(), MPI_DOUBLE, MPI_SUM, 0, world);
#else
  array_2d<double> tmp_recv;
  if( mpigetrank() == 0) {
    for(int i=1;i<world.size();++i) {
      world.recv(i, i, tmp_recv);
      for(int k=0;k<onepdm.dim1();++k)
        for(int l=0;l<onepdm.dim2();++l)
          if( abs(tmp_recv(k,l)) > NUMERICAL_ZERO )
          {
          assert ( abs(onepdm(k,l)) < NUMERICAL_ZERO );
          onepdm(k,l) = tmp_recv(k,l);
          }
    }
  }
  else {
    world.send(0, mpigetrank(), onepdm);
  }
#endif
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Onepdm_container::accumulate_spatial_npdm()
{
#ifndef SERIAL
  mpi::communicator world;
#ifdef NDEBUG
  if( mpigetrank() == 0)
    MPI_Reduce(MPI_IN_PLACE, spatial_onepdm.data(),  spatial_onepdm.size(), MPI_DOUBLE, MPI_SUM, 0, world);
  else
    MPI_Reduce(spatial_onepdm.data(), spatial_onepdm.data(),  spatial_onepdm.size(), MPI_DOUBLE, MPI_SUM, 0, world);
#else
  array_2d<double> tmp_recv;
  if( mpigetrank() == 0) {
    for(int i=1;i<world.size();++i) {
      world.recv(i, i, tmp_recv);
      for(int k=0;k<spatial_onepdm.dim1();++k)
        for(int l=0;l<spatial_onepdm.dim2();++l)
          if( abs(tmp_recv(k,l)) > NUMERICAL_ZERO )
          {
          assert ( abs(onepdm(k,l)) < NUMERICAL_ZERO );
          onepdm(k,l) = tmp_recv(k,l);
          }
    }
  }
  else {
    world.send(0, mpigetrank(), spatial_onepdm);
  }
#endif
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Onepdm_container::calculate_spatial_npdm()
{
  //const std::vector<int>& ro = dmrginp.reorder_vector();
  //mpi::communicator world;
  if( mpigetrank() == 0) {
    for(int k=0;k<spatial_onepdm.dim1();++k)
      for(int l=0;l<spatial_onepdm.dim2();++l)
            spatial_onepdm(k,l)= onepdm(2*k,2*l)+onepdm(2*k+1,2*l+1);

  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Onepdm_container::update_full_spin_array( std::vector< std::pair< std::vector<int>, double > >& spin_batch )
{
  const std::vector<int>& ro = dmrginp.reorder_vector();
  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    int i0 = (it->first)[0];
    int j0 = (it->first)[1];
    int i= ro.at(i0/2)*2 + i0%2;
    int j= ro.at(j0/2)*2 + j0%2;

    double val = it->second;
//    if ( abs(val) > 1e-8 ) pout << "so-onepdm val: i,j = " << i << "," << j << "\t\t" << val << endl;
//    pout << "so-onepdm val: i,j = " << i << "," << j << "\t\t" << val << endl;

    // Test for duplicates
    if ( abs(val) > NUMERICAL_ZERO ) onepdm(i,j) = val;
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// This routine assumes that no spin-orbital indices are generated more than once

void Onepdm_container::update_full_spatial_array( std::vector< std::pair< std::vector<int>, double > >& spin_batch )
{
  // Take into account orbital reordering
  const std::vector<int>& ro = dmrginp.reorder_vector();

  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    assert( (it->first).size() == 2 );

    // Store significant elements only
    if ( abs(it->second) > NUMERICAL_ZERO ) {
      // Spin indices
      int i = (it->first)[0];
      int j = (it->first)[1];
      spatial_onepdm( ro.at(i), ro.at(j) ) += it->second;
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//
//void Onepdm_container::build_spatial_elements( std::map< std::vector<int>, double >& spin_batch, std::map< std::vector<int>, double >& spatial_batch )
//{
//  double factor = 1.0;
//
//  for (auto it = spatial_batch.begin(); it != spatial_batch.end(); ++it) {
//    int i = (it->first)[0];
//    int j = (it->first)[1];
//    // Sum over spin indices 
//    double val = 0.0;
//    for (int s=0; s<2; s++) {
//      std::vector<int> idx = { 2*i+s, 2*j+s };
//      val += spin_batch[ idx ];
//    }
//    // Store significant elements only
//    if ( abs(val) > NUMERICAL_ZERO ) {
//      if ( abs( spatial_onepdm(i,j) ) > NUMERICAL_ZERO ) { 
//        pout << "repeated spatial indices!\n";
//        abort();
//      }
//      spatial_onepdm(i,j) = factor * val;
//    }
//  }
//
//}
//
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Onepdm_container::store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements)
{
  Onepdm_permutations perm;
  if(dmrginp.spinAdapted())
  {
    std::vector< std::pair< std::vector<int>, double > > spatial_batch;
    perm.get_spatial_batch(new_spin_orbital_elements,spatial_batch);
    update_full_spatial_array(spatial_batch);
    if( store_full_spin_array_)
    {
      std::vector< std::pair< std::vector<int>, double > > spin_batch;
      perm.process_new_elements( new_spin_orbital_elements, nonredundant_elements, spin_batch );
      update_full_spin_array( spin_batch );
    }
  }
  else{
    std::vector< std::pair< std::vector<int>, double > > spin_batch;
    perm.process_new_elements( new_spin_orbital_elements, nonredundant_elements, spin_batch );
    update_full_spin_array( spin_batch );
  }
}

//===========================================================================================================================================================

}
}

