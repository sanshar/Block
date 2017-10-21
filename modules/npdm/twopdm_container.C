/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "global.h"
#include <boost/format.hpp>
#include "twopdm_container.h"
#include "npdm_permutations.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

Twopdm_container::Twopdm_container( int sites )
{

  if ( store_full_spin_array_ ) {
    if(dmrginp.spinAdapted())
      twopdm.resize(2*sites,2*sites,2*sites,2*sites);
    else
      twopdm.resize(sites,sites,sites,sites);
    twopdm.Clear();
  }
  if ( store_full_spatial_array_ ) {
    if(dmrginp.spinAdapted())
      spatial_twopdm.resize(sites,sites,sites,sites);
    else
      spatial_twopdm.resize(sites/2,sites/2,sites/2,sites/2);
    spatial_twopdm.Clear();
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
  if ( store_full_spin_array_ ) {
    accumulate_npdm();
    save_npdm_binary(i, j);
    save_npdm_text(i, j);
  }
  if ( store_full_spatial_array_ ) {
    if(dmrginp.spinAdapted())
      accumulate_spatial_npdm();
    else
      calculate_spatial_npdm();
    save_spatial_npdm_binary(i, j);
    save_spatial_npdm_text(i, j);
  }

#ifndef SERIAL
  world.barrier();
#endif
  ecpu = timer.elapsedcputime();ewall=timer.elapsedwalltime();
  p3out << "2PDM save full array time " << ewall << " " << ecpu << endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_container::save_npdm_text(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/twopdm.", i, j,".txt");
    ofstream ofs(file);
    ofs << twopdm.dim1() << endl;
    double trace = 0.0;
    for(int k=0; k<twopdm.dim1(); ++k)
      for(int l=0; l<twopdm.dim2(); ++l)
        for(int m=0; m<twopdm.dim3(); ++m)
          for(int n=0; n<twopdm.dim4(); ++n) {
            ofs << boost::format("%d %d %d %d %20.14e\n") % k % l % m % n % twopdm(k,l,m,n);
            if ( (k==n) && (l==m) ) trace += twopdm(k,l,m,n);
          }
    ofs.close();
    pout << "Spin-orbital 2PDM trace = " << trace << "\n";
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_container::save_spatial_npdm_text(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_twopdm.", i, j,".txt");
    ofstream ofs(file);
    ofs << spatial_twopdm.dim1() << endl;
    double trace = 0.0;
    for(int k=0; k<spatial_twopdm.dim1(); ++k)
      for(int l=0; l<spatial_twopdm.dim2(); ++l)
        for(int m=0; m<spatial_twopdm.dim3(); ++m)
          for(int n=0; n<spatial_twopdm.dim4(); ++n) {
            ofs << boost::format("%d %d %d %d %20.14e\n") % k % l % m % n % spatial_twopdm(k,l,m,n);
            if ( (k==n) && (l==m) ) trace += spatial_twopdm(k,l,m,n);
          }
    ofs.close();
    pout << "Spatial      2PDM trace = " << trace << "\n";
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_container::save_npdm_binary(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
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
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_twopdm.", i, j,".bin");
#ifndef MOLCAS
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << spatial_twopdm;
    ofs.close();
#else
    FILE* f = fopen(file,"wb");
    fwrite(spatial_twopdm.data(),sizeof(double),spatial_twopdm.size(),f);
    fclose(f);
#endif
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_container::accumulate_npdm()
{
#ifndef SERIAL
  mpi::communicator world;
#ifdef NDEBUG
  if( mpigetrank() == 0)
    MPI_Reduce(MPI_IN_PLACE, twopdm.data(),  twopdm.size(), MPI_DOUBLE, MPI_SUM, 0, world);
  else
    MPI_Reduce(twopdm.data(), twopdm.data(),  twopdm.size(), MPI_DOUBLE, MPI_SUM, 0, world);
#else
  array_4d<double> tmp_recv;
  if( mpigetrank() == 0) {
    for(int i=1;i<world.size();++i) {
      world.recv(i, i, tmp_recv);
      for(int k=0;k<twopdm.dim1();++k)
        for(int l=0;l<twopdm.dim2();++l)
          for(int m=0;m<twopdm.dim3();++m)
            for(int n=0;n<twopdm.dim4();++n)
              if ( abs(tmp_recv(k,l,m,n)) > NUMERICAL_ZERO) {
                // Test for duplicates
                assert(abs(twopdm(k,l,m,n)) < NUMERICAL_ZERO ) ;
                twopdm(k,l,m,n) = tmp_recv(k,l,m,n);
              }
	 }
  }
  else {
    world.send(0, mpigetrank(), twopdm);
  }
#endif
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_container::accumulate_spatial_npdm()
{
#ifndef SERIAL
  mpi::communicator world;
#ifdef NDEBUG
  if( mpigetrank() == 0)
    MPI_Reduce(MPI_IN_PLACE, spatial_twopdm.data(),  spatial_twopdm.size(), MPI_DOUBLE, MPI_SUM, 0, world);
  else
    MPI_Reduce(spatial_twopdm.data(), spatial_twopdm.data(),  spatial_twopdm.size(), MPI_DOUBLE, MPI_SUM, 0, world);
#else
  array_4d<double> tmp_recv;
  if( mpigetrank() == 0) {
    for(int i=1;i<world.size();++i) {
      world.recv(i, i, tmp_recv);
      for(int k=0;k<spatial_twopdm.dim1();++k)
        for(int l=0;l<spatial_twopdm.dim2();++l)
          for(int m=0;m<spatial_twopdm.dim3();++m)
            for(int n=0;n<spatial_twopdm.dim4();++n)
              if ( abs(tmp_recv(k,l,m,n)) > NUMERICAL_ZERO) {
                // Test for duplicates
                assert(abs(spatial_twopdm(k,l,m,n)) < NUMERICAL_ZERO ) ;
                spatial_twopdm(k,l,m,n) = tmp_recv(k,l,m,n);
              }
	 }
  }
  else {
    world.send(0, mpigetrank(), spatial_twopdm);
  }
#endif
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_container::update_full_spin_array( std::vector< std::pair< std::vector<int>, double > >& spin_batch )
{
  // Take into account orbital reordering
  const std::vector<int>& ro = dmrginp.reorder_vector();
  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    double val = it->second;
    if ( abs(val) < NUMERICAL_ZERO ) continue;
    int i0 = (it->first)[0];
    int j0 = (it->first)[1];
    int k0 = (it->first)[2];
    int l0 = (it->first)[3];
    int i = ro.at(i0/2)*2 + i0%2;
    int j = ro.at(j0/2)*2 + j0%2;
    int k = ro.at(k0/2)*2 + k0%2;
    int l = ro.at(l0/2)*2 + l0%2;

    //if ( abs(val) > 1e-8 ) pout << "so-twopdm val: i,j,k,l = " << i << "," << j << "," << k << "," << l << "\t\t" << val << endl;
    //pout << "so-twopdm val: i,j,k,l = " << i << "," << j << "," << k << "," << l << "\t\t" << val << endl;

    // Test for duplicates
    twopdm(i,j,k,l) = val;
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_container::calculate_spatial_npdm()
{
  double factor = 0.5;
  //mpi::communicator world;
  if( mpigetrank() == 0) {
    for(int k=0;k<spatial_twopdm.dim1();++k)
      for(int l=0;l<spatial_twopdm.dim2();++l)
        for(int m=0;m<spatial_twopdm.dim3();++m)
          for(int n=0;n<spatial_twopdm.dim4();++n)
            spatial_twopdm(k,l,m,n)= factor*(twopdm(2*k,2*l,2*m,2*n)+ twopdm(2*k+1,2*l,2*m,2*n+1)+ twopdm(2*k,2*l+1,2*m+1,2*n)+ twopdm(2*k+1,2*l+1,2*m+1,2*n+1));

  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// This routine assumes that no spin-orbital indices are generated more than once

void Twopdm_container::update_full_spatial_array( std::vector< std::pair< std::vector<int>, double > >& spin_batch )
{
  // Take into account orbital reordering
  const std::vector<int>& ro = dmrginp.reorder_vector();

  // Note we multiply the spatial 2PDM by a factor of 1/2 to be consistent with the old BLOCK code, but this doesn't seem conventional?
  double factor = 0.5;
  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    assert( (it->first).size() == 4 );

    // Store significant elements only
    if ( abs(it->second) > NUMERICAL_ZERO ) {
      // Spin indices
      int i = (it->first)[0];
      int j = (it->first)[1];
      int k = (it->first)[2];
      int l = (it->first)[3];

      spatial_twopdm( ro.at(i), ro.at(j), ro.at(k), ro.at(l) ) += factor * (it->second);
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_container::store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements)
{
  Twopdm_permutations perm;
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

