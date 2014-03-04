/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "onepdm_container.h"

namespace SpinAdapted{

//===========================================================================================================================================================

Onepdm_container::Onepdm_container( int sites )
{
  onepdm.resize(2*sites,2*sites);
  spatial_onepdm.resize(sites,sites);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
  
void Onepdm_container::save_npdms(const int& i, const int& j)
{
  boost::mpi::communicator world;
  world.barrier();
  Timer timer;
  accumulate_npdm();
  save_npdm_binary(i, j);
  save_npdm_text(i, j);
  accumulate_spatial_npdm();
  save_spatial_npdm_binary(i, j);
  save_spatial_npdm_text(i, j);
  world.barrier();
  pout << "1PDM save full array time " << timer.elapsedwalltime() << " " << timer.elapsedcputime() << endl;
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
    std::cout << "Spin-orbital 1PDM trace = " << trace << "\n";
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
    std::cout << "Spatial      1PDM trace = " << trace << "\n";
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
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << spatial_onepdm;
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Onepdm_container::accumulate_npdm()
{
#ifndef SERIAL
  array_2d<double> tmp_recv;
  mpi::communicator world;
  if( mpigetrank() == 0) {
    for(int i=1;i<world.size();++i) {
      world.recv(i, i, tmp_recv);
      for(int k=0;k<onepdm.dim1();++k)
        for(int l=0;l<onepdm.dim2();++l)
          if(tmp_recv(k,l) != 0.) onepdm(k,l) = tmp_recv(k,l);
    }
  }
  else {
    world.send(0, mpigetrank(), onepdm);
  }
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Onepdm_container::accumulate_spatial_npdm()
{
#ifndef SERIAL
  array_2d<double> tmp_recv;
  mpi::communicator world;
  if( mpigetrank() == 0) {
    for(int i=1;i<world.size();++i) {
      world.recv(i, i, tmp_recv);
      for(int k=0;k<spatial_onepdm.dim1();++k)
        for(int l=0;l<spatial_onepdm.dim2();++l)
          if(tmp_recv(k,l) != 0.) spatial_onepdm(k,l) = tmp_recv(k,l);
    }
  }
  else {
    world.send(0, mpigetrank(), spatial_onepdm);
  }
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Onepdm_container::update_full_spin_array( std::vector< std::pair< std::vector<int>, double > >& spin_batch )
{
  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    int i = (it->first)[0];
    int j = (it->first)[1];

    double val = it->second;
//    if ( abs(val) > 1e-8 ) pout << "so-onepdm val: i,j = " << i << "," << j << "\t\t" << val << endl;
//    pout << "so-onepdm val: i,j = " << i << "," << j << "\t\t" << val << endl;

    // Test for duplicates
    if ( onepdm(i, j) != 0.0 ) {
      cout << "WARNING: Already calculated "<<i<<" "<<j<<endl;
      cout << "earlier value: "<<onepdm(i,j)<<endl<< "new value:     "<<val<<endl;
      assert( false );
    }
    if ( abs(val) > 1e-14 ) onepdm(i,j) = val;
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// This routine assumes that no spin-orbital indices are generated more than once

void Onepdm_container::update_full_spatial_array( std::vector< std::pair< std::vector<int>, double > >& spin_batch )
{
  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    assert( (it->first).size() == 2 );

    // Store significant elements only
    if ( abs(it->second) > 1e-14 ) {
      // Spin indices
      int i = (it->first)[0];
      int j = (it->first)[1];
      if ( i%2 != j%2 ) continue;
      spatial_onepdm(i/2,j/2) += it->second;
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
//    if ( abs(val) > 1e-14 ) {
//      if ( abs( spatial_onepdm(i,j) ) > 1e-14 ) { 
//        cout << "repeated spatial indices!\n";
//        assert(false);
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
  assert( new_spin_orbital_elements.size() == 2 );
  Onepdm_permutations perm;
  std::vector< std::pair< std::vector<int>, double > > spin_batch;
  // Work with the non-redundant elements only, and get all unique spin-permutations as a by-product
  perm.process_new_elements( new_spin_orbital_elements, nonredundant_elements, spin_batch );

  update_full_spin_array( spin_batch );
  update_full_spatial_array( spin_batch );
}

//===========================================================================================================================================================

}

