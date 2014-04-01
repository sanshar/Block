/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "global.h"
#include <boost/format.hpp>
#include "threepdm_container.h"
#include "npdm_permutations.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

Threepdm_container::Threepdm_container( int sites )
{
  store_nonredundant_spin_elements_ = true;
  store_full_spin_array_ = false;
  store_full_spatial_array_ = true;

  if ( store_full_spin_array_ ) {
    threepdm.resize(2*sites,2*sites,2*sites,2*sites,2*sites,2*sites);
    threepdm.Clear();
  } 
  if ( store_full_spatial_array_ ) {
    spatial_threepdm.resize(sites,sites,sites,sites,sites,sites);
    spatial_threepdm.Clear();
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
  if ( store_full_spin_array_ ) {
    accumulate_npdm();
    save_npdm_binary(i, j);
    save_npdm_text(i, j);
  }
  if ( store_full_spatial_array_ ) {
    accumulate_spatial_npdm();
    save_spatial_npdm_binary(i, j);
    save_spatial_npdm_text(i, j);
  }

#ifndef SERIAL
  world.barrier();
#endif
  pout << "3PDM save full array time " << timer.elapsedwalltime() << " " << timer.elapsedcputime() << endl;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::save_npdm_text(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/threepdm.", i, j,".txt");
    ofstream ofs(file);
    ofs << threepdm.dim1() << endl;

    double trace = 0.0;
    for(int i=0; i<threepdm.dim1(); ++i)
      for(int j=0; j<threepdm.dim2(); ++j)
        for(int k=0; k<threepdm.dim3(); ++k)
          for(int l=0; l<threepdm.dim4(); ++l)
            for(int m=0; m<threepdm.dim5(); ++m)
              for(int n=0; n<threepdm.dim6(); ++n) {
                if ( abs(threepdm(i,j,k,l,m,n)) > 1e-14 ) {
                  ofs << boost::format("%d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % threepdm(i,j,k,l,m,n);
                  if ( (i==n) && (j==m) && (k==l) ) trace += threepdm(i,j,k,l,m,n);
                }
              }
    ofs.close();
    std::cout << "Spin-orbital 3PDM trace = " << trace << "\n";
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::save_spatial_npdm_text(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.", i, j,".txt");
    ofstream ofs(file);
    ofs << spatial_threepdm.dim1() << endl;

    double trace = 0.0;
    for(int i=0; i<spatial_threepdm.dim1(); ++i)
      for(int j=0; j<spatial_threepdm.dim2(); ++j)
        for(int k=0; k<spatial_threepdm.dim3(); ++k)
          for(int l=0; l<spatial_threepdm.dim4(); ++l)
            for(int m=0; m<spatial_threepdm.dim5(); ++m)
              for(int n=0; n<spatial_threepdm.dim6(); ++n) {
                if ( abs(spatial_threepdm(i,j,k,l,m,n)) > 1e-14 ) {
                  ofs << boost::format("%d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % spatial_threepdm(i,j,k,l,m,n);
                  if ( (i==n) && (j==m) && (k==l) ) trace += spatial_threepdm(i,j,k,l,m,n);
                }
              }
    ofs.close();
    std::cout << "Spatial      3PDM trace = " << trace << "\n";
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::save_npdm_binary(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/threepdm.", i, j,".bin");
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << threepdm;
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::save_spatial_npdm_binary(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.", i, j,".bin");
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << spatial_threepdm;
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::load_npdm_binary(const int &i, const int &j) { assert(false); }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::accumulate_npdm()
{
#ifndef SERIAL
  array_6d<double> tmp_recv;
  mpi::communicator world;
  if( mpigetrank() == 0)
  {
    for(int p=1; p<world.size(); ++p) {
      world.recv(p, p, tmp_recv);
      for(int i=0; i<threepdm.dim1(); ++i)
        for(int j=0; j<threepdm.dim2(); ++j)
          for(int k=0; k<threepdm.dim3(); ++k)
            for(int l=0; l<threepdm.dim4(); ++l)
              for(int m=0; m<threepdm.dim5(); ++m)
                for(int n=0; n<threepdm.dim6(); ++n) {
                  if ( abs(tmp_recv(i,j,k,l,m,n)) > 1e-15 ) {
                    if ( abs(threepdm(i,j,k,l,m,n)) > 1e-14 ) assert(false);
                    threepdm(i,j,k,l,m,n) = tmp_recv(i,j,k,l,m,n);
                  }
                }
    }
  }
  else 
  {
    world.send(0, mpigetrank(), threepdm);
  }
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::accumulate_spatial_npdm()
{
#ifndef SERIAL
  array_6d<double> tmp_recv;
  mpi::communicator world;
  if( mpigetrank() == 0)
  {
    for(int p=1; p<world.size(); ++p) {
      world.recv(p, p, tmp_recv);
      for(int i=0; i<spatial_threepdm.dim1(); ++i)
        for(int j=0; j<spatial_threepdm.dim2(); ++j)
          for(int k=0; k<spatial_threepdm.dim3(); ++k)
            for(int l=0; l<spatial_threepdm.dim4(); ++l)
              for(int m=0; m<spatial_threepdm.dim5(); ++m)
                for(int n=0; n<spatial_threepdm.dim6(); ++n) {
                  if( abs(tmp_recv(i,j,k,l,m,n)) > 1e-15 ) {
                    // Test if any duplicate elements built on different processors
                    if ( abs(spatial_threepdm(i,j,k,l,m,n)) > 1e-14 ) assert(false);
                    spatial_threepdm(i,j,k,l,m,n) = tmp_recv(i,j,k,l,m,n);
                  }
                }
    }
  }
  else 
  {
    world.send(0, mpigetrank(), spatial_threepdm);
  }
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::update_full_spin_array( std::vector< std::pair< std::vector<int>, double > >& spin_batch )
{

  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    double val = it->second;
    if ( abs(val) < 1e-15 ) continue;

    assert( (it->first).size() == 6 );
    int i = (it->first)[0];
    int j = (it->first)[1];
    int k = (it->first)[2];
    int l = (it->first)[3];
    int m = (it->first)[4];
    int n = (it->first)[5];

    //if ( abs(val) > 1e-8 ) {
    //  cout << "so-threepdm val: i,j,k,l,m,n = " 
    //       << i << "," << j << "," << k << "," << l << "," << m << "," << n
    //       << "\t\t" << val << endl;
    //}

    // Test for duplicates
    if ( abs(threepdm(i,j,k,l,m,n)) != 0.0 ) {
      cout << "WARNING: Already calculated "<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<m<<" "<<n<<endl;
      cout << "earlier value: " << threepdm(i,j,k,l,m,n) << endl << "new value:     " <<val<<endl;
      assert( false );
    }
    threepdm(i,j,k,l,m,n) = val;
  }

}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// This routine assumes that no spin-orbital indices are generated more than once

void Threepdm_container::update_full_spatial_array( std::vector< std::pair< std::vector<int>, double > >& spin_batch )
{
  // Take into account orbital reordering
  const std::vector<int>& ro = dmrginp.reorder_vector();

  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    assert( (it->first).size() == 6 );

    // Store significant elements only
    if ( abs(it->second) > 1e-14 ) {
      // Spin indices
      int i = (it->first)[0];
      int j = (it->first)[1];
      int k = (it->first)[2];
      int l = (it->first)[3];
      int m = (it->first)[4];
      int n = (it->first)[5];

      if ( i%2 != n%2 ) continue;
      if ( j%2 != m%2 ) continue;
      if ( k%2 != l%2 ) continue;

      spatial_threepdm( ro.at(i/2), ro.at(j/2), ro.at(k/2), ro.at(l/2), ro.at(m/2), ro.at(n/2) ) += it->second;
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//
//void Threepdm_container::build_spatial_elements( std::map< std::vector<int>, double >& spin_batch, 
//                                                 std::map< std::vector<int>, double >& spatial_batch )
//{
//  double factor = 1.0;
//
//  for (auto it = spatial_batch.begin(); it != spatial_batch.end(); ++it) {
//    assert( (it->first).size() == 6 );
//    int i = (it->first)[0];
//    int j = (it->first)[1];
//    int k = (it->first)[2];
//    int l = (it->first)[3];
//    int m = (it->first)[4];
//    int n = (it->first)[5];
//    // Sum over spin indices
//    double val = 0.0;
//    for (int s=0; s<2; s++) {
//      for (int t=0; t<2; t++) {
//        for (int u=0; u<2; u++) {
//          std::vector<int> idx = { 2*i+s, 2*j+t, 2*k+u, 2*l+u, 2*m+t, 2*n+s };
//          val += spin_batch[ idx ];
//        }
//      }
//    }
//    // Store significant elements only
//    if ( abs(val) > 1e-14 ) {
//      //cout << "i,j,k,l,m,n = " << i << "," << j << "," << k << "," << l << "," << m << "," << n << "\t\t" << val << endl;
//      if ( store_sparse_spatial_array_ ) sparse_spatial_pdm[ it->first ] = factor * val;
//      if ( store_full_spatial_array_ ) {
//        if ( abs( spatial_threepdm(i,j,k,l,m,n) ) > 1e-14 ) {
//          cout << "repeated spatial indices!\n";
//          assert(false);
//        }
//        spatial_threepdm(i,j,k,l,m,n) = factor * val;
//      }
//    }
//  }
//
//}
//
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements)
{
  assert( new_spin_orbital_elements.size() == 20 );
  Threepdm_permutations perm;
  std::vector< std::pair< std::vector<int>, double > > spin_batch;
  // Work with the non-redundant elements only, and get all unique spin-permutations as a by-product
  perm.process_new_elements( new_spin_orbital_elements, nonredundant_elements, spin_batch );

  //FIXME add options to dump to disk if memory becomes bottleneck
  if ( ! store_nonredundant_spin_elements_ ) nonredundant_elements.clear();
  if ( store_full_spin_array_ ) update_full_spin_array( spin_batch );
  if ( store_full_spatial_array_ ) update_full_spatial_array( spin_batch );
}

//===========================================================================================================================================================

}
}


