/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include <algorithm>
#include "threepdm_container.h"

namespace SpinAdapted{

//===========================================================================================================================================================

Threepdm_container::Threepdm_container( int sites )
{
  bool store_full_spin_array_ = true;
  bool store_full_spatial_array_ = true;

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
//  save_sparse_array(i,j);

//  if ( use_full_array_ ) {

//FIXME
boost::mpi::communicator world;
world.barrier();
    Timer timer;
    // Combine NPDM elements from all mpi ranks and dump to file
    accumulate_npdm();
    save_npdm_text(i, j);
    save_npdm_binary(i, j);
    save_spatial_npdm_text(i, j);
    save_spatial_npdm_binary(i, j);
world.barrier();
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
//cout << "threepdm size = " << threepdm.get_size() << " ; rank = " << mpigetrank() << endl;
//cout.flush();
//threepdm.resize(26,26,26,26,26,26);
  if( mpigetrank() == 0)
  {
    for(int p=1; p<world.size(); ++p) {
//FIXME
//assert(false); // MEMORY USE LARGE HERE!!!  Is this why it crashes sometimes????
      world.recv(p, p, tmp_recv);
      for(int i=0; i<threepdm.dim1(); ++i)
        for(int j=0; j<threepdm.dim2(); ++j)
          for(int k=0; k<threepdm.dim3(); ++k)
            for(int l=0; l<threepdm.dim4(); ++l)
              for(int m=0; m<threepdm.dim5(); ++m)
                for(int n=0; n<threepdm.dim6(); ++n) 
                  if(tmp_recv(i,j,k,l,m,n) != 0.0) threepdm(i,j,k,l,m,n) = tmp_recv(i,j,k,l,m,n);
    }
  }
  else 
  {
    world.send(0, mpigetrank(), threepdm);
  }
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::update_full_spin_array()
{
  for (auto it = sparse_spin_pdm.begin(); it != sparse_spin_pdm.end(); ++it) {
    int i = (it->first)[0];
    int j = (it->first)[1];
    int k = (it->first)[2];
    int l = (it->first)[3];
    int m = (it->first)[4];
    int n = (it->first)[5];

    double val = it->second;
    //if ( abs(val) > 1e-8 ) {
    //  cout << "so-threepdm val: i,j,k,l,m,n = " 
    //       << i << "," << j << "," << k << "," << l << "," << m << "," << n
    //       << "\t\t" << val << endl;
    //}

    // Test for duplicates
    if ( threepdm(i,j,k,l,m,n) != 0.0 && abs(threepdm(i,j,k,l,m,n)-val) > 1e-6) {
      void *array[10];
      size_t size;
      size = backtrace(array, 10);
      cout << "WARNING: Already calculated "<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<m<<" "<<n<<endl;
      //backtrace_symbols_fd(array, size, 2);
      cout << "earlier value: " << threepdm(i,j,k,l,m,n) << endl << "new value:     " <<val<<endl;
      assert( false );
    }
  
    if ( abs(val) < 1e-14 ) return;
  
    // If indices are not all unique, then all elements should be zero
    std::vector<int> v = {i,j,k};
    std::sort( v.begin(), v.end() );
    if ( (v[0]==v[1]) || (v[1]==v[2]) ) { if (abs(val) > 1e-15) { std::cout << abs(val) << std::endl; assert(false); } }
    std::vector<int> w = {l,m,n};
    std::sort( w.begin(), w.end() );
    if ( (w[0]==w[1]) || (w[1]==w[2]) ) { if (abs(val) > 1e-15) { std::cout << abs(val) << std::endl; assert(false); } }

    threepdm(i,j,k,l,m,n) = val;
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::update_full_spatial_array()
{
  for (auto it = sparse_spatial_pdm.begin(); it != sparse_spatial_pdm.end(); ++it) {
    int i = (it->first)[0];
    int j = (it->first)[1];
    int k = (it->first)[2];
    int l = (it->first)[3];
    int m = (it->first)[4];
    int n = (it->first)[5];
    spatial_threepdm(i,j,k,l,m,n) = it->second;
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::build_spatial_elements( std::map< std::vector<int>, double >& spin_batch, 
                                                 std::map< std::vector<int>, double >& spatial_batch )
{
  double factor = 1.0;

  for (auto it = spatial_batch.begin(); it != spatial_batch.end(); ++it) {
    int i = (it->first)[0];
    int j = (it->first)[1];
    int k = (it->first)[2];
    int l = (it->first)[3];
    int m = (it->first)[4];
    int n = (it->first)[5];
    // Sum over spin indices
    double val = 0.0;
    for (int s=0; s<2; s++) {
      for (int t=0; t<2; t++) {
        for (int u=0; u<2; u++) {
          std::vector<int> idx = { 2*i+s, 2*j+t, 2*k+u, 2*l+u, 2*m+t, 2*n+s };
          val += spin_batch[ idx ];
        }
      }
    }
    // Store significant elements only
    if ( abs(val) > 1e-14 ) sparse_spatial_pdm[ it->first ] = factor * val;
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

std::map< std::vector<int>, int > Threepdm_container::get_spin_permutations( const std::vector<int>& indices )
{
  assert( indices.size() == 6 );
  std::map< std::vector<int>, int > perms;
  std::vector<int> idx;
  int i = indices[0];
  int j = indices[1];
  int k = indices[2];
  int l = indices[3];
  int m = indices[4];
  int n = indices[5];

  // The number of possible combinations is (3!)**2
  //------------------------------------------------

  idx = { i, j, k, l, m, n }; perms[ idx ] =  1;
  idx = { i, j, k, l, n, m }; perms[ idx ] = -1;
  idx = { i, j, k, n, m, l }; perms[ idx ] = -1;
  idx = { i, j, k, n, l, m }; perms[ idx ] =  1;
  idx = { i, j, k, m, l, n }; perms[ idx ] = -1;
  idx = { i, j, k, m, n, l }; perms[ idx ] =  1;

  idx = { i, k, j, l, m, n }; perms[ idx ] = -1;
  idx = { i, k, j, l, n, m }; perms[ idx ] =  1;
  idx = { i, k, j, n, m, l }; perms[ idx ] =  1;
  idx = { i, k, j, n, l, m }; perms[ idx ] = -1;
  idx = { i, k, j, m, l, n }; perms[ idx ] =  1;
  idx = { i, k, j, m, n, l }; perms[ idx ] = -1;

  idx = { j, i, k, l, m, n }; perms[ idx ] = -1;
  idx = { j, i, k, l, n, m }; perms[ idx ] =  1;
  idx = { j, i, k, n, m, l }; perms[ idx ] =  1;
  idx = { j, i, k, n, l, m }; perms[ idx ] = -1;
  idx = { j, i, k, m, l, n }; perms[ idx ] =  1;
  idx = { j, i, k, m, n, l }; perms[ idx ] = -1;

  idx = { j, k, i, l, m, n }; perms[ idx ] =  1;
  idx = { j, k, i, l, n, m }; perms[ idx ] = -1;
  idx = { j, k, i, n, m, l }; perms[ idx ] = -1;
  idx = { j, k, i, n, l, m }; perms[ idx ] =  1;
  idx = { j, k, i, m, l, n }; perms[ idx ] = -1;
  idx = { j, k, i, m, n, l }; perms[ idx ] =  1;

  idx = { k, j, i, l, m, n }; perms[ idx ] = -1;
  idx = { k, j, i, l, n, m }; perms[ idx ] =  1;
  idx = { k, j, i, n, m, l }; perms[ idx ] =  1;
  idx = { k, j, i, n, l, m }; perms[ idx ] = -1;
  idx = { k, j, i, m, l, n }; perms[ idx ] =  1;
  idx = { k, j, i, m, n, l }; perms[ idx ] = -1;

  idx = { k, i, j, l, m, n }; perms[ idx ] =  1;
  idx = { k, i, j, l, n, m }; perms[ idx ] = -1;
  idx = { k, i, j, n, m, l }; perms[ idx ] = -1;
  idx = { k, i, j, n, l, m }; perms[ idx ] =  1;
  idx = { k, i, j, m, l, n }; perms[ idx ] = -1;
  idx = { k, i, j, m, n, l }; perms[ idx ] =  1;

  // Get transpose elements with same parity factors
  std::map< std::vector<int>, int > trans_perms;
  for (auto it = perms.begin(); it != perms.end(); ++it) {
    std::vector<int> indices = it->first;
    std::reverse( indices.begin(), indices.end() );
    trans_perms[ indices ] = it->second;
  }

  // Now bundle them togther
  for (auto it = trans_perms.begin(); it != trans_perms.end(); ++it) {
    perms[ it->first ] = it->second;
  }

  return perms;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements)
{
  assert( new_spin_orbital_elements.size() == 20 );
  // Temporary batches of npdm elements
  std::map< std::vector<int>, double > spin_batch;
  std::map< std::vector<int>, double > spatial_batch;

  for (int idx=0; idx < new_spin_orbital_elements.size(); ++idx) {
    // Get all spin-index permutations
    std::map< std::vector<int>, int > spin_indices = get_spin_permutations( new_spin_orbital_elements[idx].first );
    double val = new_spin_orbital_elements[idx].second;
    for (auto it = spin_indices.begin(); it != spin_indices.end(); ++it) {
      // Initialize spatial indices
      std::vector<int> vec;
      for (int i=0; i < (it->first).size(); ++i)
        vec.push_back( (it->first)[i]/2 );
      spatial_batch[ vec ] = 0.0;
      // Assign temporary batch of spin-orbital elements
      spin_batch[ it->first ] = it->second * val;
      // Store significant elements only
      if ( abs(val) > 1e-14 ) sparse_spin_pdm[ it->first ] = it->second * val;
    }
  }

  // Build and store new spatial elements
  build_spatial_elements( spin_batch, spatial_batch );

}

//===========================================================================================================================================================

}


