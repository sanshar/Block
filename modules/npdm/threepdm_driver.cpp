/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include <algorithm>
#include "threepdm_driver.h"

namespace SpinAdapted{

//===========================================================================================================================================================

Threepdm_driver::Threepdm_driver( int sites ) : Npdm_driver(3) 
{
  store_full_spin_array_ = true;
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

void Threepdm_driver::save_npdms(const int& i, const int& j)
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

void Threepdm_driver::save_npdm_text(const int &i, const int &j)
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

void Threepdm_driver::save_spatial_npdm_text(const int &i, const int &j)
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

void Threepdm_driver::save_npdm_binary(const int &i, const int &j)
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

void Threepdm_driver::save_spatial_npdm_binary(const int &i, const int &j)
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

void Threepdm_driver::load_npdm_binary(const int &i, const int &j) { assert(false); }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_driver::accumulate_npdm()
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

void Threepdm_driver::update_full_spin_array()
{
  for (auto it = spin_map.begin(); it != spin_map.end(); ++it) {
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

void Threepdm_driver::update_full_spatial_array()
{
  for (auto it = spatial_map.begin(); it != spatial_map.end(); ++it) {
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

void Threepdm_driver::build_spatial_elements()
{
//FIXME parallelization!
//  if( mpigetrank() == 0)
//  {
    double factor = 1.0;

    for (auto it = spatial_map.begin(); it != spatial_map.end(); ++it) {
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
            val += spin_map[ idx ];
          }
        }
      }
      it->second = factor * val;
    }
// }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

std::map< std::vector<int>, int > Threepdm_driver::get_spin_permutations( std::vector<int>& indices )
{

  std::map< std::vector<int>, int > perms;
  std::vector<int> idx;
  int i = indices[0];
  int j = indices[1];
  int k = indices[2];
  int l = indices[3];
  int m = indices[4];
  int n = indices[5];

  // The number of possible combinations is (3!)**2 before transpose
  //-----------------------------------------------------------------

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

void Threepdm_driver::screen_sparse_arrays()
{
  // Screen spatial elements
  for (auto it = spatial_map.begin(); it != spatial_map.end(); ) {
    if ( abs(it->second) < 1e-14 ) {
      spatial_map.erase(it++); // Note postfix operator
    }
    else { ++it; }
  }

  // Screen spin-orbital elements
  for (auto it = spin_map.begin(); it != spin_map.end(); ) {
    if ( abs(it->second) < 1e-14 ) {
      spin_map.erase(it++); // Note postfix operator
    }
    else { ++it; }
  }
//FIXME add prints to see if it makes any difference?

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_driver::store_npdm_elements( std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements)
{
  assert( new_spin_orbital_elements.size() == 20 );

  for (int idx=0; idx < new_spin_orbital_elements.size(); ++idx) {
    // Get all spin-index permutations
    std::map< std::vector<int>, int > spin_indices = get_spin_permutations( new_spin_orbital_elements[idx].first );
    // Assign batch of spin-orbital 3PDM values
    double val = new_spin_orbital_elements[idx].second;
    for (auto it = spin_indices.begin(); it != spin_indices.end(); ++it) {
      spin_map[ it->first ] = it->second * val;
    }
    // Initialize spatial 3PDM indices
    for (auto it = spin_indices.begin(); it != spin_indices.end(); ++it) {
      std::vector<int> idx;
      for (int i=0; i < (it->first).size(); ++i)
         idx.push_back( (it->first)[i]/2 );
      spatial_map[ idx ] = 0.0;
    }
  }

  // Build and store spatial 3PDM elements at this sweep position
  build_spatial_elements();

  // Screen away small or zero elements
  screen_sparse_arrays();
}

//===========================================================================================================================================================

}


