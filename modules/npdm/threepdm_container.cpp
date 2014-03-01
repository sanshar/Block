/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "threepdm_container.h"

//FIXME parallel
namespace SpinAdapted{

//===========================================================================================================================================================

Threepdm_container::Threepdm_container( int sites )
{
  store_full_spin_array_ = false;
  store_full_spatial_array_ = true;
  store_sparse_spin_array_ = false;
  store_sparse_spatial_array_ = true;

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
  boost::mpi::communicator world;
  world.barrier();
  Timer timer;
  if ( store_full_spin_array_ ) {
    accumulate_npdm();
    save_npdm_text(i, j);
    save_npdm_binary(i, j);
  }
  if ( store_full_spatial_array_ ) {
    accumulate_spatial_npdm();
    save_spatial_npdm_text(i, j);
    save_spatial_npdm_binary(i, j);
  }

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
  if( mpigetrank() == 0)
  {
    for(int p=1; p<world.size(); ++p) {
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
                for(int n=0; n<spatial_threepdm.dim6(); ++n) 
                  if(tmp_recv(i,j,k,l,m,n) != 0.0) spatial_threepdm(i,j,k,l,m,n) = tmp_recv(i,j,k,l,m,n);
    }
  }
  else 
  {
    world.send(0, mpigetrank(), spatial_threepdm);
  }
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::update_full_spin_array( std::map< std::vector<int>, double >& spin_batch )
{

  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    assert( (it->first).size() == 6 );
    int i = (it->first)[0];
    int j = (it->first)[1];
    int k = (it->first)[2];
    int l = (it->first)[3];
    int m = (it->first)[4];
    int n = (it->first)[5];

    double val = it->second;
//    if ( abs(val) > 1e-8 ) {
//      cout << "so-threepdm val: i,j,k,l,m,n = " 
//           << i << "," << j << "," << k << "," << l << "," << m << "," << n
//           << "\t\t" << val << endl;
//    }

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
  
    if ( abs(val) < 1e-14 ) continue;
  
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

void Threepdm_container::build_spatial_elements( std::map< std::vector<int>, double >& spin_batch, 
                                                 std::map< std::vector<int>, double >& spatial_batch )
{
  double factor = 1.0;

  for (auto it = spatial_batch.begin(); it != spatial_batch.end(); ++it) {
    assert( (it->first).size() == 6 );
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
    if ( abs(val) > 1e-14 ) {
      //cout << "i,j,k,l,m,n = " << i << "," << j << "," << k << "," << l << "," << m << "," << n << "\t\t" << val << endl;
      if ( store_sparse_spatial_array_ ) sparse_spatial_pdm[ it->first ] = factor * val;
      if ( store_full_spatial_array_ ) {
        if ( abs( spatial_threepdm(i,j,k,l,m,n) ) > 1e-14 ) {
          cout << "repeated spatial indices!\n";
          assert(false);
        }
        spatial_threepdm(i,j,k,l,m,n) = factor * val;
      }
    }
  }

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
    Threepdm_permutations p;
    std::map< std::vector<int>, int > spin_indices = p.get_spin_permutations( new_spin_orbital_elements[idx].first );
    double val = new_spin_orbital_elements[idx].second;
    for (auto it = spin_indices.begin(); it != spin_indices.end(); ++it) {
      // Initialize spatial indices
      std::vector<int> vec;
      for (int i=0; i < (it->first).size(); ++i)
        vec.push_back( (it->first)[i]/2 );
      spatial_batch[ vec ] = 0.0;
      // Assign temporary batch of spin-orbital elements
      spin_batch[ it->first ] = it->second * val;
      if ( store_sparse_spin_array_ && (abs(val) > 1e-14) ) sparse_spin_pdm[ it->first ] = it->second * val;
    }
  }

  // Build and store new spatial elements
  build_spatial_elements( spin_batch, spatial_batch );
  if ( store_full_spin_array_ ) update_full_spin_array( spin_batch );
}

//===========================================================================================================================================================

}


