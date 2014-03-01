/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "fourpdm_container.h"

namespace SpinAdapted{

//===========================================================================================================================================================

Fourpdm_container::Fourpdm_container( int sites )
{
  store_full_spin_array_ = false;
  store_full_spatial_array_ = true;
  store_sparse_spin_array_ = false;
  store_sparse_spatial_array_ = true;

  if ( store_full_spin_array_ ) {
    fourpdm.resize(2*sites,2*sites,2*sites,2*sites,2*sites,2*sites,2*sites,2*sites);
    fourpdm.Clear();
  } 
  if ( store_full_spatial_array_ ) {
    spatial_fourpdm.resize(sites,sites,sites,sites,sites,sites,sites,sites);
    spatial_fourpdm.Clear();
  } 

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::save_npdms(const int& i, const int& j)
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
  pout << "4PDM save full array time " << timer.elapsedwalltime() << " " << timer.elapsedcputime() << endl;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::save_npdm_text(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/fourpdm.", i, j,".txt");
    ofstream ofs(file);
    ofs << fourpdm.dim1() << endl;

    double trace = 0.0;
    for(int i=0; i<fourpdm.dim1(); ++i)
      for(int j=0; j<fourpdm.dim2(); ++j)
        for(int k=0; k<fourpdm.dim3(); ++k)
          for(int l=0; l<fourpdm.dim4(); ++l)
            for(int m=0; m<fourpdm.dim5(); ++m)
              for(int n=0; n<fourpdm.dim6(); ++n)
                for(int p=0; p<fourpdm.dim7(); ++p)
                  for(int q=0; q<fourpdm.dim8(); ++q) {
                    if ( abs(fourpdm(i,j,k,l,m,n,p,q)) > 1e-14 ) {
                      ofs << boost::format("%d %d %d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % p % q % fourpdm(i,j,k,l,m,n,p,q);
                      if ( (i==q) && (j==p) && (k==n) && (l==m) ) trace += fourpdm(i,j,k,l,m,n,p,q);
                    }
                  }
    ofs.close();
    std::cout << "Spin-orbital 4PDM trace = " << trace << "\n";
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::save_spatial_npdm_text(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_fourpdm.", i, j,".txt");
    ofstream ofs(file);
    ofs << spatial_fourpdm.dim1() << endl;

    double trace = 0.0;
    for(int i=0; i<spatial_fourpdm.dim1(); ++i)
      for(int j=0; j<spatial_fourpdm.dim2(); ++j)
        for(int k=0; k<spatial_fourpdm.dim3(); ++k)
          for(int l=0; l<spatial_fourpdm.dim4(); ++l)
            for(int m=0; m<spatial_fourpdm.dim5(); ++m)
              for(int n=0; n<spatial_fourpdm.dim6(); ++n)
                for(int p=0; p<spatial_fourpdm.dim7(); ++p)
                  for(int q=0; q<spatial_fourpdm.dim8(); ++q) {
                    if ( abs(spatial_fourpdm(i,j,k,l,m,n,p,q)) > 1e-14 ) {
                      ofs << boost::format("%d %d %d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % p % q % spatial_fourpdm(i,j,k,l,m,n,p,q);
                      if ( (i==q) && (j==p) && (k==n) && (l==m) ) trace += spatial_fourpdm(i,j,k,l,m,n,p,q);
                    }
                  }
    ofs.close();
    std::cout << "Spatial      4PDM trace = " << trace << "\n";
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::save_npdm_binary(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/fourpdm.", i, j,".bin");
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << fourpdm;
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::save_spatial_npdm_binary(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_fourpdm.", i, j,".bin");
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << spatial_fourpdm;
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::load_npdm_binary(const int &i, const int &j) { assert(false); }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::accumulate_npdm()
{
#ifndef SERIAL
  array_8d<double> tmp_recv;
  mpi::communicator world;
  if( mpigetrank() == 0)
  {
    for(int u=1; u<world.size(); ++u) {
      world.recv(u, u, tmp_recv);
      for(int i=0; i<fourpdm.dim1(); ++i)
        for(int j=0; j<fourpdm.dim2(); ++j)
          for(int k=0; k<fourpdm.dim3(); ++k)
            for(int l=0; l<fourpdm.dim4(); ++l)
              for(int m=0; m<fourpdm.dim5(); ++m)
                for(int n=0; n<fourpdm.dim6(); ++n) 
                  for(int p=0; p<fourpdm.dim7(); ++p) 
                    for(int q=0; q<fourpdm.dim8(); ++q) {
                      if ( abs(tmp_recv(i,j,k,l,m,n,p,q)) > 1e-15 ) {
                        if ( abs(fourpdm(i,j,k,l,m,n,p,q)) > 1e-14 ) assert(false);
                        fourpdm(i,j,k,l,m,n,p,q) = tmp_recv(i,j,k,l,m,n,p,q);
                      }
                    }
    }
  }
  else 
  {
    world.send(0, mpigetrank(), fourpdm);
  }
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::accumulate_spatial_npdm()
{
#ifndef SERIAL
  array_8d<double> tmp_recv;
  mpi::communicator world;
  if( mpigetrank() == 0)
  {
    for(int u=1; u<world.size(); ++u) {
      world.recv(u, u, tmp_recv);
      for(int i=0; i<spatial_fourpdm.dim1(); ++i)
        for(int j=0; j<spatial_fourpdm.dim2(); ++j)
          for(int k=0; k<spatial_fourpdm.dim3(); ++k)
            for(int l=0; l<spatial_fourpdm.dim4(); ++l)
              for(int m=0; m<spatial_fourpdm.dim5(); ++m)
                for(int n=0; n<spatial_fourpdm.dim6(); ++n) 
                  for(int p=0; p<spatial_fourpdm.dim7(); ++p) 
                    for(int q=0; q<spatial_fourpdm.dim8(); ++q) {
                      if ( abs(tmp_recv(i,j,k,l,m,n,p,q)) > 1e-15 ) {
                        if ( abs(spatial_fourpdm(i,j,k,l,m,n,p,q)) > 1e-14 ) assert(false);
                        spatial_fourpdm(i,j,k,l,m,n,p,q) = tmp_recv(i,j,k,l,m,n,p,q);
                      }
                    }
    }
  }
  else 
  {
    world.send(0, mpigetrank(), spatial_fourpdm);
  }
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::update_full_spin_array( std::map< std::vector<int>, double >& spin_batch )
{
  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    assert( (it->first).size() == 8 );
    int i = (it->first)[0];
    int j = (it->first)[1];
    int k = (it->first)[2];
    int l = (it->first)[3];
    int m = (it->first)[4];
    int n = (it->first)[5];
    int p = (it->first)[6];
    int q = (it->first)[7];

    double val = it->second;

//if ( abs(val) > 1e-8 ) {
//  pout << "so-fourpdm val: i,j,k,l,m,n,p,q = " 
//       << i << "," << j << "," << k << "," << l << "," << m << "," << n << "," << p << "," << q
//       << "\t\t" << val << endl;
//}

    // Test for duplicates
    if ( fourpdm(i,j,k,l,m,n,p,q) != 0.0 && abs(fourpdm(i,j,k,l,m,n,p,q)-val) > 1e-6) {
      void *array[10];
      size_t size;
      size = backtrace(array, 10);
      cout << "WARNING: Already calculated "<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<m<<" "<<n<<" "<<p<<" "<<q<<endl;
      //backtrace_symbols_fd(array, size, 2);
      cout << "earlier value: " << fourpdm(i,j,k,l,m,n,p,q) << endl << "new value:     " <<val<<endl;
      cout.flush();
      assert( false );
    }
  
    if ( abs(val) < 1e-14 ) continue;
  
    // If indices are not all unique, then all elements should be zero (and next_even_permutation fails)
    std::vector<int> v = {i,j,k,l};
    std::sort( v.begin(), v.end() );
    if ( (v[0]==v[1]) || (v[1]==v[2]) || (v[2]==v[3]) ) { if (abs(val) > 1e-15) { std::cout << abs(val) << std::endl; assert(false); } }
    std::vector<int> w = {m,n,p,q};
    std::sort( w.begin(), w.end() );
    if ( (w[0]==w[1]) || (w[1]==w[2]) || (w[2]==w[3]) ) { if (abs(val) > 1e-15) { std::cout << abs(val) << std::endl; assert(false); } }

    fourpdm(i,j,k,l,m,n,p,q) = val;
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::build_spatial_elements( std::map< std::vector<int>, double >& spin_batch, 
                                                std::map< std::vector<int>, double >& spatial_batch )
{
  double factor = 1.0;

  for (auto it = spatial_batch.begin(); it != spatial_batch.end(); ++it) {
    assert( (it->first).size() == 8 );
    int i = (it->first)[0];
    int j = (it->first)[1];
    int k = (it->first)[2];
    int l = (it->first)[3];
    int m = (it->first)[4];
    int n = (it->first)[5];
    int p = (it->first)[6];
    int q = (it->first)[7];
    // Sum over spin indices
    double val = 0.0;
    for (int s=0; s<2; s++) {
      for (int t=0; t<2; t++) {
        for (int u=0; u<2; u++) {
          for (int v=0; v<2; v++) {
            std::vector<int> idx = { 2*i+s, 2*j+t, 2*k+u, 2*l+v, 2*m+v, 2*n+u, 2*p+t, 2*q+s };
            val += spin_batch[ idx ];
          }
        }
      }
    }
    // Store significant elements only
    if ( abs(val) > 1e-14 ) {
      if ( store_sparse_spatial_array_ ) sparse_spatial_pdm[ it->first ] = factor * val;
      if ( store_full_spatial_array_ ) {
        if ( abs( spatial_fourpdm(i,j,k,l,m,n,p,q) ) > 1e-14 ) {
          cout << "repeated spatial indices!\n";
          assert(false);
        }
        spatial_fourpdm(i,j,k,l,m,n,p,q) = factor * val;
      }
    }
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Fourpdm_container::store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements)
{
  assert( new_spin_orbital_elements.size() == 70 );
  // Temporary batches of npdm elements
  std::map< std::vector<int>, double > spin_batch;
  std::map< std::vector<int>, double > spatial_batch;

  for (int idx=0; idx < new_spin_orbital_elements.size(); ++idx) {
    // Get all spin-index permutations
    Fourpdm_permutations p;
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
