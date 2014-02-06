/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "twopdm_driver.h"

namespace SpinAdapted{

//===========================================================================================================================================================

Twopdm_driver::Twopdm_driver( int sites ) : Npdm_driver(2) 
{
  store_full_spin_array_ = true;
  store_full_spatial_array_ = true;

  if ( store_full_spin_array_ ) {
    twopdm.resize(2*sites,2*sites,2*sites,2*sites);
    twopdm.Clear();
  }
  if ( store_full_spatial_array_ ) {
    spatial_twopdm.resize(sites,sites,sites,sites);
    spatial_twopdm.Clear();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
  
void Twopdm_driver::save_npdms(const int& i, const int& j)
{
//  save_sparse_array(i,j);

//  if ( store_full_array_ ) {

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
    pout << "2PDM save full array time " << timer.elapsedwalltime() << " " << timer.elapsedcputime() << endl;
//  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_driver::save_npdm_text(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/twopdm.", i, j,".txt");
    ofstream ofs(file);
    ofs << twopdm.dim1() << endl;
    double trace = 0.0;
    for(int k=0;k<twopdm.dim1();++k)
      for(int l=0;l<twopdm.dim2();++l)
        for(int m=0;m<twopdm.dim3();++m)
          for(int n=0;n<twopdm.dim4();++n) {
            ofs << boost::format("%d %d %d %d %20.14e\n") % k % l % m % n % twopdm(k,l,m,n);
            if ( (k==n) && (l==m) ) trace += twopdm(k,l,m,n);
          }
    ofs.close();
    std::cout << "Spin-orbital 2PDM trace = " << trace << "\n";
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_driver::save_spatial_npdm_text(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_twopdm.", i, j,".txt");
    ofstream ofs(file);
    ofs << spatial_twopdm.dim1() << endl;
    double trace = 0.0;
    for(int k=0;k<spatial_twopdm.dim1();++k)
      for(int l=0;l<spatial_twopdm.dim2();++l)
        for(int m=0;m<spatial_twopdm.dim3();++m)
          for(int n=0;n<spatial_twopdm.dim4();++n) {
            ofs << boost::format("%d %d %d %d %20.14e\n") % k % l % m % n % spatial_twopdm(k,l,m,n);
            if ( (k==n) && (l==m) ) trace += spatial_twopdm(k,l,m,n);
          }
    ofs.close();
    std::cout << "Spatial      2PDM trace = " << trace << "\n";
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_driver::save_npdm_binary(const int &i, const int &j)
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

void Twopdm_driver::save_spatial_npdm_binary(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_twopdm.", i, j,".bin");
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << spatial_twopdm;
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_driver::load_npdm_binary(const int &i, const int &j)
{
assert(false); // <<<< CAN WE RETHINK USE OF DISK FOR NPDM?
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/twopdm.", i, j,".bin");
    std::ifstream ifs(file, std::ios::binary);
    boost::archive::binary_iarchive load(ifs);
    load >> twopdm;
    ifs.close();
  }
//#ifndef SERIAL
//  mpi::communicator world;
//  mpi::broadcast(world,twopdm,0);
//FIXME this is a contradiction --- BUT IS IT GOOD TO CLEAR TO SAVE MEMORY???
//  if( mpigetrank() != 0)
//    twopdm.Clear();
//#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_driver::accumulate_npdm()
{
#ifndef SERIAL
  array_4d<double> tmp_recv;
  mpi::communicator world;
  if( mpigetrank() == 0) {
    for(int i=1;i<world.size();++i) {
      world.recv(i, i, tmp_recv);
      for(int k=0;k<twopdm.dim1();++k)
        for(int l=0;l<twopdm.dim2();++l)
          for(int m=0;m<twopdm.dim3();++m)
            for(int n=0;n<twopdm.dim4();++n)
              if(tmp_recv(k,l,m,n) != 0.) twopdm(k,l,m,n) = tmp_recv(k,l,m,n);
	 }
  }
  else {
    world.send(0, mpigetrank(), twopdm);
  }
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_driver::update_full_spin_array()
{
  for (auto it = spin_map.begin(); it != spin_map.end(); ++it) {
    int i = std::get<0>( it->first );
    int j = std::get<1>( it->first );
    int k = std::get<2>( it->first );
    int l = std::get<3>( it->first );

    double val = it->second;
    //if ( abs(val) > 1e-8 ) pout << "so-twopdm val: i,j,k,l = " << i << "," << j << "," << k << "," << l << "\t\t" << val << endl;
    //pout << "so-twopdm val: i,j,k,l = " << i << "," << j << "," << k << "," << l << "\t\t" << val << endl;

    // Test for duplicates
    if ( twopdm(i, j, k, l) != 0.0 && abs(twopdm(i,j,k,l)-val) > 1e-6) {
      void *array[10];
      size_t size;
      size = backtrace(array, 10);
      cout << "WARNING: Already calculated "<<i<<" "<<j<<" "<<k<<" "<<l<<endl;
      //backtrace_symbols_fd(array, size, 2);
      cout << "earlier value: "<<twopdm(i,j,k,l)<<endl<< "new value:     "<<val<<endl;
      assert( false );
    }
    twopdm(i,j,k,l) = val;
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_driver::update_full_spatial_array()
{
  for (auto it = spatial_map.begin(); it != spatial_map.end(); ++it) {
    int i = std::get<0>( it->first );
    int j = std::get<1>( it->first );
    int k = std::get<2>( it->first );
    int l = std::get<3>( it->first );
    spatial_twopdm(i,j,k,l) = it->second;
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_driver::build_spatial_elements()
{
  //Note we multiply the spatial 2PDM by a factor of 1/2 to be consistent with the old BLOCK code, but this doesn't seem conventional?
//FIXME parallelization!
//  if( mpigetrank() == 0)
//  {
    double factor = 0.5;

    for (auto it = spatial_map.begin(); it != spatial_map.end(); ++it) {
      int i = std::get<0>( it->first );
      int j = std::get<1>( it->first );
      int k = std::get<2>( it->first );
      int l = std::get<3>( it->first );

      double val = 0.0;
      for (int s=0; s<2; s++) {
        for (int t =0; t<2; t++) {
          val += spin_map[ std::make_tuple( 2*i+s, 2*j+t, 2*k+t, 2*l+s ) ];
        }
      }
      it->second = factor * val;
    }
// }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

std::map< std::tuple<int,int,int,int>, int > Twopdm_driver::get_spin_permutations( std::vector<int>& indices )
{

  std::map< std::tuple<int,int,int,int>, int > perms;
  int i = indices[0];
  int j = indices[1];
  int k = indices[2];
  int l = indices[3];

  // 8 permutations
  //--------------------------
  perms[std::make_tuple(i, j, k, l)] = 1;
  perms[std::make_tuple(i, j, l, k)] = -1;
  perms[std::make_tuple(j, i, k, l)] = -1;
  perms[std::make_tuple(j, i, l, k)] = 1;

  perms[std::make_tuple(l, k, j, i)] = 1;
  perms[std::make_tuple(k, l, j, i)] = -1;
  perms[std::make_tuple(l, k, i, j)] = -1;
  perms[std::make_tuple(k, l, i, j)] = 1;

  return perms;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_driver::screen_sparse_arrays()
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

void Twopdm_driver::store_npdm_elements( std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements)
{
  assert( new_spin_orbital_elements.size() == 6 );

  for (int idx=0; idx < new_spin_orbital_elements.size(); ++idx) {
    // Get all spin-index permutations
    std::map< std::tuple<int,int,int,int>, int > spin_indices = get_spin_permutations( new_spin_orbital_elements[idx].first );
    // Assign batch of spin-orbital 2PDM values
    double val = new_spin_orbital_elements[idx].second;
    for (auto it = spin_indices.begin(); it != spin_indices.end(); ++it) {
      spin_map[ it->first ] = it->second * val;
    }
    // Initialize spatial 2PDM indices
    for (auto it = spin_indices.begin(); it != spin_indices.end(); ++it) {
      int i = std::get<0>( it->first );
      int j = std::get<1>( it->first );
      int k = std::get<2>( it->first );
      int l = std::get<3>( it->first );
      spatial_map[ std::make_tuple(i/2, j/2, k/2, l/2) ] = 0.0;
    }
  }

  // Build and store batch of spatial 2PDM elements
  build_spatial_elements();

  // Screen away small or zero elements
  screen_sparse_arrays();

}

//===========================================================================================================================================================

}

