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

//FIXME parallel    
boost::mpi::communicator world;
world.barrier();
    Timer timer;
    save_npdm_text(i, j);
    save_npdm_binary(i, j);
    save_spatial_npdm_text(i, j);
    save_spatial_npdm_binary(i, j);
world.barrier();
    pout << "1PDM save full array time " << timer.elapsedwalltime() << " " << timer.elapsedcputime() << endl;
//  }
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

void Onepdm_container::load_npdm_binary(const int &i, const int &j)
{
assert(false); // <<<< CAN WE RETHINK USE OF DISK FOR NPDM?
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/onepdm.", i, j,".bin");
    std::ifstream ifs(file, std::ios::binary);
    boost::archive::binary_iarchive load(ifs);
    load >> onepdm;
    ifs.close();
  }
//#ifndef SERIAL
//  mpi::communicator world;
//  mpi::broadcast(world,onepdm,0);
//FIXME this is a contradiction --- BUT IS IT GOOD TO CLEAR TO SAVE MEMORY???
//  if( mpigetrank() != 0)
//    onepdm.Clear();
//#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Onepdm_container::accumulate_npdm()
{
  assert(false);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Onepdm_container::update_full_spin_array( std::map< std::vector<int>, double >& spin_batch ) 
{
  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    int i = (it->first)[0];
    int j = (it->first)[1];

    double val = it->second;
//    if ( abs(val) > 1e-8 ) pout << "so-onepdm val: i,j = " << i << "," << j << "\t\t" << val << endl;
//    pout << "so-onepdm val: i,j = " << i << "," << j << "\t\t" << val << endl;


    // Test for duplicates
    if ( onepdm(i, j) != 0.0 && abs(onepdm(i,j)-val) > 1e-6) {
      void *array[10];
      size_t size;
      size = backtrace(array, 10);
      cout << "WARNING: Already calculated "<<i<<" "<<j<<endl;
      //backtrace_symbols_fd(array, size, 2);
      cout << "earlier value: "<<onepdm(i,j)<<endl<< "new value:     "<<val<<endl;
      assert( false );
    }
    if ( abs(val) > 1e-14 ) onepdm(i,j) = val;
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Onepdm_container::build_spatial_elements( std::map< std::vector<int>, double >& spin_batch, std::map< std::vector<int>, double >& spatial_batch )
{
  double factor = 1.0;

  for (auto it = spatial_batch.begin(); it != spatial_batch.end(); ++it) {
    int i = (it->first)[0];
    int j = (it->first)[1];
    // Sum over spin indices 
    double val = 0.0;
    for (int s=0; s<2; s++) {
      std::vector<int> idx = { 2*i+s, 2*j+s };
      val += spin_batch[ idx ];
    }
    // Store significant elements only
    if ( abs(val) > 1e-14 ) {
      if ( abs( spatial_onepdm(i,j) ) > 1e-14 ) { 
        cout << "repeated spatial indices!\n";
        assert(false);
      }
      spatial_onepdm(i,j) = factor * val;
    }
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Onepdm_container::store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements)
{
  assert( new_spin_orbital_elements.size() == 2 );
  // Temporary batches of npdm elements
  std::map< std::vector<int>, double > spin_batch;
  std::map< std::vector<int>, double > spatial_batch;

  for (int idx=0; idx < new_spin_orbital_elements.size(); ++idx) {
    // Get all spin-index permutations
    Onepdm_permutations p;
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
    }
  }

  // Build and store new spatial elements
  build_spatial_elements( spin_batch, spatial_batch );
  update_full_spin_array( spin_batch );

}

//===========================================================================================================================================================

}

