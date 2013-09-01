/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/
#include <boost/format.hpp>
#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif

#include "spinblock.h"
#include "npdm_sparse_array.h"

namespace SpinAdapted{

//===========================================================================================================================================================

void Npdm_sparse_array::insert( std::vector< std::pair< std::vector<int>, double > > & new_data ) 
{
  sparse_array_.insert( sparse_array_.end(), new_data.begin(), new_data.end() );
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_sparse_array::dump_file(int i, int j)
{
  char file[5000];
  sprintf( file, "%s%s%d.%d.p%d", dmrginp.save_prefix().c_str(),"/nonredundant_npdm.", i, j, mpigetrank() );
  std::ofstream ofs(file, std::ios::binary);
  boost::archive::binary_oarchive save(ofs);
  save << sparse_array_;
  ofs.close();
}
  
//===========================================================================================================================================================

}

