/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "global.h"
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/map.hpp>
#include "npdm_array_buffer.h"

namespace SpinAdapted{
namespace Npdm{

const unsigned int MAX_BUFFER_SIZE = 10000000;

//===========================================================================================================================================================

double Npdm_array_buffer::operator()(int i, int j, int k, int l, int m, int n) const
{
  assert((0 <= i) && (i < dim_));
  assert((0 <= j) && (j < dim_));
  assert((0 <= k) && (k < dim_));
  assert((0 <= l) && (l < dim_));
  assert((0 <= m) && (m < dim_));
  assert((0 <= n) && (n < dim_));
  std::vector<int> vec = { i,j,k,l,m,n };
  return data_.at( vec );
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

double& Npdm_array_buffer::operator()(int i, int j, int k, int l, int m, int n)
{
  assert((0 <= i) && (i < dim_));
  assert((0 <= j) && (j < dim_));
  assert((0 <= k) && (k < dim_));
  assert((0 <= l) && (l < dim_));
  assert((0 <= m) && (m < dim_));
  assert((0 <= n) && (n < dim_));
  std::vector<int> vec = { i,j,k,l,m,n };

  // Note if element doesn't exit, it is now zero-initialized (according to C++11 standard)
  double& ref = data_[ vec ];
  // Check buffer size in case new element created
  if ( data_.size() > MAX_BUFFER_SIZE ) {
    close_buffer();
    bufferID_++; 
    data_.clear();
    return data_[ vec ];
  }
  return ref;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_array_buffer::close_buffer()
{
  char file[5000];
  sprintf( file, "%s%s%d.%d.p%d.n%d%s", dmrginp.save_prefix().c_str(),"/partial_A16_matrix.", 0, 0, mpigetrank(), bufferID_, ".bin" );
  std::ofstream ofs(file, std::ios::binary);
  boost::archive::binary_oarchive save(ofs);
  save << data_;
  ofs.close();
}
  
//===========================================================================================================================================================

}
}

