/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NPDM_SPARSE_ARRAY_HEADER_H
#define NPDM_SPARSE_ARRAY_HEADER_H

#include <vector>

//FIXME NPDM namespace!!!
//FIXME destructors!!!!
namespace SpinAdapted{

//===========================================================================================================================================================
// This is a real simple class now, but can obviously be extended (e.g. finite-size buffer, auto dump to disk when full etc...)
//===========================================================================================================================================================

class Npdm_sparse_array {

  public:
    Npdm_sparse_array() { sparse_array_.clear(); }

    void insert( std::vector< std::pair< std::vector<int>, double > > & new_data );
    void dump_file( int i, int j );

  private:
    std::vector< std::pair< std::vector<int>, double > > sparse_array_;

};
  
//===========================================================================================================================================================

}

#endif

