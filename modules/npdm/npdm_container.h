/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NPDM_CONTAINER_H
#define NPDM_CONTAINER_H

#include <tuple>
#include <boost/mpi.hpp>
#include <vector>
#include <multiarray.h>
#include "spinblock.h"
#include "wavefunction.h"
#include "BaseOperator.h"

namespace SpinAdapted{

//===========================================================================================================================================================

class Npdm_container {

  public:
    Npdm_container() {}
    virtual ~Npdm_container() {}
  
    virtual void store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements ) = 0;
    virtual void clear_sparse_arrays() = 0;
    virtual void save_npdms(const int &i, const int &j) = 0;

};

//===========================================================================================================================================================

}

#endif
