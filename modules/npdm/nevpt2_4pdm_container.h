/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NEVPT2_4PDM_CONTAINER_H
#define NEVPT2_4PDM_CONTAINER_H

#include "npdm_container.h"

namespace SpinAdapted{

//===========================================================================================================================================================

class Nevpt2_4pdm_container : public Npdm_container {

  public:
    Nevpt2_4pdm_container( array_6d<double>& a16_matrix ) : a16_matrix_(a16_matrix) {};
    ~Nevpt2_4pdm_container() {};
  
    void save_npdms(const int &i, const int &j) { };
    void clear_sparse_arrays() { };
    void store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements );

  private:
    array_6d<double>& a16_matrix_;

    void store_a16_contribution( std::map< std::vector<int>, double >& spatial_batch );
    void build_spatial_elements( std::map< std::vector<int>, double >& spin_batch, std::map< std::vector<int>, double >& spatial_batch );

};

//===========================================================================================================================================================

}

#endif
