/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NEVPT2_A16_MATRIX_H
#define NEVPT2_A16_MATRIX_H

#include "multiarray.h"
#include "npdm_container.h"
#include "npdm_permutations.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

class Nevpt2_A16_matrix : public Npdm_container {

  public:
    explicit Nevpt2_A16_matrix( int sites );
    ~Nevpt2_A16_matrix() {};
  
    void save_npdms(const int &i, const int &j);
    void store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements );

  private:
    array_6d<double> a16_matrix_;

    void save_A16_matrix_text();

    void build_spatial_2pdm_elements( std::map< std::vector<int>, double >& spin_batch, std::map< std::vector<int>, double >& spatial_batch );
    void build_spatial_3pdm_elements( std::map< std::vector<int>, double >& spin_batch, std::map< std::vector<int>, double >& spatial_batch );
    void build_spatial_4pdm_elements( std::map< std::vector<int>, double >& spin_batch, std::map< std::vector<int>, double >& spatial_batch );

    void store_A16_2pdm_contribution( std::map< std::vector<int>, double >& spatial_batch );
    void store_A16_3pdm_contribution( std::map< std::vector<int>, double >& spatial_batch );
    void store_A16_4pdm_contribution( std::map< std::vector<int>, double >& spatial_batch );

    void build_npdm_batch( Npdm_permutations& p, const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements,
                           std::map< std::vector<int>, double >& spatial_batch );
};

//===========================================================================================================================================================

}
}

#endif
