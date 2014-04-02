/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NEVPT2_A16_MATRIX_H
#define NEVPT2_A16_MATRIX_H

//#define DEBUG_A16_FULL_MATRIX

#include "multiarray.h"
#include "npdm_container.h"
#include "npdm_permutations.h"
#include "npdm_array_buffer.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================
// This class computes the 3 <EEEE> contributions to the A-matrix as given in eqn.(A16) of JCP 117, 9138 (2002)

class Nevpt2_A16_matrix : public Npdm_container {

  public:
    explicit Nevpt2_A16_matrix( int sites );
    ~Nevpt2_A16_matrix() {};
  
    void save_npdms(const int &i, const int &j);
    void store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements );

  private:
#ifdef DEBUG_A16_FULL_MATRIX
    array_6d<double> a16_matrix_;
    void accumulate_full_array();
#else
    Npdm_array_buffer a16_matrix_;
    void accumulate_files();
#endif
    void save_full_array_text( array_6d<double>& matrix );
    void update_1pdm_contribution( std::vector< std::pair< std::vector<int>, double > >& spin_batch );
    void update_2pdm_contribution( std::vector< std::pair< std::vector<int>, double > >& spin_batch );
    void update_3pdm_contribution( std::vector< std::pair< std::vector<int>, double > >& spin_batch );
    void update_4pdm_contribution( std::vector< std::pair< std::vector<int>, double > >& spin_batch );

};

//===========================================================================================================================================================

}
}

#endif
