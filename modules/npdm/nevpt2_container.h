/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NEVPT2_CONTAINER_H
#define NEVPT2_CONTAINER_H


#include "multiarray.h"
#include "npdm_container.h"
#include "npdm_permutations.h"
#include <vector>

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

class Nevpt2_container : public Npdm_container {

  public:
    explicit Nevpt2_container( int sites );
    ~Nevpt2_container() {};
  
    void save_npdms(const int &i, const int &j);
    void store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements );
    void clear() { a16_matrix_.Clear(); a22_matrix_.Clear();}

  private:
    array_6d<double> a16_matrix_;
    array_6d<double> a22_matrix_;

    void save_full_array_text( array_6d<double>& matrix , char file[]);
    void save_full_array_text( array_4d<double>& matrix , char file[]);
    void save_full_array_text( array_2d<double>& matrix , char file[]);
    void save_full_array_text( std::vector<double>& matrix , char file[]);
    void update_1pdm_contribution( std::vector< std::pair< std::vector<int>, double > >& spin_batch );
    void update_2pdm_contribution( std::vector< std::pair< std::vector<int>, double > >& spin_batch );
    void update_3pdm_contribution( std::vector< std::pair< std::vector<int>, double > >& spin_batch );
    void update_4pdm_contribution( std::vector< std::pair< std::vector<int>, double > >& spin_batch );
    void accumulate_full_array(array_6d<double>& matrix);
    void accumulate_full_array(array_4d<double>& matrix);
    void accumulate_full_array(array_2d<double>& matrix);
    void accumulate_full_array(std::vector<double>& matrix);

};

//===========================================================================================================================================================

}
}

#endif
