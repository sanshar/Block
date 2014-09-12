/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef PAIRPDM_CONTAINER_H
#define PAIRPDM_CONTAINER_H

#include "multiarray.h"
#include "npdm_container.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

class Pairpdm_container : public Npdm_container {

  public:
    Pairpdm_container( int sites );
    ~Pairpdm_container() {};
  
    void save_npdms(const int &i, const int &j);
    void store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements );
    void clear() { pairpdm.Clear(); spatial_pairpdm.Clear(); nonredundant_elements.clear(); }

    array_2d<double>& get_spatial_pairpdm() { return spatial_pairpdm; }

  private:
    std::vector< std::pair< std::vector<int>, double > > nonredundant_elements;
    array_2d<double> pairpdm;
    array_2d<double> spatial_pairpdm;

    void save_npdm_text(const int &i, const int &j);
    void save_npdm_binary(const int &i, const int &j);
    void save_spatial_npdm_text(const int &i, const int &j);
    void save_spatial_npdm_binary(const int &i, const int &j);
    void load_npdm_binary(const int &i, const int &j);
    void accumulate_npdm();
    void accumulate_spatial_npdm();
  
    void update_full_spin_array( std::vector< std::pair< std::vector<int>, double > >& spin_batch );
    void update_full_spatial_array( std::vector< std::pair< std::vector<int>, double > >& spin_batch );
    void calculate_spatial_npdm();

};

//===========================================================================================================================================================

}
}

#endif
