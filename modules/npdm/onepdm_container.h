/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef ONEPDM_CONTAINER_H
#define ONEPDM_CONTAINER_H

#include "multiarray.h"
#include "npdm_container.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

class Onepdm_container : public Npdm_container {

  public:
    Onepdm_container( int sites );
    ~Onepdm_container() {};
  
    void save_npdms(const int &i, const int &j);
    void store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements );
    void clear() { onepdm.Clear(); spatial_onepdm.Clear(); nonredundant_elements.clear(); }

    array_2d<double>& get_spatial_onepdm() { return spatial_onepdm; }

  private:
    std::vector< std::pair< std::vector<int>, double > > nonredundant_elements;
    array_2d<double> onepdm;
    array_2d<double> spatial_onepdm;

    bool store_full_spin_array_ = true;
    bool store_nonredundant_spin_elements_ = false;

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
