/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef FOURPDM_CONTAINER_H
#define FOURPDM_CONTAINER_H

#include "multiarray.h"
#include "npdm_container.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

class Fourpdm_container : public Npdm_container {

  public:
    Fourpdm_container( int sites );
//FIXME destructor?
  
    void save_npdms(const int &i, const int &j);
    void store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements );

    array_8d<double>& get_spatial_fourpdm() { assert(store_full_spatial_array_); return spatial_fourpdm; }

  private:
    // Vector to store nonredundant spin-orbital elements only
    std::vector< std::pair< std::vector<int>, double > > nonredundant_elements;
    // Optional arrays to store the full spin and/or spatial PDMs in core if memory allows.
    array_8d<double> fourpdm;
    array_8d<double> spatial_fourpdm;

    bool store_full_spin_array_;
    bool store_full_spatial_array_;
    bool store_nonredundant_spin_elements_;

    void save_npdm_text(const int &i, const int &j);
    void save_npdm_binary(const int &i, const int &j);
    void save_spatial_npdm_text(const int &i, const int &j);
    void save_spatial_npdm_binary(const int &i, const int &j);
    void load_npdm_binary(const int &i, const int &j);
    void accumulate_npdm();
    void accumulate_spatial_npdm();
  
    void update_full_spin_array( std::vector< std::pair< std::vector<int>, double > >& spin_batch );
    void update_full_spatial_array( std::vector< std::pair< std::vector<int>, double > >& spin_batch );

};

//===========================================================================================================================================================

}
}

#endif
