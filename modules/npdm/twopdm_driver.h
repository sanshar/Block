/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef TWOPDM_DRIVER_HEADER_H
#define TWOPDM_DRIVER_HEADER_H

#include "npdm_driver.h"

namespace SpinAdapted{

//===========================================================================================================================================================

class Twopdm_driver : public Npdm_driver {

  public:
    Twopdm_driver( int sites );
    ~Twopdm_driver() {};
  
    // These maps are designed to hold elements computed at one sweep position only, but could still be memory-intensive.
    std::map< std::tuple<int,int,int,int>, double > spin_map;
    std::map< std::tuple<int,int,int,int>, double > spatial_map;
    // Optional arrays to store the full spin and/or spatial 2PDMs in core if memory allows.
    array_4d<double> twopdm;
    array_4d<double> spatial_twopdm;

    void save_npdms(const int &i, const int &j);

  private:
    void save_npdm_text(const int &i, const int &j);
    void save_npdm_binary(const int &i, const int &j);
    void save_spatial_npdm_text(const int &i, const int &j);
    void save_spatial_npdm_binary(const int &i, const int &j);
    void load_npdm_binary(const int &i, const int &j);
    void accumulate_npdm();
  
    void update_full_spin_array();
    void update_full_spatial_array();
    void build_spatial_elements();
    void screen_sparse_arrays();
    std::map< std::tuple<int,int,int,int>, int > get_spin_permutations( std::vector<int>& indices );
    void store_npdm_elements(std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements);
    void clear_sparse_arrays() { spin_map.clear(); spatial_map.clear(); }

};

//===========================================================================================================================================================

}

#endif
