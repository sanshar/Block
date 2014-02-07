/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef THREEPDM_CONTAINER_H
#define THREEPDM_CONTAINER_H

#include "npdm_container.h"

namespace SpinAdapted{

//===========================================================================================================================================================

class Threepdm_container : public Npdm_container {

  public:
    Threepdm_container( int sites );
    ~Threepdm_container() {};

    void save_npdms(const int &i, const int &j);
    void clear_sparse_arrays() { sparse_spin_pdm.clear(); sparse_spatial_pdm.clear(); }
    void store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements );

    array_6d<double>& get_spatial_threepdm() { assert(store_full_spatial_array_); return spatial_threepdm; }

  private:
    // These maps are designed to hold elements computed at one sweep position only, but could still be memory-intensive.
    std::map< std::vector<int>, double > sparse_spin_pdm;
    std::map< std::vector<int>, double > sparse_spatial_pdm;
    // Optional arrays to store the full spin and/or spatial PDMs in core if memory allows.
    array_6d<double> threepdm;
    array_6d<double> spatial_threepdm;

    bool store_full_spin_array_;
    bool store_full_spatial_array_;
    bool store_sparse_spin_array_;
    bool store_sparse_spatial_array_;

    void save_npdm_text(const int &i, const int &j);
    void save_npdm_binary(const int &i, const int &j);
    void save_spatial_npdm_text(const int &i, const int &j);
    void save_spatial_npdm_binary(const int &i, const int &j);
    void load_npdm_binary(const int &i, const int &j);
    void accumulate_npdm();
  
    void update_full_spin_array( std::map< std::vector<int>, double >& spin_batch );
    std::map< std::vector<int>, int > get_spin_permutations( const std::vector<int>& indices );

    void build_spatial_elements( std::map< std::vector<int>, double >& spin_batch, 
                                 std::map< std::vector<int>, double >& spatial_batch );

};

//===========================================================================================================================================================

}

#endif
