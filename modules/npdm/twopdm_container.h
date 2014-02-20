/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef TWOPDM_CONTAINER_H
#define TWOPDM_CONTAINER_H

#include "npdm_container.h"

namespace SpinAdapted{

//===========================================================================================================================================================

class Twopdm_container : public Npdm_container {

  public:
    Twopdm_container( int sites );
    ~Twopdm_container() {};
  
    void clear_sparse_arrays() { sparse_spin_pdm.clear(); sparse_spatial_pdm.clear(); }
    void store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements );
    void save_npdms(const int &i, const int &j);

    array_4d<double>& get_spatial_twopdm() { assert(store_full_spatial_array_); return spatial_twopdm; }
    std::map< std::vector<int>, double >& get_sparse_spatial_pdm() { assert(store_sparse_spatial_array_); return sparse_spatial_pdm; }

  private:
    // These maps are designed to hold elements computed at one sweep position only, but could still be memory-intensive.
    std::map< std::vector<int>, double > sparse_spin_pdm;
    std::map< std::vector<int>, double > sparse_spatial_pdm;
    // Optional arrays to store the full spin and/or spatial PDMs in core if memory allows.
    array_4d<double> twopdm;
    array_4d<double> spatial_twopdm;

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
    void build_spatial_elements( std::map< std::vector<int>, double >& spin_batch, std::map< std::vector<int>, double >& spatial_batch );

};

//===========================================================================================================================================================

}

#endif
