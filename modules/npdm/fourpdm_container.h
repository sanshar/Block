/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef FOURPDM_CONTAINER_H
#define FOURPDM_CONTAINER_H

#include "npdm_container.h"
#include "npdm.h"

namespace SpinAdapted{

//===========================================================================================================================================================

class Fourpdm_container : public Npdm_container {

  public:
    Fourpdm_container( int sites );
    ~Fourpdm_container() {};
  
    // These maps are designed to hold elements computed at one sweep position only, but could still be memory-intensive.
    std::map< std::vector<int>, double > sparse_spin_pdm;
    std::map< std::vector<int>, double > sparse_spatial_pdm;
    // Optional arrays to store the full spin and/or spatial PDMs in core if memory allows.
    array_8d<double> fourpdm;
    array_8d<double> spatial_fourpdm;

    void save_npdms(const int &i, const int &j);
    void update_full_spin_array();
    void update_full_spatial_array();
    void clear_sparse_arrays() { sparse_spin_pdm.clear(); sparse_spatial_pdm.clear(); }
    void store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements );

  private:
    void save_npdm_text(const int &i, const int &j);
    void save_npdm_binary(const int &i, const int &j);
    void save_spatial_npdm_text(const int &i, const int &j);
    void save_spatial_npdm_binary(const int &i, const int &j);
    void load_npdm_binary(const int &i, const int &j);
    void accumulate_npdm();
  
    void get_even_and_odd_perms( const std::vector<int> mnpq, 
                                 std::vector< std::vector<int> > & even_perms, 
                                 std::vector< std::vector<int> > & odd_perms );
    std::map< std::vector<int>, int > get_spin_permutations( const std::vector<int>& indices );

    void build_spatial_elements( std::map< std::vector<int>, double >& spin_batch, 
                                 std::map< std::vector<int>, double >& spatial_batch );

};

//===========================================================================================================================================================

}

#endif
