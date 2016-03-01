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
#include "externalsort.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

class Fourpdm_container : public Npdm_container {

  public:
    Fourpdm_container( int sites );
//FIXME destructor?
  
    void save_npdms(const int &i, const int &j);
    void store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements );
    void clear() { fourpdm.Clear(); spatial_fourpdm.Clear(); nonredundant_elements.clear(); }

    array_8d<double>& get_spatial_fourpdm() { assert(dmrginp.store_spinpdm()); return spatial_fourpdm; }

  private:
    // Vector to store nonredundant spin-orbital elements only
    std::vector< std::pair< std::vector<int>, double > > nonredundant_elements;
    // Optional arrays to store the full spin and/or spatial PDMs in core if memory allows.
    array_8d<double> fourpdm;
    array_8d<double> spatial_fourpdm;
    std::vector<int> elements_stride_;
    std::vector<Sortpdm::batch_index> nonspin_batch;
    FILE* spatpdm_disk;
    FILE* batch_index_file;

    void save_npdm_text(const int &i, const int &j);
    void save_npdm_binary(const int &i, const int &j);
    void save_spatial_npdm_text(const int &i, const int &j);
    void save_spatial_npdm_binary(const int &i, const int &j);
    void load_npdm_binary(const int &i, const int &j);
    void accumulate_npdm();
    void accumulate_spatial_npdm();
    void external_sort_index(const int &i, const int &j);
  
    void update_full_spin_array( std::vector< std::pair< std::vector<int>, double > >& spin_batch );
    void update_full_spatial_array( std::vector< std::pair< std::vector<int>, double > >& spin_batch );
    long oneindex_spin(const std::vector<int> & orbital_element_index);
    void dump_to_disk(std::vector< std::pair< std::vector<int>, double > > & spin_batch);

    void dump_text_to_disk(std::vector< std::pair< std::vector<int>, double > > & spin_batch);
    void dump_binary_to_disk(std::vector< std::pair< std::vector<int>, double > > & spin_batch);
    void merge_diskfile(const int &i, const int &j);
};

//===========================================================================================================================================================

}
}

#endif
