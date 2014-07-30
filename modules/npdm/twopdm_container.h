/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef TWOPDM_CONTAINER_H
#define TWOPDM_CONTAINER_H

#include "npdm_container.h"
#include "npdm_symmetric_array.hpp"
#include "npdm_symmetric_spatial_array.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

class Twopdm_container : public Npdm_container {

  public:
    Twopdm_container( int sites );
    //FIXME destructor? -> depends on how this is used, so it should be defined.
    ~Twopdm_container() { }

    void save_npdms(const int &i, const int &j);
    void store_npdm_elements(const std::vector<std::pair<std::vector<int>, double> > & in);
//  void clear() { twopdm.Clear(); spatial_twopdm.Clear(); nonredundant_elements.clear(); }
    void clear() { twopdm.fill(0.0); spatial_twopdm.fill(0.0); }

    symmetric_spatial_array<double, 2>& get_spatial_twopdm() { assert(store_full_spatial_array_); return spatial_twopdm; }

    void update_array_component();

  private:

    bool store_full_spin_array_;

    bool store_full_spatial_array_;

    /// Temporary storage to store non-redundant elements 
    std::vector<std::pair<std::vector<int>, double> > tmp_store_; 

    /// Spin-orbital 2PDM, with full permutation symmetry
    /// only stored on master process
    symmetric_array<double, 2> twopdm;

    /// Spatial 2PDM, with full permutation symmetry
    /// this is optional, to be computed if store_full_spatial_array_ = true
    symmetric_spatial_array<double, 2> spatial_twopdm;

    void build_full_spatial_array();

    void save_npdm_text(const int &i, const int &j);
    void save_npdm_binary(const int &i, const int &j);

    void save_spatial_npdm_text(const int &i, const int &j);
    void save_spatial_npdm_binary(const int &i, const int &j);

    void load_npdm_binary(const int &i, const int &j);
  
};

//===========================================================================================================================================================

}
}

#endif
