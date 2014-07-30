/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef ONEPDM_CONTAINER_H
#define ONEPDM_CONTAINER_H

#include "npdm_container.h"
#include "npdm_symmetric_array.hpp"
#include "npdm_symmetric_spatial_array.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

class Onepdm_container : public Npdm_container {

  public:
    Onepdm_container( int sites );
    //FIXME destructor? -> depends on how this is used, so it should be defined.
    ~Onepdm_container() { }

    void save_npdms(const int &i, const int &j);
    void store_npdm_elements(const std::vector<std::pair<std::vector<int>, double> > & in);
    void clear() { onepdm.fill(0.0); spatial_onepdm.fill(0.0); }

    symmetric_spatial_array<double, 1>& get_spatial_onepdm() { assert(store_full_spatial_array_); return spatial_onepdm; }

    void update_array_component();

  private:

    bool store_full_spin_array_;

    bool store_full_spatial_array_;

    /// Temporary storage to store non-redundant elements 
    std::vector<std::pair<std::vector<int>, double> > tmp_store_; 

    /// Spin-orbital 1PDM, with full permutation symmetry
    /// only stored on master process
    symmetric_array<double, 1> onepdm;

    /// Spatial 1PDM, with full permutation symmetry
    /// this is optional, to be computed if store_full_spatial_array_ = true
    symmetric_spatial_array<double, 1> spatial_onepdm;

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
