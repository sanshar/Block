/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef TWOPDM_DRIVER_HEADER_H
#define TWOPDM_DRIVER_HEADER_H

#include <boost/format.hpp>
#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include "execinfo.h"

#include <multiarray.h>
#include <vector>
#include "npdm_driver.h"

namespace SpinAdapted{

//===========================================================================================================================================================

class Twopdm_driver : public Npdm_driver {

  public:
    Twopdm_driver( int sites );
    ~Twopdm_driver() {};
  
    array_4d<double> twopdm;
    array_4d<double> spatial_twopdm;

  private:
    void save_npdm_text(const int &i, const int &j);
    void save_npdm_binary(const int &i, const int &j);
    void build_spatial_npdm(const int &i, const int &j);
    void save_spatial_npdm_text(const int &i, const int &j);
    void save_spatial_npdm_binary(const int &i, const int &j);
    void load_npdm_binary(const int &i, const int &j);
  
    void resize_npdm_array(int dim) { twopdm.resize(dim,dim,dim,dim); spatial_twopdm.resize(dim/2,dim/2,dim/2,dim/2); }
    void clear_npdm_array() { twopdm.Clear(); spatial_twopdm.Clear(); }
    void accumulate_npdm();
  
  //  void save_averaged_twopdm(const int &nroots);
    void assign_npdm_antisymmetric(const int i, const int j, const int k, const int l, const double val);
    void store_npdm_elements(std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements);
    void calcenergy(int state);
  

};

//===========================================================================================================================================================

}

#endif
