/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef THREEPDM_DRIVER_HEADER_H
#define THREEPDM_DRIVER_HEADER_H

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

class Threepdm_driver : public Npdm_driver {

  public:
    Threepdm_driver() : Npdm_driver(3) { }
  
  private:
    void save_npdm_text(const int &i, const int &j);
    void save_npdm_binary(const int &i, const int &j);
    void save_spatial_npdm_text(const int &i, const int &j);
    void save_spatial_npdm_binary(const int &i, const int &j);
    void load_npdm_binary(const int &i, const int &j);
    void resize_npdm_array(int dim) { threepdm.resize(dim,dim,dim,dim,dim,dim); }
    void clear_npdm_array() { threepdm.Clear(); }
    void accumulate_npdm();

    void assign_npdm_antisymmetric(const int i, const int j, const int k, const int l, const int m, const int n, const double val);
    void assign_npdm_elements(std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements);
  //  void calcenergy(int state);
  
//FIXME this gets too big for practical calcs.  Need to consider a sparse storage, or direct to disk?
    array_6d<double> threepdm;

};

//===========================================================================================================================================================

}

#endif
