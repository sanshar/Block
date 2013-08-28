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
//#include "npdm_patterns.h"
//#include "npdm_expectations.h"
//#include "npdm_operator_wrappers.h"

namespace SpinAdapted{

//===========================================================================================================================================================

class Twopdm_driver : public Npdm_driver {

public:
  Twopdm_driver() : Npdm_driver(2) { }
  void save_npdm_text(const int &i, const int &j);
  void save_npdm_binary(const int &i, const int &j);
  void save_spatial_npdm_text(const int &i, const int &j);
  void save_spatial_npdm_binary(const int &i, const int &j);
  void load_npdm_binary(const int &i, const int &j);
  void resize_array(int dim) { twopdm.resize(dim,dim,dim,dim); }
  void clear_array() { twopdm.Clear(); }

protected:
//  void save_averaged_twopdm(const int &nroots);
  void accumulate_npdm();
  void assign_npdm_antisymmetric(const int i, const int j, const int k, const int l, const double val);
  void assign_npdm_elements(std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements);
  void calcenergy(int state);

private:
  array_4d<double> twopdm;

};

//===========================================================================================================================================================

}

#endif
