/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef FOURPDM_DRIVER_HEADER_H
#define FOURPDM_DRIVER_HEADER_H

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

class Fourpdm_driver : public Npdm_driver {

public:
  Fourpdm_driver() : Npdm_driver(4) { }
  void save_npdm_text(const int &i, const int &j);
  void save_npdm_binary(const int &i, const int &j);
  void save_spatial_npdm_text(const int &i, const int &j);
  void save_spatial_npdm_binary(const int &i, const int &j);
  void load_npdm_binary(const int &i, const int &j);
  void npdm_resize_array(int dim) { fourpdm.resize(dim,dim,dim,dim,dim,dim,dim,dim); }
  void npdm_clear_array() { fourpdm.Clear(); }

protected:
  void accumulate_npdm();
  void assign_npdm_antisymmetric(const int i, const int j, const int k, const int l, 
                                 const int m, const int n, const int p, const int q, const double val);
  void assign_npdm_elements(std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements);
  void get_even_and_odd_perms( const std::vector<int> mnpq, std::vector< std::vector<int> > & even_perms, std::vector< std::vector<int> > & odd_perms );
  void calcenergy(int state);

private:
  array_8d<double> fourpdm;

};

//===========================================================================================================================================================

}

#endif
