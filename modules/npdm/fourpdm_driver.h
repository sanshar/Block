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
    Fourpdm_driver( int sites );
    ~Fourpdm_driver() {};
  
    array_8d<double> fourpdm;
    array_8d<double> spatial_fourpdm;

  private:
    void save_npdm_text(const int &i, const int &j);
    void save_npdm_binary(const int &i, const int &j);
    void build_spatial_npdm(const int &i, const int &j);
    void save_spatial_npdm_text(const int &i, const int &j);
    void save_spatial_npdm_binary(const int &i, const int &j);
    void load_npdm_binary(const int &i, const int &j);
  
    void resize_npdm_array(int dim) { fourpdm.resize(dim,dim,dim,dim,dim,dim,dim,dim); 
                                      spatial_fourpdm.resize(dim/2,dim/2,dim/2,dim/2,dim/2,dim/2,dim/2,dim/2); }
    void clear_npdm_array() { fourpdm.Clear(); spatial_fourpdm.Clear(); }
    void accumulate_npdm();
  
    void assign_npdm_antisymmetric(const int i, const int j, const int k, const int l, 
                                   const int m, const int n, const int p, const int q, const double val);
    void store_npdm_elements(std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements);
    void get_even_and_odd_perms( const std::vector<int> mnpq, std::vector< std::vector<int> > & even_perms, std::vector< std::vector<int> >& odd_perms );
    void calcenergy(int state);
  
};

//===========================================================================================================================================================

}

#endif
