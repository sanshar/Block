/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NPDM_PERMUTATIONS_H
#define NPDM_PERMUTATIONS_H

#include <vector>
#include <multiarray.h>
#include "spinblock.h"
#include "wavefunction.h"
#include "BaseOperator.h"
#include <algorithm>
#include "npdm_epermute.h"

namespace SpinAdapted{

//===========================================================================================================================================================

class Npdm_permutations {
  public:
    Npdm_permutations() {}
    virtual ~Npdm_permutations() {}

    void process_new_elements( const std::vector< std::pair< std::vector<int>, double > >& in, 
                               std::vector< std::pair< std::vector<int>, double > >& nonredundant_elements,
                               std::vector< std::pair< std::vector<int>, double > >& spin_batch );
  private:
    virtual void get_spin_permutations( std::vector<std::pair<std::vector<int>,double> >& spin_batch, 
                                        const std::vector<int>& indices, const double& val ) = 0;
};

//===========================================================================================================================================================

class Onepdm_permutations : public Npdm_permutations {
  private:
    void get_spin_permutations( std::vector<std::pair<std::vector<int>,double> >& spin_batch, const std::vector<int>& indices, const double& val );
};

//===========================================================================================================================================================

class Twopdm_permutations : public Npdm_permutations {
  private:
    void get_spin_permutations( std::vector<std::pair<std::vector<int>,double> >& spin_batch, const std::vector<int>& indices, const double& val );
};

//===========================================================================================================================================================

class Threepdm_permutations : public Npdm_permutations {
  private:
    void get_spin_permutations( std::vector<std::pair<std::vector<int>,double> >& spin_batch, const std::vector<int>& indices, const double& val );
};

//===========================================================================================================================================================

class Fourpdm_permutations : public Npdm_permutations {
  private:
    void get_even_and_odd_perms( const std::vector<int> mnpq, std::vector<std::vector<int> > & even_perms, std::vector<std::vector<int> > & odd_perms );
    void get_spin_permutations( std::vector<std::pair<std::vector<int>,double> >& spin_batch, const std::vector<int>& indices, const double& val );
    void get_spin_permutations_general(std::vector<std::pair<std::vector<int>,double> >& spin_batch, const std::vector<int>& indices, const double& val);
};

//===========================================================================================================================================================

}

#endif
