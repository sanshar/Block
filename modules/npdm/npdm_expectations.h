/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NPDM_EXPECT_H
#define NPDM_EXPECT_H

#include "npdm_patterns.h"
#include "npdm_expectations.h"
#include "npdm_operator_wrappers.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

class Npdm_expectations {
 
  public:
    Npdm_expectations( Wavefunction & wavefunction, 
                       const SpinBlock & big, 
                       NpdmSpinOps & lhsOps,
                       NpdmSpinOps & dotOps,
                       NpdmSpinOps & rhsOps );

    std::vector< std::pair< std::vector<int>, double > > get_nonspin_adapted_expectations( int dim );
//    void transform_spin_adapt_to_nonspin_adapt( array_4d<double> & twopdm );
//    void old_transform_spin_adapt_to_nonspin_adapt( array_4d<double> & twopdm );

  private:
    std::vector< double > expectations_;
    Wavefunction & wavefunction_; 
    const SpinBlock & big_; 
    NpdmSpinOps & lhsOps_;
    NpdmSpinOps & dotOps_;
    NpdmSpinOps & rhsOps_;

    void build_spin_adapted_singlet_expectations();
    double contract_spin_adapted_operators( int ilhs, int idot, int irhs );
    bool test_for_singlet( int lhs_mult, int dot_mult, int rhs_mult );
//    Oporder old_parse_build_pattern( std::vector<char> build_pattern );
    std::string get_op_string();

};

//===========================================================================================================================================================

}
}

#endif
