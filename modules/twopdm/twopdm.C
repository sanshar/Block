/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "twopdm.h"
#include "npdm_patterns.h"
#include "npdm_expectations.h"
#include "npdm_operator_wrappers.h"

namespace SpinAdapted{

// Forward declaration
boost::shared_ptr<NpdmSpinOps> select_op_wrapper( SpinBlock * spinBlock, std::vector<Npdm::CD> & cd_type );

//===========================================================================================================================================================

void loop_over_block_operators( Wavefunction & wavefunction, 
                                const SpinBlock & big, 
                                std::vector<Npdm::CD> & lhs_cd_type,
                                std::vector<Npdm::CD> & dot_cd_type,
                                std::vector<Npdm::CD> & rhs_cd_type,
                                array_4d<double> & twopdm )
{

pout << "-------------------------------\n";
pout << "CD pattern (0=Cre, 1=Des):\n";
for (auto it = lhs_cd_type.begin(); it != lhs_cd_type.end(); ++it) {
  pout << *it;
}
pout << ",";
for (auto it = dot_cd_type.begin(); it != dot_cd_type.end(); ++it) {
  pout << *it;
}
pout << ",";
for (auto it = rhs_cd_type.begin(); it != rhs_cd_type.end(); ++it) {
  pout << *it;
}
pout << std::endl;

  SpinBlock* rhsBlock = big.get_rightBlock();
  SpinBlock* lhsdotBlock = big.get_leftBlock();
  SpinBlock* lhsBlock = lhsdotBlock->get_leftBlock();
  SpinBlock* dotBlock = lhsdotBlock->get_rightBlock();
//  SpinBlock* lhsBlock = big.get_leftBlock();
//  SpinBlock* dotBlock = lhsBlock->get_rightBlock();

  boost::shared_ptr<NpdmSpinOps> lhsOps = select_op_wrapper( lhsBlock, lhs_cd_type );
  boost::shared_ptr<NpdmSpinOps> dotOps = select_op_wrapper( dotBlock, dot_cd_type );
  boost::shared_ptr<NpdmSpinOps> rhsOps = select_op_wrapper( rhsBlock, rhs_cd_type );

  Npdm_expectations npdm_expectations( wavefunction, big, *lhsOps, *dotOps, *rhsOps );
  Npdm_spin_adaptation npdm_spin_adapt( *lhsOps, *dotOps, *rhsOps, twopdm );

  // Only one spatial combination on the dot block
  assert( dotOps->size() == 1 );
  dotOps->set_local_ops( 0 );
  assert( dotOps->mults_.size() == dotOps->opReps_.size() );

pout << "lhsOps->size()" << lhsOps->size() << std::endl;
pout << "dotOps->size()" << dotOps->size() << std::endl;
pout << "rhsOps->size()" << rhsOps->size() << std::endl;
pout << "------\n";
  // Many spatial combinations on left block
  for ( int ilhs = 0; ilhs < lhsOps->size(); ++ilhs ) {
    lhsOps->set_local_ops( ilhs );
pout << "lhsOps->indices_.size()" << lhsOps->indices_.size() << std::endl;
pout << "lhsOps->opReps_.size()" << lhsOps->opReps_.size() << std::endl;
pout << "dotOps->indices_.size()" << dotOps->indices_.size() << std::endl;
pout << "dotOps->opReps_.size()" << dotOps->opReps_.size() << std::endl;
    if ( lhsOps->opReps_.size() > 0 ) assert( lhsOps->mults_.size() == lhsOps->opReps_.size() );

    // Many spatial combinations on right block
    for ( int irhs = 0; irhs < rhsOps->size(); ++irhs ) {
      rhsOps->set_local_ops( irhs );
pout << "----------------------------------------\n";
pout << "rhsOps->indices_.size()" << rhsOps->indices_.size() << std::endl;
pout << "rhsOps->opReps_.size()" << rhsOps->opReps_.size() << std::endl;
pout << "------\n";
pout << "spatial: ilhs, irhs = " << ilhs << "," << irhs << std::endl;
      if ( rhsOps->opReps_.size() > 0 ) assert( rhsOps->mults_.size() == rhsOps->opReps_.size() );

      // Contract these spin-adapted spatial operators and build expectation values
      std::vector<double> vals = npdm_expectations.get_expectations();
      // Transform from spin-adapated expectations to non-spin-adapted and store as twopdm elements
      assert( vals.size() > 0 );
      npdm_spin_adapt.to_nonspin_adapt( vals ); 

    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void compute_twopdm_sweep(std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int state, int sweepPos, int endPos)
{

pout << "Sweep position = "<< sweepPos << " (" << endPos << ")\n";

  // 2pdm array built so far
  array_4d<double> twopdm(big.size()*2, big.size()*2, big.size()*2, big.size()*2);
  load_twopdm_binary(twopdm, state, state);
  
  // Loop over NPDM operator patterns (here we initialize for a 2PDM)
  // FIXME change this so the patterns are altered if we're at the edges of the sweep!
  // Npdm_patterns npdm_patterns( 2, sweepPos, endPos );
  Npdm::Npdm_patterns npdm_patterns( 2 );

  for (auto pattern = npdm_patterns.ldr_cd_begin(); pattern != npdm_patterns.ldr_cd_end(); ++pattern) {
    std::vector<Npdm::CD> lhs_cd_type = pattern->at('l');
    std::vector<Npdm::CD> dot_cd_type = pattern->at('d');
    std::vector<Npdm::CD> rhs_cd_type = pattern->at('r');
    // Compute all irreducible pdm elements generated by this block operator pattern at this sweep position
//FIXME only do some patterns
if ( lhs_cd_type.size() == 1 &&
     dot_cd_type.size() == 3 &&
     rhs_cd_type.size() == 0 )
    loop_over_block_operators( wavefunctions.at(0), big, lhs_cd_type, dot_cd_type, rhs_cd_type, twopdm );
  }
  
  // Combine NPDM elements from this sweep point with others
  accumulate_twopdm(twopdm);
  save_twopdm_binary(twopdm, state, state);

}

//===========================================================================================================================================================

}

