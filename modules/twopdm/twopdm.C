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

void assign_twopdm_elements( std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements, array_4d<double> & twopdm )
{
  for (int i=0; i < new_spin_orbital_elements.size(); ++i) {
    int ix = new_spin_orbital_elements[i].first[0];
    int jx = new_spin_orbital_elements[i].first[1];
    int kx = new_spin_orbital_elements[i].first[2];
    int lx = new_spin_orbital_elements[i].first[3];
    double x = new_spin_orbital_elements[i].second;
    assign_antisymmetric(twopdm, ix, jx, kx, lx, x);
  }

//FIXME is the transpose always needed?
  for (int i=0; i < new_spin_orbital_elements.size(); ++i) {
    int ix = new_spin_orbital_elements[i].first[3];
    int jx = new_spin_orbital_elements[i].first[2];
    int kx = new_spin_orbital_elements[i].first[1];
    int lx = new_spin_orbital_elements[i].first[0];
    double x = new_spin_orbital_elements[i].second;
    assign_antisymmetric(twopdm, ix, jx, kx, lx, x);
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void loop_over_block_operators( Wavefunction & wavefunction, 
                                const SpinBlock & big, 
                                std::vector<Npdm::CD> & lhs_cd_type,
                                std::vector<Npdm::CD> & dot_cd_type,
                                std::vector<Npdm::CD> & rhs_cd_type,
                                array_4d<double> & twopdm )
{
  SpinBlock* rhsBlock = big.get_rightBlock();
  SpinBlock* lhsdotBlock = big.get_leftBlock();
  SpinBlock* lhsBlock = lhsdotBlock->get_leftBlock();
  SpinBlock* dotBlock = lhsdotBlock->get_rightBlock();

  boost::shared_ptr<NpdmSpinOps> lhsOps = select_op_wrapper( lhsBlock, lhs_cd_type );
  boost::shared_ptr<NpdmSpinOps> dotOps = select_op_wrapper( dotBlock, dot_cd_type );
  boost::shared_ptr<NpdmSpinOps> rhsOps = select_op_wrapper( rhsBlock, rhs_cd_type );

  Npdm::Npdm_expectations npdm_expectations( wavefunction, big, *lhsOps, *dotOps, *rhsOps );
pout << "lhsOps->size()" << lhsOps->size() << std::endl;
pout << "dotOps->size()" << dotOps->size() << std::endl;
pout << "rhsOps->size()" << rhsOps->size() << std::endl;
pout << "------\n";

  // Only one spatial combination on the dot block
  assert( dotOps->size() == 1 );
  bool skip = dotOps->set_local_ops( 0 );
  if (skip) return;
  if ( lhsOps->opReps_.size() > 0 ) assert( dotOps->mults_.size() == dotOps->opReps_.size() );

  // Many spatial combinations on left block
  for ( int ilhs = 0; ilhs < lhsOps->size(); ++ilhs ) {
    skip = lhsOps->set_local_ops( ilhs );
    if (skip) continue;
pout << "lhsOps->indices_.size() " << lhsOps->indices_.size() << std::endl;
pout << "lhsOps->opReps_.size() " << lhsOps->opReps_.size() << std::endl;
pout << "dotOps->indices_.size() " << dotOps->indices_.size() << std::endl;
pout << "dotOps->opReps_.size() " << dotOps->opReps_.size() << std::endl;
    if ( lhsOps->opReps_.size() > 0 ) assert( lhsOps->mults_.size() == lhsOps->opReps_.size() );

    // Many spatial combinations on right block
    for ( int irhs = 0; irhs < rhsOps->size(); ++irhs ) {
      skip = rhsOps->set_local_ops( irhs );
      if (skip) continue;
pout << "rhsOps->indices_.size() " << rhsOps->indices_.size() << std::endl;
pout << "rhsOps->opReps_.size() " << rhsOps->opReps_.size() << std::endl;
pout << "-------------------------------------------------------------------------------------------\n";
pout << "spatial: ilhs, irhs = " << ilhs << ", " << irhs << std::endl;
      if ( rhsOps->opReps_.size() > 0 ) assert( rhsOps->mults_.size() == rhsOps->opReps_.size() );

      // Get non-spin-adapated 2PDM elements after building spin-adapted elements
      // FIXME magic number 6
      std::vector< std::pair< std::vector<int>, double > > new_spin_orbital_elements = npdm_expectations.get_nonspin_adapted_expectations( 6 );

      // Assign twopdm elements
      assign_twopdm_elements( new_spin_orbital_elements, twopdm );
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
  
  // Loop over NPDM operator patterns (here we initialize for 2PDM)
  Npdm::Npdm_patterns npdm_patterns( 2, sweepPos, endPos );

  for (auto pattern = npdm_patterns.ldr_cd_begin(); pattern != npdm_patterns.ldr_cd_end(); ++pattern) {
    std::vector<Npdm::CD> lhs_cd_type = pattern->at('l');
    std::vector<Npdm::CD> dot_cd_type = pattern->at('d');
    std::vector<Npdm::CD> rhs_cd_type = pattern->at('r');
//FIXME
std::vector<Npdm::CD> foo = {Npdm::DESTRUCTION,Npdm::CREATION,Npdm::DESTRUCTION};
if ( (lhs_cd_type == foo) || (rhs_cd_type == foo) || (dot_cd_type == foo) ) continue;
foo = {Npdm::DESTRUCTION,Npdm::DESTRUCTION,Npdm::CREATION};
if ( (lhs_cd_type == foo) || (rhs_cd_type == foo) || (dot_cd_type == foo) ) continue;
//foo = {Npdm::CREATION,Npdm::DESTRUCTION,Npdm::CREATION};
//if ( (lhs_cd_type == foo) || (rhs_cd_type == foo) || (dot_cd_type == foo) ) continue;
    // Compute all irreducible PDM elements generated by this block operator pattern at this sweep position
//if ( lhs_cd_type.size() == 2  &&
//     dot_cd_type.size() == 1  &&
//     rhs_cd_type.size() == 1 )
    loop_over_block_operators( wavefunctions.at(0), big, lhs_cd_type, dot_cd_type, rhs_cd_type, twopdm );
  }
  
  // Combine NPDM elements from this sweep point with others
  accumulate_twopdm(twopdm);
  save_twopdm_binary(twopdm, state, state);

}

//===========================================================================================================================================================

}

