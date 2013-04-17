/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "twopdm.h"
#include "npdm_expectations.h"
#include "npdm_operators.h"
#include "npdm_patterns.h"

namespace SpinAdapted{

//===========================================================================================================================================================

Npdm_spin_adaptation::Npdm_spin_adaptation( NpdmSpinOps & lhsOps,
                                            NpdmSpinOps & dotOps,
                                            NpdmSpinOps & rhsOps,
                                            array_4d<double> & twopdm )
: lhsOps_(lhsOps),
  dotOps_(dotOps),
  rhsOps_(rhsOps),
  twopdm_(twopdm)
{ }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

Oporder Npdm_spin_adaptation::parse_build_pattern( std::vector<char> build_pattern )
{

  std::vector<char> test;

  // 2,2,0
  test = { '(','(','C','C',')','(','D','D',')',')','(',')' };
  if ( build_pattern == test ) return CC_DD;
  test = { '(','(','C','D',')','(','C','D',')',')','(',')' };
  if ( build_pattern == test ) return CD_CD;

  // 2,1,1
  test = { '(','(','C','C',')','(','D',')',')','(','D',')' };
  if ( build_pattern == test ) return CC_D_D;
  test = { '(','(','C','D',')','(','C',')',')','(','D',')' };
  if ( build_pattern == test ) return CD_CD;  // ???
  test = { '(','(','C','D',')','(','D',')',')','(','C',')' };
  if ( build_pattern == test ) return CD_D_C; 

  // 1,3,0
  test = { '(','(','C',')','(','C','(','D','D',')',')',')','(',')' };
  if ( build_pattern == test ) return CC_DD;  // ???

  // 1,2,1
  test = { '(','(','C',')','(','C','D',')',')','(','D',')' };
  if ( build_pattern == test ) return C_CD_D;
  test = { '(','(','C',')','(','D','D',')',')','(','C',')' };
  if ( build_pattern == test ) return D_CC_D; // ????

  // 0,4,0
  test = { '(','(',')','(','(','C','C',')','(','D','D',')',')',')','(',')' };
  if ( build_pattern == test ) return CC_DD; // ????

  // 0,3,1
  test = { '(','(',')','(','C','C',')','(','D',')',')','(','D',')' };
  if ( build_pattern == test ) return CC_D_D;

  for (auto i = build_pattern.begin(); i != build_pattern.end(); ++i) {
    std::cout << *i;
  }
  std::cout << std::endl;
  assert( false );

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_spin_adaptation::to_nonspin_adapt( std::vector<double> & spin_adapted_vals )
{

  // Set-up npdm element indices
  std::vector<int> indices;
  indices.reserve( lhsOps_.indices_.size() + dotOps_.indices_.size() + rhsOps_.indices_.size() );
  indices.insert( indices.end(), lhsOps_.indices_.begin(), lhsOps_.indices_.end() );
  indices.insert( indices.end(), dotOps_.indices_.begin(), dotOps_.indices_.end() );
  indices.insert( indices.end(), rhsOps_.indices_.begin(), rhsOps_.indices_.end() );
  assert (indices.size() == 4);

  // Set-up how tensor operator is constructed from (compound) block operators
  std::vector<char> build_pattern = { '(' };
  build_pattern.reserve( lhsOps_.build_pattern_.size() + dotOps_.build_pattern_.size() + rhsOps_.build_pattern_.size() + 2 );
  build_pattern.insert( build_pattern.end(), lhsOps_.build_pattern_.begin(), lhsOps_.build_pattern_.end() );
  build_pattern.insert( build_pattern.end(), dotOps_.build_pattern_.begin(), dotOps_.build_pattern_.end() );
  build_pattern.push_back( ')' );
  build_pattern.insert( build_pattern.end(), rhsOps_.build_pattern_.begin(), rhsOps_.build_pattern_.end() );

  // Translate our format for the build pattern into the types used by the old twopdm implementation
  Oporder build_pattern_type = parse_build_pattern( build_pattern );  

  // Call old code for transforming spin-adapted expectatio values and twopdm update
  spin_to_nonspin( indices, spin_adapted_vals, twopdm_, build_pattern_type, true );

}

//===========================================================================================================================================================

Npdm_expectations::Npdm_expectations( Wavefunction & wavefunction, 
                                      const SpinBlock & big, 
                                      NpdmSpinOps & lhsOps,
                                      NpdmSpinOps & dotOps,
                                      NpdmSpinOps & rhsOps )
: wavefunction_(wavefunction),
  big_(big),
  lhsOps_(lhsOps),
  dotOps_(dotOps),
  rhsOps_(rhsOps)
{ }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_expectations::contract_spin_operators( int ilhs, int idot, int irhs )
{

//////  SpinBlock* leftBlock = big_.get_leftBlock();
//////  SpinBlock* rightBlock = big_.get_rightBlock();
//////  SpinBlock* dotBlock = leftBlock->get_rightBlock();
//////
//////  boost::shared_ptr<SparseMatrix> lhsOp = leftBlock->get_leftBlock()->get_op_array(CRE).get_local_element(0)[0];
//////  boost::shared_ptr<SparseMatrix> dotOp = dotBlock->get_op_array(CRE_DES).get_local_element(0)[0];
//////  boost::shared_ptr<SparseMatrix> rhsOp = rightBlock->get_op_array(CRE).get_local_element(0)[0];
//////
//////  Transposeview rhsOpTr = Transposeview(*rhsOp);
//////  spinExpectation(wavefunction_, wavefunction_, *lhsOp, *dotOp, rhsOpTr, big_, expectations_, false);


  boost::shared_ptr<SparseMatrix> lhsPtr, dotPtr, rhsPtr, lhsOp, dotOp, rhsOp;

  // Pointers to the numerical operator representations or their tranposes
  if ( lhsOps_.opReps_.size() > 0 ) {
pout << "hello1  " << ilhs << std::endl;
    assert( lhsOps_.mults_.at(ilhs) -1 == lhsOps_.opReps_.at(ilhs)->get_deltaQuantum().totalSpin );
    lhsOp = lhsOps_.opReps_.at(ilhs);
//    lhsPtr = lhsOps_.opReps_.at(ilhs);
//    Transposeview lhsOpTr = Transposeview(*lhsPtr);
//    if ( lhsOps_.transpose_ ) lhsOp = boost::shared_ptr<SparseMatrix>( &lhsOpTr, boostutils::null_deleter() );
//    else lhsOp = lhsPtr;
  }
  if ( dotOps_.opReps_.size() > 0 ) {
pout << "hello2  " << idot << std::endl;
    assert( dotOps_.mults_.at(idot) -1 == dotOps_.opReps_.at(idot)->get_deltaQuantum().totalSpin );
    dotOp = dotOps_.opReps_.at(idot);
//    dotPtr = dotOps_.opReps_.at(idot);
//    Transposeview dotOpTr = Transposeview(*dotPtr);
//    if ( dotOps_.transpose_ ) dotOp = boost::shared_ptr<SparseMatrix>( &dotOpTr, boostutils::null_deleter() );
//    else dotOp = dotPtr;
  }
  if ( rhsOps_.opReps_.size() > 0 ) {
pout << "hello3  " << irhs << std::endl;
    assert( rhsOps_.mults_.at(irhs) -1 == rhsOps_.opReps_.at(irhs)->get_deltaQuantum().totalSpin );
    rhsOp = rhsOps_.opReps_.at(irhs);
//    rhsPtr = rhsOps_.opReps_.at(irhs);
//    Transposeview rhsOpTr = Transposeview(*rhsPtr);
//    if ( rhsOps_.transpose_ ) rhsOp = boost::shared_ptr<SparseMatrix>( &rhsOpTr, boostutils::null_deleter() );
//    else rhsOp = rhsPtr;
  }

  // Push-back into expectations the singlet component of this spin-operator contraction
  int index_begin = expectations_.size();

  // We need to distinguish cases where one or more blocks has an empty operator string
  assert ( dotOps_.opReps_.size() != 0 );
  SparseMatrix* null = 0; 
  if ( (lhsOps_.opReps_.size() == 0) && (rhsOps_.opReps_.size() == 0) )
    // 0_X_0 case
    spinExpectation(wavefunction_, wavefunction_, *null, *dotOp, *null, big_, expectations_, false);
  else if ( (lhsOps_.opReps_.size() == 0) && (rhsOps_.opReps_.size() > 0) )
    // 0_X_X case
    spinExpectation(wavefunction_, wavefunction_, *null, *dotOp, *rhsOp, big_, expectations_, false);
  else if ( (lhsOps_.opReps_.size() > 0) && (rhsOps_.opReps_.size() == 0) )
    // X_X_0 case
    spinExpectation(wavefunction_, wavefunction_, *lhsOp, *dotOp, *null, big_, expectations_, false);
  else {
    // X_X_X case
pout << "hello4\n";
//    spinExpectation(wavefunction_, wavefunction_, *lhsOp, *dotOp, *rhsOp, big_, expectations_, false);
    Transposeview rhsOpTr = Transposeview(*rhsOp);
    spinExpectation(wavefunction_, wavefunction_, *lhsOp, *dotOp, rhsOpTr, big_, expectations_, false);
  }


  // Modify new elements with sign factors
  double factor = lhsOps_.factor_ * dotOps_.factor_ * rhsOps_.factor_;
  for (int i = index_begin; i < expectations_.size(); ++i) {
    expectations_[i] *= factor;
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_expectations::test_for_singlet( int lhs_mult, int dot_mult, int rhs_mult )
{
  // Work with 2*S
  int lhs2S = lhs_mult-1;
  int dot2S = dot_mult-1;
  int rhs2S = rhs_mult-1;

  // Couple LHS and Dot spin angular momenta and see if any equal RHS  
  for (int s = lhs2S + dot2S; s <= std::abs( lhs2S - dot2S ); s -= 2) {
    if ( s == rhs2S ) return true;
  }
  return false;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Filter operator combinations that do not combine to give a singlet
void Npdm_expectations::build_singlet_expectations()
{
  expectations_.clear();
  bool singlet;

  // We need to distinguish cases where one or more blocks has an empty operator string
  assert ( dotOps_.opReps_.size() > 0 );

  if ( (lhsOps_.opReps_.size() == 0) && (rhsOps_.opReps_.size() == 0) ) {
    // 0_X_0 case
    for (int idot = 0; idot < dotOps_.opReps_.size(); ++idot) {
      if ( dotOps_.mults_.at(idot) == 1) contract_spin_operators( 0, idot, 0 );
    }
  }
  else if ( (lhsOps_.opReps_.size() == 0) && (rhsOps_.opReps_.size() > 0) ) {
    // 0_X_X case
    for (int idot = 0; idot < dotOps_.opReps_.size(); ++idot) {
      for (int irhs = 0; irhs < rhsOps_.opReps_.size(); ++irhs) {
        singlet = test_for_singlet( 0, dotOps_.mults_.at(idot), rhsOps_.mults_.at(irhs) );
        if ( singlet ) contract_spin_operators( 0, idot, irhs );
      }
    }
  }
  else if ( (lhsOps_.opReps_.size() > 0) && (rhsOps_.opReps_.size() == 0) ){
    // X_X_0 case
    for (int ilhs = 0; ilhs < lhsOps_.opReps_.size(); ++ilhs) {
      for (int idot = 0; idot < dotOps_.opReps_.size(); ++idot) {
        singlet = test_for_singlet( lhsOps_.mults_.at(ilhs), dotOps_.mults_.at(idot), 0 );
        if ( singlet ) contract_spin_operators( ilhs, idot, 0 );
      }
    }
  }
  else {
    // X_X_X case
    for (int ilhs = 0; ilhs < lhsOps_.opReps_.size(); ++ilhs) {
      for (int idot = 0; idot < dotOps_.opReps_.size(); ++idot) {
        for (int irhs = 0; irhs < rhsOps_.opReps_.size(); ++irhs) {
pout << "spin comp: ilhs, idot, irhs = " << ilhs << idot << irhs << std::endl;
  
          bool singlet = test_for_singlet( lhsOps_.mults_.at(ilhs), dotOps_.mults_.at(idot), rhsOps_.mults_.at(irhs) );
          if ( singlet ) contract_spin_operators( ilhs, idot, irhs );
        }
      }
    }
  }

}

//===========================================================================================================================================================

}

