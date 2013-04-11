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
  dotOps_(lhsOps),
  rhsOps_(lhsOps),
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
  test = { '(','(',')','(','C','C',')','(','D','D',')',')','(',')' };
  if ( build_pattern == test ) return CC_DD; // ????

  // 0,3,1
  test = { '(','(',')','(','C','C',')','(','D',')',')','(','D',')' };
  if ( build_pattern == test ) return CC_D_D;


  for (auto i = build_pattern.begin(); i != build_pattern.end(); i++) {
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
  indices.reserve( lhsOps_.indices.size() + dotOps_.indices.size() + rhsOps_.indices.size() );
  indices.insert( indices.end(), lhsOps_.indices.begin(), lhsOps_.indices.end() );
  indices.insert( indices.end(), dotOps_.indices.begin(), dotOps_.indices.end() );
  indices.insert( indices.end(), rhsOps_.indices.begin(), rhsOps_.indices.end() );
  assert (indices.size() == 4);

  // Set-up how tensor operator is constructed from (compound) block operators
  std::vector<char> build_pattern = { '(' };
  build_pattern.reserve( lhsOps_.build_pattern.size() + dotOps_.build_pattern.size() + rhsOps_.build_pattern.size() + 2 );
  build_pattern.insert( build_pattern.end(), lhsOps_.build_pattern.begin(), lhsOps_.build_pattern.end() );
  build_pattern.insert( build_pattern.end(), dotOps_.build_pattern.begin(), dotOps_.build_pattern.end() );
  build_pattern.push_back( ')' );
  build_pattern.insert( build_pattern.end(), rhsOps_.build_pattern.begin(), rhsOps_.build_pattern.end() );

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

  // Pointers to the numerical operator representations
  boost::shared_ptr<SparseMatrix> lhsPtr = lhsOps_.opReps.at(ilhs);
  boost::shared_ptr<SparseMatrix> dotPtr = dotOps_.opReps.at(idot);
  boost::shared_ptr<SparseMatrix> rhsPtr = rhsOps_.opReps.at(irhs);

  // Pointers to the transposes (may not all be needed)
  Transposeview lhsOpTr = Transposeview(lhsPtr);
  Transposeview dotOpTr = Transposeview(dotPtr);
  Transposeview rhsOpTr = Transposeview(rhsPtr);

  // Set actual pointers we'll contract
  boost::shared_ptr<SparseMatrix> lhsOp, dotOp, rhsOp;
  // LHS
  if ( lhsOps_.transpose )
    lhsOp = boost::shared_ptr<SparseMatrix>( &lhsOpTr, boostutils::null_deleter() );
  else
    lhsOp = lhsPtr;
  // Dot
  if ( dotOps_.transpose )
    dotOp = boost::shared_ptr<SparseMatrix>( &dotOpTr, boostutils::null_deleter() );
  else
    dotOp = dotPtr;
  // RHS
  if ( rhsOps_.transpose )
    rhsOp = boost::shared_ptr<SparseMatrix>( &rhsOpTr, boostutils::null_deleter() );
  else
    rhsOp = rhsPtr;

  // Push back into expectations the singlet component of this spin-operator contraction
  int index_begin = expectations_.size();
  spinExpectation(wavefunction_, wavefunction_, *lhsOp, *dotOp, *rhsOp, big_, expectations_, false);

  // Modify new elements with sign factors
  double factor = lhsOps_.factor * dotOps_.factor * rhsOps_.factor;
  for (int i = index_begin; i < expectations_.size(); i++) {
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
  for (int s = lhs2S+dot2S; s <= std::abs( lhs2S-dot2S ); s-=2) {
    if ( s == rhs2S ) return true;
  }
  return false;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_expectations::build_singlet_expectations()
{
  expectations_.clear();

  for (int ilhs = 0; ilhs < lhsOps_.opReps.size(); ilhs++) {
    for (int idot = 0; idot < dotOps_.opReps.size(); idot++) {
      for (int irhs = 0; irhs < rhsOps_.opReps.size(); irhs++) {
  
        // .mults should be redundant!  Test we're doing what we think we're doing
        assert( lhsOps_.mults.at(ilhs)-1 == lhsOps_.opReps.at(ilhs)->get_deltaQuantum().totalSpin );
        assert( dotOps_.mults.at(idot)-1 == dotOps_.opReps.at(idot)->get_deltaQuantum().totalSpin );
        assert( rhsOps_.mults.at(irhs)-1 == rhsOps_.opReps.at(irhs)->get_deltaQuantum().totalSpin );

        // If this combination allows a singlet, compute it
        if ( test_for_singlet( lhsOps_.mults.at(ilhs), dotOps_.mults.at(idot), rhsOps_.mults.at(irhs) ) ) {
          contract_spin_operators( ilhs, idot, irhs );
        }

      }
    }
  }
}

//===========================================================================================================================================================

}

