/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include <boost/lexical_cast.hpp>
#include "MatrixBLAS.h"
#include "pario.h"
//FIXME use forward declaration for spinExpectation
#include "npdm_expectations_engine.h"
#include "npdm_expectations.h"
#include "npdm_patterns.h"

namespace SpinAdapted{
namespace Npdm{

// Forward declaration
void npdm_set_up_linear_equations(std::string& s, std::vector<double>& b0, Matrix& A, ColumnVector& b, std::vector< std::vector<int> >& so_indices);

//===========================================================================================================================================================

Npdm_expectations::Npdm_expectations( Npdm_patterns& npdm_patterns, const int order, Wavefunction & wavefunction, const SpinBlock & big )
: npdm_patterns_(npdm_patterns),
  npdm_order_(order),
  wavefunction_(wavefunction),
  big_(big)
{ }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// The pattern generator for the non-redundant NPDM elements leads to duplicates if indices are repeated, so some can be skipped explicitly.
// If the normal-ordered string is not of non-redundant form, then the original string produces duplicates when permutations are applied.

bool Npdm_expectations::screen_op_string_for_duplicates( std::string& op )
{
  std::vector<int> indices;  
  std::string CD;
  for (auto it = op.begin(); it != op.end(); ++it) {
    if ( (*it == '(') || (*it == ')') ) { continue; }
    else if ( (*it == 'C') || (*it == 'D') ) { CD.push_back(*it); }
    else { indices.push_back(*it); }
  }

  if ( indices.size() == 2 ) {
    // 1PDM case
//    return npdm_patterns_.screen_1pdm_strings( indices, CD ); 
    return false;
  }
  else if ( indices.size() == 4 ) {
    // 2PDM case
    return npdm_patterns_.screen_2pdm_strings( indices, CD ); 
  }
  else if ( indices.size() == 6 ) {
    // 3PDM case
    return npdm_patterns_.screen_3pdm_strings( indices, CD ); 
  }
  else if ( indices.size() == 8 ) {
    // 4PDM case
    return npdm_patterns_.screen_4pdm_strings( indices, CD ); 
  }
  else assert(false);

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

std::string Npdm_expectations::get_full_op_string( NpdmSpinOps_base & lhsOps, NpdmSpinOps_base & rhsOps, NpdmSpinOps_base & dotOps )
{

  // Set up npdm element indices
  std::vector<int> indices;
  indices.reserve( lhsOps.indices_.size() + dotOps.indices_.size() + rhsOps.indices_.size() );
  indices.insert( indices.end(), lhsOps.indices_.begin(), lhsOps.indices_.end() );
  indices.insert( indices.end(), dotOps.indices_.begin(), dotOps.indices_.end() );
  indices.insert( indices.end(), rhsOps.indices_.begin(), rhsOps.indices_.end() );
  //cout << "lhs indices = "; for (auto it = lhsOps.indices_.begin(); it != lhsOps.indices_.end(); ++it) { cout << *it << " "; } cout << std::endl;
  //cout << "dot indices = "; for (auto it = dotOps.indices_.begin(); it != dotOps.indices_.end(); ++it) { cout << *it << " "; } cout << std::endl;
  //cout << "rhs indices = "; for (auto it = rhsOps.indices_.begin(); it != rhsOps.indices_.end(); ++it) { cout << *it << " "; } cout << std::endl;
  //cout << "spatial indices = "; for (auto it = indices.begin(); it != indices.end(); ++it) { cout << *it << " "; } cout << std::endl;

  // Set up how tensor operator is constructed from (compound) block operators
  std::string build_pattern = "(";
  build_pattern.reserve( lhsOps.build_pattern_.size() + dotOps.build_pattern_.size() + rhsOps.build_pattern_.size() + 2 );
  build_pattern.insert( build_pattern.end(), lhsOps.build_pattern_.begin(), lhsOps.build_pattern_.end() );
  build_pattern.insert( build_pattern.end(), dotOps.build_pattern_.begin(), dotOps.build_pattern_.end() );
  build_pattern.push_back( ')' );
  build_pattern.insert( build_pattern.end(), rhsOps.build_pattern_.begin(), rhsOps.build_pattern_.end() );

  // Combine indices and build_pattern into one string
  std::string op_string;
  for (auto it = build_pattern.begin(); it != build_pattern.end(); ++it) {
    op_string.push_back(*it);
    if ( (*it == 'C') || (*it == 'D') ) {
      string index = boost::lexical_cast<string>( indices.at(0) );
      op_string.append(index);
      indices.erase( indices.begin() );  
    }
  }
  cout << op_string << std::endl;

  return op_string;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

double Npdm_expectations::contract_spin_adapted_operators( int ilhs, int idot, int irhs, 
                                                           NpdmSpinOps_base& lhsOps, NpdmSpinOps_base& rhsOps, NpdmSpinOps_base& dotOps )
{
//FIXME is null_deleter() necessary below?
  SparseMatrix* null = 0; 
  double expectation;

  boost::shared_ptr<SparseMatrix> lhsOp, dotOp, rhsOp;
  if ( lhsOps.opReps_.size() > 0 ) lhsOp = lhsOps.opReps_.at(ilhs);
  if ( dotOps.opReps_.size() > 0 ) dotOp = dotOps.opReps_.at(idot);
  if ( rhsOps.opReps_.size() > 0 ) rhsOp = rhsOps.opReps_.at(irhs);

  // We need to distinguish cases where one or more blocks has an empty operator string
  // X_X_X
  if ( (lhsOps.opReps_.size() > 0) && (dotOps.opReps_.size() > 0) && (rhsOps.opReps_.size() > 0) ) {
    Transposeview lhsOpTr = Transposeview(*lhsOp);
    if ( lhsOps.transpose_ ) lhsOp = boost::shared_ptr<SparseMatrix>( &lhsOpTr, boostutils::null_deleter() );
    Transposeview dotOpTr = Transposeview(*dotOp);
    if ( dotOps.transpose_ ) dotOp = boost::shared_ptr<SparseMatrix>( &dotOpTr, boostutils::null_deleter() );
    Transposeview rhsOpTr = Transposeview(*rhsOp);
    if ( rhsOps.transpose_ ) rhsOp = boost::shared_ptr<SparseMatrix>( &rhsOpTr, boostutils::null_deleter() );
    expectation = spinExpectation(wavefunction_, wavefunction_, *lhsOp, *dotOp, *rhsOp, big_);
  }
  // X_X_0
  else if ( (lhsOps.opReps_.size() > 0) && (dotOps.opReps_.size() > 0) && (rhsOps.opReps_.size() == 0) ) {
    Transposeview lhsOpTr = Transposeview(*lhsOp);
    if ( lhsOps.transpose_ ) lhsOp = boost::shared_ptr<SparseMatrix>( &lhsOpTr, boostutils::null_deleter() );
    Transposeview dotOpTr = Transposeview(*dotOp);
    if ( dotOps.transpose_ ) dotOp = boost::shared_ptr<SparseMatrix>( &dotOpTr, boostutils::null_deleter() );
    expectation = spinExpectation(wavefunction_, wavefunction_, *lhsOp, *dotOp, *null, big_);
  }
  // X_0_X
  else if ( (lhsOps.opReps_.size() > 0) && (dotOps.opReps_.size() == 0) && (rhsOps.opReps_.size() > 0) ) {
    Transposeview lhsOpTr = Transposeview(*lhsOp);
    if ( lhsOps.transpose_ ) lhsOp = boost::shared_ptr<SparseMatrix>( &lhsOpTr, boostutils::null_deleter() );
    Transposeview rhsOpTr = Transposeview(*rhsOp);
    if ( rhsOps.transpose_ ) rhsOp = boost::shared_ptr<SparseMatrix>( &rhsOpTr, boostutils::null_deleter() );
    expectation = spinExpectation(wavefunction_, wavefunction_, *lhsOp, *null, *rhsOp, big_);
  }
  // 0_X_X
  else if ( (lhsOps.opReps_.size() == 0) && (dotOps.opReps_.size() > 0) && (rhsOps.opReps_.size() > 0) ) {
    Transposeview dotOpTr = Transposeview(*dotOp);
    if ( dotOps.transpose_ ) dotOp = boost::shared_ptr<SparseMatrix>( &dotOpTr, boostutils::null_deleter() );
    Transposeview rhsOpTr = Transposeview(*rhsOp);
    if ( rhsOps.transpose_ ) rhsOp = boost::shared_ptr<SparseMatrix>( &rhsOpTr, boostutils::null_deleter() );
    expectation = spinExpectation(wavefunction_, wavefunction_, *null, *dotOp, *rhsOp, big_);
  }
  // X_0_0
  else if ( (lhsOps.opReps_.size() > 0) && (dotOps.opReps_.size() == 0) && (rhsOps.opReps_.size() == 0) ) {
    Transposeview lhsOpTr = Transposeview(*lhsOp);
    if ( lhsOps.transpose_ ) lhsOp = boost::shared_ptr<SparseMatrix>( &lhsOpTr, boostutils::null_deleter() );
    expectation = spinExpectation(wavefunction_, wavefunction_, *lhsOp, *null, *null, big_);
  }
  // 0_X_0
  else if ( (lhsOps.opReps_.size() == 0) && (dotOps.opReps_.size() > 0) && (rhsOps.opReps_.size() == 0) ) {
    Transposeview dotOpTr = Transposeview(*dotOp);
    if ( dotOps.transpose_ ) dotOp = boost::shared_ptr<SparseMatrix>( &dotOpTr, boostutils::null_deleter() );
    expectation = spinExpectation(wavefunction_, wavefunction_, *null, *dotOp, *null, big_);
  }
  // 0_0_X
  else if ( (lhsOps.opReps_.size() == 0) && (dotOps.opReps_.size() == 0) && (rhsOps.opReps_.size() > 0) ) {
    Transposeview rhsOpTr = Transposeview(*rhsOp);
    if ( rhsOps.transpose_ ) rhsOp = boost::shared_ptr<SparseMatrix>( &rhsOpTr, boostutils::null_deleter() );
    expectation = spinExpectation(wavefunction_, wavefunction_, *null, *null, *rhsOp, big_);
  }
  else assert(false);

  // Modify new element with sign factors and return
  double factor = lhsOps.factor_ * dotOps.factor_ * rhsOps.factor_;
//pout << "expectation, factor = " << expectation << ", " << factor  << std::endl;
  return expectation*factor;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_expectations::test_for_singlet( int lhs_mult, int dot_mult, int rhs_mult )
{
  // Work with 2*S instead of multiplicities
  int lhs2S = lhs_mult -1;
  int dot2S = dot_mult -1;
  int rhs2S = rhs_mult -1;

  // Couple LHS and Dot spin angular momenta and see if any equal RHS  
  for (int twoS = std::abs(lhs2S - dot2S); twoS <= ( lhs2S + dot2S ); twoS += 2 ) {
    if ( twoS == rhs2S ) {
//cout << "singlet!\n";
      return true;
    }
  }
  return false;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_expectations::build_spin_adapted_singlet_expectations( NpdmSpinOps_base & lhsOps, NpdmSpinOps_base & rhsOps, NpdmSpinOps_base & dotOps )
{
  expectations_.clear();

  // IMPORTANT: generate spin-components in the same order as RHS of linear equation solver in npdm_set_up_linear_equations routine
  // i.e. in accordance with the operator string build_pattern
  for (int irhs = 0; irhs < rhsOps.mults_.size(); ++irhs) {
//if (rhsOps.opReps_.size() > 0) { pout << "rhs: " << rhsOps.mults_.at(irhs) << std::endl; }
    for (int idot = 0; idot < dotOps.mults_.size(); ++idot) {
//if (dotOps.opReps_.size() > 0) { pout << "dot: " << dotOps.mults_.at(idot) << std::endl; }
      for (int ilhs = 0; ilhs < lhsOps.mults_.size(); ++ilhs) {
//if (lhsOps.opReps_.size() > 0) { pout << "lhs: " << lhsOps.mults_.at(ilhs) << std::endl; }
//pout << "spin comp: ilhs, idot, irhs = " << ilhs << idot << irhs << std::endl;

//FIXME get rid of mults_; just use get_deltaQuantum directly
        // Check that the spin multiplicities of the actual operators we've got are what we think they are!
        if ( lhsOps.opReps_.size() > 0 ) assert( lhsOps.mults_.at(ilhs) -1 == lhsOps.opReps_.at(ilhs)->get_deltaQuantum().totalSpin );
        if ( dotOps.opReps_.size() > 0 ) assert( dotOps.mults_.at(idot) -1 == dotOps.opReps_.at(idot)->get_deltaQuantum().totalSpin );
        if ( rhsOps.opReps_.size() > 0 ) assert( rhsOps.mults_.at(irhs) -1 == rhsOps.opReps_.at(irhs)->get_deltaQuantum().totalSpin );

        // Screen operator combinations that do not combine to give a singlet
        bool singlet = test_for_singlet( lhsOps.mults_.at(ilhs), dotOps.mults_.at(idot), rhsOps.mults_.at(irhs) );
        if ( singlet ) expectations_.push_back( contract_spin_adapted_operators( ilhs, idot, irhs, lhsOps, rhsOps, dotOps ) );
      }
    }
  }

  assert (expectations_.size() > 0);
//cout << "---------------------------------\n";
//cout << "spin-adapted expectations =\n";
////cout << "mpirank = " << mpigetrank() << endl;
//for (auto it = expectations_.begin(); it != expectations_.end(); ++it) {
//  cout << *it << std::endl;
//}
//cout << "---------------------------------\n";

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

std::vector< std::pair< std::vector<int>, double > > 
Npdm_expectations::get_nonspin_adapted_expectations( NpdmSpinOps_base & lhsOps, NpdmSpinOps_base & rhsOps, NpdmSpinOps_base & dotOps )
{
  // Initialize dimension of spin-adapted to non-spin-adapted transformation
  int dim;
  if ( npdm_order_ == 1 ) dim = 2;
  else if ( npdm_order_ == 2 ) dim = 6;
  else if ( npdm_order_ == 3 ) dim = 20;
  else if ( npdm_order_ == 4 ) dim = 70;
  else assert(false);
  std::vector< std::pair< std::vector<int>, double > > new_pdm_elements;

  // Get operator build string. e.g. (C2C4)(D5D6)
  std::string op_string = get_full_op_string( lhsOps, rhsOps, dotOps );

  // Screen away unwanted strings (e.g. those that produce duplicate NPDM elements)
  if ( screen_op_string_for_duplicates(op_string) ) return new_pdm_elements;

  // Contract spin-adapted spatial operators and build singlet expectation values
  build_spin_adapted_singlet_expectations( lhsOps, rhsOps, dotOps );
  
  // Now transform to non-spin-adapted spin-orbital representation
  // b holds the spin-adapted expectation values (we only care about the singlets)
  // Note the spin-order of elements in b follows a convention
  ColumnVector x(dim), b(dim);
  // x holds the non-spin-adapted expectation values
  x=0.0;

  // Transformation matrix
  Matrix A(dim,dim);
  // Vector of spin-orbital indices ordered according to A
  std::vector< std::vector<int> > so_indices(dim);

  // Parse operator string and set up linear equations
  npdm_set_up_linear_equations(op_string, expectations_, A, b, so_indices );

  // Solve A.x = b to get non-spin-adapted expectations in x
  xsolve_AxeqB(A, b, x);

  // Package transformed elements into container and return
  for (int i=0; i < so_indices.size(); ++i) {
    new_pdm_elements.push_back( std::make_pair(so_indices[i], x(i+1)) );
  } 

//pout << "x vector:\n";
//for (int i=1; i<(dim+1); ++i) { 
//    pout << i << "\t\t" << x(i) << std::endl;
//}
  
  return new_pdm_elements;
}

//===========================================================================================================================================================

}
}
