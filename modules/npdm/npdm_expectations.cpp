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
#include "npdm_operators.h"
#include "npdm_patterns.h"

namespace SpinAdapted{
namespace Npdm{

// Forward declaration
void npdm_set_up_linear_equations(std::string& s, std::vector<double>& b0, Matrix& A, ColumnVector& b, std::vector< std::vector<int> >& so_indices);

//===========================================================================================================================================================

Npdm_expectations::Npdm_expectations( const int order, Wavefunction & wavefunction, const SpinBlock & big )
: npdm_order_(order),
  wavefunction_(wavefunction),
  big_(big)
{ }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

std::string Npdm_expectations::get_op_string( NpdmSpinOps_base & lhsOps, NpdmSpinOps_base & rhsOps, NpdmSpinOps_base & dotOps )
{

  // Set up npdm element indices
  std::vector<int> indices;
  indices.reserve( lhsOps.indices_.size() + dotOps.indices_.size() + rhsOps.indices_.size() );
  indices.insert( indices.end(), lhsOps.indices_.begin(), lhsOps.indices_.end() );
  indices.insert( indices.end(), dotOps.indices_.begin(), dotOps.indices_.end() );
  indices.insert( indices.end(), rhsOps.indices_.begin(), rhsOps.indices_.end() );
  assert( (indices.size() == 4) || (indices.size() == 6) || (indices.size() == 8) );
//  cout << "dot indices = "; for (auto it = dotOps.indices_.begin(); it != dotOps.indices_.end(); ++it) { cout << *it << " "; } cout << std::endl;
//  cout << "spatial indices = "; for (auto it = indices.begin(); it != indices.end(); ++it) { cout << *it << " "; } cout << std::endl;

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

//if ( lhsOps.opReps_.size() > 0 ) pout << "lhsOp:\n" << *lhsOp;
//if ( dotOps.opReps_.size() > 0 ) pout << "dotOp:\n" << *dotOp;
//if ( rhsOps.opReps_.size() > 0 ) pout << "rhsOp:\n" << *rhsOp;

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
//pout << "2S(lhs,dot,rhs) =   " << lhs2S << " " << dot2S << " " << rhs2S << std::endl;

//if ( (lhs2S == 0) && (dot2S == 0) && (rhs2S == 0) ) return true;
//else return false;

  // Couple LHS and Dot spin angular momenta and see if any equal RHS  
  for (int twoS = std::abs(lhs2S - dot2S); twoS <= ( lhs2S + dot2S ); twoS += 2 ) {
    if ( twoS == rhs2S ) {
//      pout << "\nsinglet found!\n";
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
    for (int idot = 0; idot < dotOps.mults_.size(); ++idot) {
      for (int ilhs = 0; ilhs < lhsOps.mults_.size(); ++ilhs) {
//pout << "spin comp: ilhs, idot, irhs = " << ilhs << idot << irhs << std::endl;

        // Check that the spin multiplicities of the actual operators we've got are what we think they are!
//if ( lhsOps.opReps_.size() > 0 ) {
//pout << lhsOps.mults_.at(ilhs) -1 << "          " << lhsOps.opReps_.at(ilhs)->get_deltaQuantum().totalSpin << std::endl;
//}
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
//cout << "mpirank = " << mpigetrank() << endl;
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
  if ( npdm_order_ == 2 ) dim = 6;
  else if ( npdm_order_ == 3 ) dim = 20;
  else if ( npdm_order_ == 4 ) dim = 70;
  else assert(false);

  // Contract spin-adapted spatial operators and build singlet expectation values
  build_spin_adapted_singlet_expectations( lhsOps, rhsOps, dotOps );
  
  // Now transform to non-spin-adapted spin-orbital representation
  // b holds the spin-adapted expectation values (we only care about the singlets)
  ColumnVector x(dim), b(dim);
  // x holds the non-spin-adapted expectation values
  x=0.0;

  // Transformation matrix
  Matrix A(dim,dim);
  // Vector of spin-orbital indices ordered according to A
  std::vector< std::vector<int> > so_indices(dim);

  // Parse operator string and set up linear equations
  std::string op_string = get_op_string( lhsOps, rhsOps, dotOps );
  npdm_set_up_linear_equations(op_string, expectations_, A, b, so_indices );

//pout << "A matrix:\n";
//for (int i=1; i<(dim+1); ++i) { 
//  pout << i << "\t\t";
//  for (int j=1; j<(dim+1); ++j) {
//    pout << "  " << A(i,j);
//  }
//  pout << std::endl;
//}

//pout << "b vector:\n";
//for (int i=1; i<(dim+1); ++i) { 
//    pout << i << "\t\t" << b(i) << std::endl;
//}

  // Solve A.x = b to get non-spin-adapted expectations in x
  xsolve_AxeqB(A, b, x);

  // Package transformed elements into container and return
  std::vector< std::pair< std::vector<int>, double > > new_pdm_elements;
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
