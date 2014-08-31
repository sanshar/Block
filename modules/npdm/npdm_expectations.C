/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include <algorithm>
#include <boost/lexical_cast.hpp>
#include "MatrixBLAS.h"
#include "pario.h"
#include "npdm_expectations_engine.h"
#include "npdm_expectations.h"
#include "pario.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

Npdm_expectations::Npdm_expectations( Npdm_spin_adaptation& spin_adaptation, Npdm_patterns& npdm_patterns, 
                                      const int order, Wavefunction & wavefunction0, Wavefunction & wavefunction1, const SpinBlock & big )
: spin_adaptation_(spin_adaptation),
  npdm_patterns_(npdm_patterns),
  npdm_order_(order),
  wavefunction_0(wavefunction0),
  wavefunction_1(wavefunction1),
  big_(big)
{ }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// The pattern generator for the non-redundant NPDM elements leads to duplicates if indices are repeated, so some can be skipped explicitly.
// If the normal-ordered string is not of non-redundant form, then the original string produces duplicates when permutations are applied.

bool Npdm_expectations::screen_op_string_for_duplicates( const std::string& op, const std::vector<int>& indices )
{
  std::string CD;
  for (auto it = op.begin(); it != op.end(); ++it) {
    if ( (*it == 'C') || (*it == 'D') ) CD.push_back(*it);
  }

  if ( indices.size() == 2 ) {
    // 1PDM case
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
  else abort();

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_expectations::get_full_op_string( NpdmSpinOps_base & lhsOps, NpdmSpinOps_base & rhsOps, NpdmSpinOps_base & dotOps, 
                                            std::string& op_string, std::vector<int>& indices )
{

  // Set up npdm element indices
  indices.clear();
  indices.insert( indices.end(), lhsOps.indices_.begin(), lhsOps.indices_.end() );
  indices.insert( indices.end(), dotOps.indices_.begin(), dotOps.indices_.end() );
  indices.insert( indices.end(), rhsOps.indices_.begin(), rhsOps.indices_.end() );
  //cout << "lhs indices = "; for (auto it = lhsOps.indices_.begin(); it != lhsOps.indices_.end(); ++it) { cout << *it << " "; } cout << std::endl;
  //cout << "dot indices = "; for (auto it = dotOps.indices_.begin(); it != dotOps.indices_.end(); ++it) { cout << *it << " "; } cout << std::endl;
  //cout << "rhs indices = "; for (auto it = rhsOps.indices_.begin(); it != rhsOps.indices_.end(); ++it) { cout << *it << " "; } cout << std::endl;
  //cout << "spatial indices = "; for (auto it = indices.begin(); it != indices.end(); ++it) { cout << *it << " "; } cout << std::endl;

  // Record how tensor operator is constructed from (compound) block operators
  std::string build_pattern = "(";
  build_pattern.reserve( lhsOps.build_pattern_.size() + dotOps.build_pattern_.size() + rhsOps.build_pattern_.size() + 2 );
  build_pattern.insert( build_pattern.end(), lhsOps.build_pattern_.begin(), lhsOps.build_pattern_.end() );
  build_pattern.insert( build_pattern.end(), dotOps.build_pattern_.begin(), dotOps.build_pattern_.end() );
  build_pattern.push_back( ')' );
  build_pattern.insert( build_pattern.end(), rhsOps.build_pattern_.begin(), rhsOps.build_pattern_.end() );

  // Combine indices and build_pattern into one string
  op_string.clear();
  std::vector<int> idx = indices;
  for (auto it = build_pattern.begin(); it != build_pattern.end(); ++it) {
    op_string.push_back(*it);
    if ( (*it == 'C') || (*it == 'D') ) {
      std::string index = boost::lexical_cast<string>( idx.at(0) );
      op_string.append(index);
      idx.erase( idx.begin() );  
    }
  }
//cout << op_string << std::endl;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

double Npdm_expectations::contract_spin_adapted_operators( int ilhs, int idot, int irhs, 
                                                           NpdmSpinOps_base& lhsOps, NpdmSpinOps_base& rhsOps, NpdmSpinOps_base& dotOps )
{
//FIXME is null_deleter() necessary below?
  SparseMatrix* null = 0; 
  double expectation = 0;

  boost::shared_ptr<SparseMatrix> lhsOp, dotOp, rhsOp;
  if ( lhsOps.opReps_.size() > 0 ) lhsOp = lhsOps.opReps_.at(ilhs);
  if ( dotOps.opReps_.size() > 0 ) dotOp = dotOps.opReps_.at(idot);
  if ( rhsOps.opReps_.size() > 0 ) rhsOp = rhsOps.opReps_.at(irhs);

  if(SpinAdapted::dmrginp.doimplicitTranspose()==false){
    if(lhsOps.transpose_ || dotOps.transpose_ || rhsOps.transpose_)
    {
      pout << "Transposeview could not be used if m_implicitTranspose is false" <<endl;
      abort();
    }
  }
    
  // We need to distinguish cases where one or more blocks has an empty operator string
  // X_X_X
  if ( (lhsOps.opReps_.size() > 0) && (dotOps.opReps_.size() > 0) && (rhsOps.opReps_.size() > 0) ) {
    Transposeview lhsOpTr = Transposeview(*lhsOp);
    if ( lhsOps.transpose_ ) lhsOp = boost::shared_ptr<SparseMatrix>( &lhsOpTr, boostutils::null_deleter() );
    Transposeview dotOpTr = Transposeview(*dotOp);
    if ( dotOps.transpose_ ) dotOp = boost::shared_ptr<SparseMatrix>( &dotOpTr, boostutils::null_deleter() );
    Transposeview rhsOpTr = Transposeview(*rhsOp);
    if ( rhsOps.transpose_ ) rhsOp = boost::shared_ptr<SparseMatrix>( &rhsOpTr, boostutils::null_deleter() );
    expectation = spinExpectation(wavefunction_0, wavefunction_1, *lhsOp, *dotOp, *rhsOp, big_);
  }
  // X_X_0
  else if ( (lhsOps.opReps_.size() > 0) && (dotOps.opReps_.size() > 0) && (rhsOps.opReps_.size() == 0) ) {
    Transposeview lhsOpTr = Transposeview(*lhsOp);
    if ( lhsOps.transpose_ ) lhsOp = boost::shared_ptr<SparseMatrix>( &lhsOpTr, boostutils::null_deleter() );
    Transposeview dotOpTr = Transposeview(*dotOp);
    if ( dotOps.transpose_ ) dotOp = boost::shared_ptr<SparseMatrix>( &dotOpTr, boostutils::null_deleter() );
    expectation = spinExpectation(wavefunction_0, wavefunction_1, *lhsOp, *dotOp, *null, big_);
  }
  // X_0_X
  else if ( (lhsOps.opReps_.size() > 0) && (dotOps.opReps_.size() == 0) && (rhsOps.opReps_.size() > 0) ) {
    Transposeview lhsOpTr = Transposeview(*lhsOp);
    if ( lhsOps.transpose_ ) lhsOp = boost::shared_ptr<SparseMatrix>( &lhsOpTr, boostutils::null_deleter() );
    Transposeview rhsOpTr = Transposeview(*rhsOp);
    if ( rhsOps.transpose_ ) rhsOp = boost::shared_ptr<SparseMatrix>( &rhsOpTr, boostutils::null_deleter() );
    expectation = spinExpectation(wavefunction_0, wavefunction_1, *lhsOp, *null, *rhsOp, big_);
  }
  // 0_X_X
  else if ( (lhsOps.opReps_.size() == 0) && (dotOps.opReps_.size() > 0) && (rhsOps.opReps_.size() > 0) ) {
    Transposeview dotOpTr = Transposeview(*dotOp);
    if ( dotOps.transpose_ ) dotOp = boost::shared_ptr<SparseMatrix>( &dotOpTr, boostutils::null_deleter() );
    Transposeview rhsOpTr = Transposeview(*rhsOp);
    if ( rhsOps.transpose_ ) rhsOp = boost::shared_ptr<SparseMatrix>( &rhsOpTr, boostutils::null_deleter() );
    expectation = spinExpectation(wavefunction_0, wavefunction_1, *null, *dotOp, *rhsOp, big_);
  }
  // X_0_0
  else if ( (lhsOps.opReps_.size() > 0) && (dotOps.opReps_.size() == 0) && (rhsOps.opReps_.size() == 0) ) {
    Transposeview lhsOpTr = Transposeview(*lhsOp);
    if ( lhsOps.transpose_ ) lhsOp = boost::shared_ptr<SparseMatrix>( &lhsOpTr, boostutils::null_deleter() );
    expectation = spinExpectation(wavefunction_0, wavefunction_1, *lhsOp, *null, *null, big_);
  }
  // 0_X_0
  else if ( (lhsOps.opReps_.size() == 0) && (dotOps.opReps_.size() > 0) && (rhsOps.opReps_.size() == 0) ) {
    Transposeview dotOpTr = Transposeview(*dotOp);
    if ( dotOps.transpose_ ) dotOp = boost::shared_ptr<SparseMatrix>( &dotOpTr, boostutils::null_deleter() );
    expectation = spinExpectation(wavefunction_0, wavefunction_1, *null, *dotOp, *null, big_);
  }
  // 0_0_X
  else if ( (lhsOps.opReps_.size() == 0) && (dotOps.opReps_.size() == 0) && (rhsOps.opReps_.size() > 0) ) {
    Transposeview rhsOpTr = Transposeview(*rhsOp);
    if ( rhsOps.transpose_ ) rhsOp = boost::shared_ptr<SparseMatrix>( &rhsOpTr, boostutils::null_deleter() );
    expectation = spinExpectation(wavefunction_0, wavefunction_1, *null, *null, *rhsOp, big_);
  }
  else abort();

  // Modify new element with sign factors and return
  double factor = lhsOps.factor_ * dotOps.factor_ * rhsOps.factor_;
//pout << "expectation, factor = " << expectation << ", " << factor  << std::endl;
  return expectation*factor;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_expectations::test_for_singlet( int ilhs, int idot, int irhs, NpdmSpinOps_base& lhsOps, NpdmSpinOps_base& rhsOps, NpdmSpinOps_base& dotOps )
{
  int lhs2S, rhs2S, dot2S;
  // LHS
  if ( lhsOps.opReps_.size() == 0 ) 
    lhs2S = 0;
  else
    lhs2S = lhsOps.opReps_.at(ilhs)->get_deltaQuantum(0).get_s().getirrep();
  // RHS
  if ( rhsOps.opReps_.size() == 0 ) 
    rhs2S = 0;
  else
    rhs2S = rhsOps.opReps_.at(irhs)->get_deltaQuantum(0).get_s().getirrep();
  // DOT
  if ( dotOps.opReps_.size() == 0 ) 
    dot2S = 0;
  else
    dot2S = dotOps.opReps_.at(idot)->get_deltaQuantum(0).get_s().getirrep();

  // Couple LHS and Dot spin angular momenta and see if any equal RHS  
  for (int twoS = std::abs(lhs2S - dot2S); twoS <= ( lhs2S + dot2S ); twoS += 2 ) {
    if ( twoS == rhs2S ) {
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
  int hilhs = lhsOps.opReps_.size();
  int hirhs = rhsOps.opReps_.size();
  int hidot = dotOps.opReps_.size();
  for (int irhs = 0; irhs < std::max(1,hirhs); ++irhs) {
    for (int idot = 0; idot < std::max(1,hidot); ++idot) {
      for (int ilhs = 0; ilhs < std::max(1,hilhs); ++ilhs) {
        // Screen operator combinations that do not combine to give a singlet
        bool singlet = test_for_singlet( ilhs, idot, irhs, lhsOps, rhsOps, dotOps );
        // Contract operators to produce spin-adapted expectation value
        if ( singlet ) expectations_.push_back( contract_spin_adapted_operators( ilhs, idot, irhs, lhsOps, rhsOps, dotOps ) );
      }
    }
  }

  assert (expectations_.size() > 0);
//cout << "---------------------------------\n";
//cout << "spin-adapted expectations =\n";
////cout << "mpirank = " << mpigetrank() << endl;
//for (auto it = expectations_.begin(); it != expectations_.end(); ++it) {
//  cout << "p" << mpigetrank() << ":  " << *it << std::endl;
//}
//cout << "---------------------------------\n";

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

std::vector< std::pair< std::vector<int>, double > > 
Npdm_expectations::get_nonspin_adapted_expectations( NpdmSpinOps_base & lhsOps, NpdmSpinOps_base & rhsOps, NpdmSpinOps_base & dotOps )
{
  // Initialize dimension of spin-adapted to non-spin-adapted transformation
  int dim = 0;
  if(dmrginp.spinAdapted()){
  if ( npdm_order_ == 1 ) dim = 2;
  else if ( npdm_order_ == 2 ) dim = 6;
  else if ( npdm_order_ == 3 ) dim = 20;
  else if ( npdm_order_ == 4 ) dim = 70;
  else abort();
  }
  else{
    //if ( npdm_order_ == 1 ) dim = 1;
    //else if ( npdm_order_ == 2 ) dim = 1;
    //else if ( npdm_order_ == 3 ) dim = 20;
   // else if ( npdm_order_ == 4 ) dim = 70;
    //else abort();

  }
  std::vector< std::pair< std::vector<int>, double > > new_pdm_elements;

  // Get operator build string. e.g. (C2C4)(D5D6)
  std::string op_string;
  std::vector<int> indices;
  get_full_op_string( lhsOps, rhsOps, dotOps, op_string, indices );

  // Screen away unwanted strings (e.g. those that produce duplicate NPDM elements)
  //cout << op_string<<endl;
  //for( int i=0;i<indices.size();i++)
  //  cout << indices[i]<<',';
  //cout <<endl;
  if ( screen_op_string_for_duplicates( op_string, indices ) ) return new_pdm_elements;


  if(!dmrginp.spinAdapted()){
    std::vector<int> cd_order;
    int k=0;
    for ( auto it = op_string.begin(); it != op_string.end(); ++it ) {
      if (*it=='C') cd_order.push_back((indices[k++]));
      if (*it=='D') cd_order.push_back(1000+indices[k++]);
    }
    //sort cd_order; bubble sort
    int parity=1;
    double tmp;
    for( int i=cd_order.size()-1;i >0;i--)
      for(int j=0; j <i;j++){
        if(cd_order[j]>cd_order[j+1]){
          tmp=cd_order[j+1];
          cd_order[j+1]=cd_order[j];
          cd_order[j]=tmp;
          parity*=-1;
          
        }
      }

    //change cd_order back to spinorbital number
    //And check whether total S_z of operator is zero
    int sz=0;
    for(int i=0;i<cd_order.size()/2;i++){
      if(cd_order[i]%2)
        sz+=1;
      else 
        sz+=-1;

    }
    for(int i=cd_order.size()/2; i<cd_order.size() ;i++){
      assert(cd_order[i]>=1000);
      cd_order[i]-=1000;
      if(cd_order[i]%2)
        sz+=-1;
      else
        sz+=1;
    }
//    cout << "nonspinadapted spinorbital ordered: ";
//    for( int i=0;i<cd_order.size();i++)
//      cout << cd_order[i]<<',';
//    cout <<endl;
//    cout <<"sz: " <<sz<<endl;
    if(sz!=0) return new_pdm_elements;
    double nonspinvalue;

   // cout <<"size of operator: " <<lhsOps.opReps_.size()<<','<<rhsOps.opReps_.size()<<','<<dotOps.opReps_.size()<<endl;
      nonspinvalue=contract_spin_adapted_operators( 0, 0, 0, lhsOps, rhsOps, dotOps );

    new_pdm_elements.push_back(std::make_pair(cd_order,nonspinvalue*parity));
    return new_pdm_elements;
  }
  
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
  spin_adaptation_.npdm_set_up_linear_equations(dim, op_string, indices, expectations_, A, b, so_indices );

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
