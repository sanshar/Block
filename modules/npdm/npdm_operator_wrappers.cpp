/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "npdm_patterns.h"
#include "npdm_operator_wrappers.h"

namespace SpinAdapted{

//===========================================================================================================================================================
// 4-INDEX COMPOUND OPERATORS (for one site only)
//===========================================================================================================================================================

Npdm_op_wrapper_compound_CCDD::Npdm_op_wrapper_compound_CCDD( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = 1;
//FIXME why do we need -1 here ??
  factor_ = -1.0;
  transpose_ = false;
  build_pattern_ = "((CC)(DD))";
  // Build singlets only here
  mults_ = { 1, 1 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_op_wrapper_compound_CCDD::set_local_ops( int idx )
{
pout << "getting CCDD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx, lx;
  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_CRE).get_local_element(idx);
  ix = twoOps.at(0)->get_orbs(0);
  jx = twoOps.at(0)->get_orbs(1);

  // Assumed single site
  assert ( ix == jx );
  indices_.push_back( ix );
  indices_.push_back( ix );
  indices_.push_back( ix );
  indices_.push_back( ix );

  opReps_.clear();
  // S=0 (+) S=0  =>  S=0
//FIXME why NOT FERMION????
  opReps_.push_back( build_compound_operator( false, twoOps.at(0), twoOps.at(0), 0, indices_, true ) );
  // S=1 (+) S=1  =>  S=0
//FIXME why NOT FERMION????
  opReps_.push_back( build_compound_operator( false, twoOps.at(1), twoOps.at(1), 0, indices_, true ) );
}

//===========================================================================================================================================================
// 3-INDEX COMPOUND OPERATORS (for one site only)
//===========================================================================================================================================================

Npdm_op_wrapper_compound_CCD::Npdm_op_wrapper_compound_CCD( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = 1;
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((CC)D)";
  // S={1/2,1/2,3/2}
  mults_ = { 2, 2, 4 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_op_wrapper_compound_CCD::set_local_ops( int idx )
{
pout << "getting compound CCD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_CRE).get_local_element(idx);
pout << "0 CC operator elements:\n";
pout << *(twoOps[0]);
pout << "1 CC operator elements:\n";
pout << *(twoOps[1]);
  ix = twoOps.at(0)->get_orbs(0);
  jx = twoOps.at(0)->get_orbs(1);
  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(CRE).get_local_element(idx);
  kx = oneOp.at(0)->get_orbs(0);

  // Assumed single site
  assert ( ix == jx );
  assert ( jx == kx );
  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
pout << "indices  " << ix << " " << jx << " " << kx << std::endl;

  opReps_.clear();
  // S=0 (+) S=1/2  =>  S=1/2
  opReps_.push_back( build_compound_operator( true, twoOps.at(0), oneOp.at(0), 0, indices_, true ) );
  // S=1 (+) S=1/2  =>  S=1/2
  opReps_.push_back( build_compound_operator( true, twoOps.at(1), oneOp.at(0), 0, indices_, true ) );
  // S=1 (+) S=1/2  =>  S=3/2
  opReps_.push_back( build_compound_operator( true, twoOps.at(1), oneOp.at(0), 1, indices_, true ) );
pout << "2a CCD operator elements:\n";
pout << *(opReps_[0]);
pout << "2b CCD operator elements:\n";
pout << *(opReps_[1]);
pout << "4  CCD operator elements:\n";
pout << *(opReps_[2]);
}

//===========================================================================================================================================================

Npdm_op_wrapper_compound_CDD::Npdm_op_wrapper_compound_CDD( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = 1;
  factor_ = 1.0;
  transpose_ = false;
//  build_pattern_ = "(C(DD))";
  build_pattern_ = "((CD)D)";
  // S={1/2,1/2,3/2}
  mults_ = { 2, 2, 4 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_op_wrapper_compound_CDD::set_local_ops( int idx )
{
pout << "getting compound CDD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

//  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_CRE).get_local_element(idx);
  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_DES).get_local_element(idx);
  ix = twoOps.at(0)->get_orbs(0);
  jx = twoOps.at(0)->get_orbs(1);
  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(CRE).get_local_element(idx);
  kx = oneOp.at(0)->get_orbs(0);

  // Assumed single site
  assert ( ix == jx );
  assert ( jx == kx );
  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );

//  "((CD)D)";
  opReps_.clear();
  // S=0 (+) S=1/2  =>  S=1/2
  opReps_.push_back( build_compound_operator( true, twoOps.at(0), oneOp.at(0), 0, indices_, true ) );
  // S=1 (+) S=1/2  =>  S=1/2
  opReps_.push_back( build_compound_operator( true, twoOps.at(1), oneOp.at(0), 0, indices_, true ) );
  // S=1 (+) S=1/2  =>  S=3/2
  opReps_.push_back( build_compound_operator( true, twoOps.at(1), oneOp.at(0), 1, indices_, true ) );

//  "(C(DD))";
//FIXME this fails... WHY??
//FIXME why do we need -1 here to make some work??
//  factor_ = -1.0;
//  opReps_.clear();
//  // S=0 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, oneOp, 0, twoOps, 0, 0, ix, true ) );
//  // S=1 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, oneOp, 0, twoOps, 1, 0, ix, true ) );
//  // S=1 (+) S=1/2  =>  S=3/2
//  opReps_.push_back( build_compound_operator( true, oneOp, 0, twoOps, 1, 1, ix, true ) );

//  "(C(DD))";
//FIXME this + transpose=true fails... WHY??
//  opReps_.clear();
//  // S=0 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, twoOps, 0, oneOp, 0, 0, ix, true ) );
//  // S=1 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, twoOps, 1, oneOp, 0, 0, ix, true ) );
//  // S=1 (+) S=1/2  =>  S=3/2
//  opReps_.push_back( build_compound_operator( true, twoOps, 1, oneOp, 0, 1, ix, true ) );
}

//===========================================================================================================================================================
// 3-INDEX OPERATORS
//===========================================================================================================================================================

Npdm_op_wrapper_CCD::Npdm_op_wrapper_CCD( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_CRE_DES).get_size();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((CC)D)";
  // S={1/2,1/2,3/2}
  mults_ = { 2, 2, 4 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_op_wrapper_CCD::set_local_ops( int idx )
{
pout << "getting CCD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  opReps_ = spinBlock_->get_op_array(CRE_CRE_DES).get_local_element(idx);
  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
pout << "indices  " << ix << " " << jx << " " << kx << std::endl;
pout << "2a CCD operator elements:\n";
pout << *(opReps_[0]);
pout << "2b CCD operator elements:\n";
pout << *(opReps_[1]);
pout << "4  CCD operator elements:\n";
pout << *(opReps_[2]);

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
}

//===========================================================================================================================================================

Npdm_op_wrapper_CDD::Npdm_op_wrapper_CDD( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_CRE_DES).get_size();
  factor_ = 1.0;
  transpose_ = false;
//FIXME
  build_pattern_ = "((CC)D)";
  // S={1/2,1/2,3/2}
  mults_ = { 2, 2, 4 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_op_wrapper_CDD::set_local_ops( int idx )
{
pout << "getting CDD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;
//FIXME
assert(false);

  opReps_ = spinBlock_->get_op_array(CRE_CRE_DES).get_local_element(idx);
  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
}

//===========================================================================================================================================================
// 2-INDEX OPERATORS
//===========================================================================================================================================================

Npdm_op_wrapper_CC::Npdm_op_wrapper_CC( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_CRE).get_size();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "(CC)";
  // S={0,1}
  mults_ = { 1, 3 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_op_wrapper_CC::set_local_ops( int idx )
{
pout << "getting CC operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx;

  opReps_ = spinBlock_->get_op_array(CRE_CRE).get_local_element(idx);
  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);

  indices_.push_back( ix );
  indices_.push_back( jx );
}

//===========================================================================================================================================================

Npdm_op_wrapper_CD::Npdm_op_wrapper_CD( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_DES).get_size();
  factor_ = 1.0;
  build_pattern_ = "(CD)";
  // S={0,1}
  mults_ = { 1, 3 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_op_wrapper_CD::set_local_ops( int idx )
{
pout << "getting CD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx;

  opReps_ = spinBlock_->get_op_array(CRE_DES).get_local_element(idx);
  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);

  transpose_ = false;
  indices_.push_back( ix );
  indices_.push_back( jx );
}

//===========================================================================================================================================================

Npdm_op_wrapper_DD::Npdm_op_wrapper_DD( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_CRE).get_size();
//FIXME why do we need -1 here ??
  factor_ = -1.0;
  transpose_ = true;
  build_pattern_ = "(DD)";
  // S={0,1}
  mults_ = { 1, 3 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_op_wrapper_DD::set_local_ops( int idx )
{
pout << "getting DD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx;

  opReps_ = spinBlock_->get_op_array(CRE_CRE).get_local_element(idx);
  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);

  // Note use of transpose means we store this as (j,i) not (i,j)
  indices_.push_back( jx );
  indices_.push_back( ix );
}

//===========================================================================================================================================================
// 1-INDEX OPERATORS
//===========================================================================================================================================================

Npdm_op_wrapper_C::Npdm_op_wrapper_C( SpinBlock * spinBlock ) 
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE).get_size();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "(C)";
  // S=1/2 only
  mults_ = { 2 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_op_wrapper_C::set_local_ops( int idx )
{
pout << "getting C operator...\n";
  indices_.clear();
  opReps_ = spinBlock_->get_op_array(CRE).get_local_element(idx);
  int ix = opReps_.at(0)->get_orbs(0);
  indices_.push_back(ix);
}

//===========================================================================================================================================================

Npdm_op_wrapper_D::Npdm_op_wrapper_D( SpinBlock * spinBlock ) 
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE).get_size();
  factor_ = 1.0;
  transpose_ = true;
  build_pattern_ = "(D)";
  // S=1/2 only
  mults_ = { 2 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_op_wrapper_D::set_local_ops( int idx )
{
pout << "getting D operator...\n";
  indices_.clear();
  opReps_ = spinBlock_->get_op_array(CRE).get_local_element(idx);
  int ix = opReps_.at(0)->get_orbs(0);
  indices_.push_back(ix);
pout << "indices  " << ix << std::endl;
pout << "D operator elements:\n";
pout << *(opReps_[0]);
}

//===========================================================================================================================================================

Npdm_op_wrapper_NULL::Npdm_op_wrapper_NULL()
{
  // Null operator
  opReps_.clear();
  indices_.clear();
  size_ = 1; // For compatibility with rest of code
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "";
  mults_ = { 1 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_op_wrapper_NULL::set_local_ops( int idx )
{
  // Null operator
  indices_.clear();
  opReps_.clear();
}


//===========================================================================================================================================================
// NOTE transpose applied to RHS operator here

boost::shared_ptr<SparseMatrix> NpdmSpinOps::build_compound_operator( bool is_fermion,
                                                                      boost::shared_ptr<SparseMatrix> lhsOp,
                                                                      boost::shared_ptr<SparseMatrix> rhsOp,
                                                                      int ispin, std::vector<int> indices, bool transpose )
{
  // Initialize new operator
  boost::shared_ptr<SparseMatrix> newOp (new Cre);
  assert( lhsOp->get_orbs().size() + rhsOp->get_orbs().size() == indices.size() );
  newOp->set_orbs() = indices;
  newOp->set_initialised() = true;
  newOp->set_fermion() = is_fermion;
  newOp->set_deltaQuantum() = ( lhsOp->get_deltaQuantum() - rhsOp->get_deltaQuantum() ).at(ispin);
  newOp->allocate( spinBlock_->get_stateInfo() );

  if (transpose) {
    // Build compound operator as product of LHS and TRANSPOSE( RHS )
    operatorfunctions::Product(spinBlock_, *lhsOp, Transposeview(*rhsOp), *newOp, 1.0 );
  }
  else {
    // Build compound operator as product of LHS and RHS
    operatorfunctions::Product(spinBlock_, *lhsOp, *rhsOp, *newOp, 1.0 );
  }

  return newOp;
}

//===========================================================================================================================================================

}
