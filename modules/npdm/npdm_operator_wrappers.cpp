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
//FIXME why do we need -1 here -- similar to (C(DD)) case
  factor_ = -1.0;
  transpose_ = false;
  build_pattern_ = "((CC)(DD))";
  // Build singlets only here
  mults_ = { 1, 1 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_compound_CCDD::set_local_ops( int idx )
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
  opReps_.push_back( build_compound_operator( false, -1, twoOps.at(0), twoOps.at(0), 0, indices_, true ) );
  // S=1 (+) S=1  =>  S=0
  opReps_.push_back( build_compound_operator( false, -1, twoOps.at(1), twoOps.at(1), 0, indices_, true ) );
  return false;
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

bool Npdm_op_wrapper_compound_CCD::set_local_ops( int idx )
{
pout << "getting compound CCD operator...\n";
assert(false);
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_CRE).get_local_element(idx);
  ix = twoOps.at(0)->get_orbs(0);
  jx = twoOps.at(0)->get_orbs(1);
  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(CRE).get_local_element(idx);
  kx = oneOp.at(0)->get_orbs(0);
//pout << "0 CC operator elements:\n";
//pout << *(twoOps[0]);
//pout << "1 CC operator elements:\n";
//pout << *(twoOps[1]);
//pout << "0 C operator elements:\n";
//pout << *(oneOp[0]);

  // Assumed single site (i=j=k)
  assert ( ix == jx );
  assert ( jx == kx );
  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
pout << "indices  " << ix << " " << jx << " " << kx << std::endl;

  opReps_.clear();
  // S=0 (+) S=1/2  =>  S=1/2
  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(0), oneOp.at(0), 0, indices_, true ) );
  // S=1 (+) S=1/2  =>  S=1/2
  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 0, indices_, true ) );
  // S=1 (+) S=1/2  =>  S=3/2
  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 1, indices_, true ) );
//pout << "2a CCD operator elements:\n";
//pout << *(opReps_[0]);
//pout << "2b CCD operator elements:\n";
//pout << *(opReps_[1]);
//pout << "4  CCD operator elements:\n";
//pout << *(opReps_[2]);
  return false;
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
//  build_pattern_ = "((CD)D)";
  // S={1/2,1/2,3/2}
  mults_ = { 2, 2, 4 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_compound_CDD::set_local_ops( int idx )
{
pout << "getting compound CDD operator...\n";
assert(false);
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_CRE).get_local_element(idx);
//  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_DES).get_local_element(idx);
  ix = twoOps.at(0)->get_orbs(0);
  jx = twoOps.at(0)->get_orbs(1);
  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(CRE).get_local_element(idx);
  kx = oneOp.at(0)->get_orbs(0);
//pout << "0 CD operator elements:\n";
//pout << *(twoOps[0]);
//pout << "1 CD operator elements:\n";
//pout << *(twoOps[1]);
//pout << "0 C operator elements:\n";
//pout << *(oneOp[0]);

  // Assumed single site (i=j=k)
  assert ( ix == jx );
  assert ( jx == kx );
  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
pout << "indices  " << ix << " " << jx << " " << kx << std::endl;

//-----------------------------
// 3 ways tested
//-----------------------------

//  "((CD)D)";
//  build_pattern_ = "((CD)D)";
//  opReps_.clear();
//  // S=0 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(0), oneOp.at(0), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=3/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 1, indices_, true ) );

//  "(C(DD))";
//  build_pattern_ = "(C(DD))";
//  factor_ = -1.0;  //FIXME note we need -1 here!
//  opReps_.clear();
//  // S=0 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, oneOp.at(0), twoOps.at(0), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, oneOp.at(0), twoOps.at(1), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=3/2
//  opReps_.push_back( build_compound_operator( true, -1, oneOp.at(0), twoOps.at(1), 1, indices_, true ) );

//  "((CC)D)" + tranpose;
//  "(C(DD))";
  transpose_ = true;
  build_pattern_ = "(C(DD))"; //FIXME note build pattern reflects structure *after* transpose is taken
  opReps_.clear();
  // S=0 (+) S=1/2  =>  S=1/2
  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(0), oneOp.at(0), 0, indices_, true ) );
  // S=1 (+) S=1/2  =>  S=1/2
  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 0, indices_, true ) );
  // S=1 (+) S=1/2  =>  S=3/2
  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 1, indices_, true ) );

// PRINT
//pout << "2a CDD operator elements:\n";
//pout << *(opReps_[0]);
//pout << "2b CDD operator elements:\n";
//pout << *(opReps_[1]);
//pout << "4  CDD operator elements:\n";
//pout << *(opReps_[2]);
  return false;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

Npdm_op_wrapper_compound_CDC::Npdm_op_wrapper_compound_CDC( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = 1;
  factor_ = 1.0;
  transpose_ = false;
//  build_pattern_ = "((CD)C)";
  // S={1/2,1/2,3/2}
  mults_ = { 2, 2, 4 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_compound_CDC::set_local_ops( int idx )
{
pout << "getting compound CDC operator...\n";
assert(false);
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_DES).get_local_element(idx);
  ix = twoOps.at(0)->get_orbs(0);
  jx = twoOps.at(0)->get_orbs(1);
  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(CRE).get_local_element(idx);
  kx = oneOp.at(0)->get_orbs(0);
//pout << "0 CD operator elements:\n";
//pout << *(twoOps[0]);
//pout << "1 CD operator elements:\n";
//pout << *(twoOps[1]);
//pout << "0 C operator elements:\n";
//pout << *(oneOp[0]);

  // Assumed single site (i=j=k)
  assert ( ix == jx );
  assert ( jx == kx );
  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
pout << "indices  " << ix << " " << jx << " " << kx << std::endl;

//----------
// 1st way
//----------
//
//  build_pattern_ = "((CD)C)";
//  opReps_.clear();
//  // S=0 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, 1, twoOps.at(0), oneOp.at(0), 0, indices_, false ) );
//  // S=1 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, 1, twoOps.at(1), oneOp.at(0), 0, indices_, false ) );
//  // S=1 (+) S=1/2  =>  S=3/2
//  opReps_.push_back( build_compound_operator( true, 1, twoOps.at(1), oneOp.at(0), 1, indices_, false ) );
//
//----------
// 2nd way
//----------
//
//  build_pattern_ = "(C(DC))";
//  factor_ = 1.0;
//  opReps_.clear();
//  //FIXME what should second argument be?
//  // S=0 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, oneOp.at(0), twoOps.at(0), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, oneOp.at(0), twoOps.at(1), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=3/2
//  opReps_.push_back( build_compound_operator( true, -1, oneOp.at(0), twoOps.at(1), 1, indices_, true ) );

//----------
// 3rd way
//----------

  // "((CD)D)" + transpose
  build_pattern_ = "(C(DC))";
  transpose_ = true;
  factor_ = 1.0;
  opReps_.clear();
  // S=0 (+) S=1/2  =>  S=1/2
  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(0), oneOp.at(0), 0, indices_, true ) );
  // S=1 (+) S=1/2  =>  S=1/2
  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 0, indices_, true ) );
  // S=1 (+) S=1/2  =>  S=3/2
  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 1, indices_, true ) );

//---------------------
// 1st way re-written
//---------------------
//
//  build_pattern_ = "((CD)C)";
//  factor_ = 1.0;
//  opReps_.clear();
//
//  // Singlet or triplet CD operator
//  boost::shared_ptr<SparseMatrix> cdOp (new CreDes);
//  cdOp->set_orbs() = indices_;
//  cdOp->set_orbs().pop_back();
//  cdOp->set_initialised() = true;
//  cdOp->set_fermion() = false;
//
//  // CD s=0
//  cdOp->set_deltaQuantum() = ( oneOp.at(0)->get_deltaQuantum() - oneOp.at(0)->get_deltaQuantum() ).at(0);
//  cdOp->allocate( spinBlock_->get_stateInfo() );
//  operatorfunctions::Product(spinBlock_, *oneOp.at(0), Transposeview(*oneOp.at(0)), *cdOp, 1.0 );
//pout << "0 CD operator elements (built by me):\n";
//pout << *cdOp;
//  // S=0 (+) S=1/2  =>  S=1/2
//  boost::shared_ptr<SparseMatrix> cdcOp (new Cre);
//  cdcOp->set_orbs() = indices_;
//  cdcOp->set_initialised() = true;
//  cdcOp->set_fermion() = true;
//  cdcOp->set_deltaQuantum() = ( cdOp->get_deltaQuantum() + oneOp.at(0)->get_deltaQuantum() ).at(0);
//  cdcOp->allocate( spinBlock_->get_stateInfo() );
//  operatorfunctions::Product(spinBlock_, *cdOp, *oneOp.at(0), *cdcOp, 1.0 );
//  opReps_.push_back( cdcOp );
////pout << *cdcOp;
//
//  // CD s=1
//  cdOp->set_deltaQuantum() = ( oneOp.at(0)->get_deltaQuantum() - oneOp.at(0)->get_deltaQuantum() ).at(1);
//  cdOp->allocate( spinBlock_->get_stateInfo() );
//  operatorfunctions::Product(spinBlock_, *oneOp.at(0), Transposeview(*oneOp.at(0)), *cdOp, 1.0 );
//pout << "1 CD operator elements (built by me):\n";
//pout << *cdOp;
//  // S=1 (+) S=1/2  =>  S=1/2
//  boost::shared_ptr<SparseMatrix> cdcOp2 (new Cre);
//  cdcOp2->set_orbs() = indices_;
//  cdcOp2->set_initialised() = true;
//  cdcOp2->set_fermion() = true;
//  cdcOp2->set_deltaQuantum() = ( cdOp->get_deltaQuantum() + oneOp.at(0)->get_deltaQuantum() ).at(0);
//  cdcOp2->allocate( spinBlock_->get_stateInfo() );
//  operatorfunctions::Product(spinBlock_, *cdOp, *oneOp.at(0), *cdcOp2, 1.0 );
//  opReps_.push_back( cdcOp2 );
//  // S=1 (+) S=1/2  =>  S=3/2
//  boost::shared_ptr<SparseMatrix> cdcOp3 (new Cre);
//  cdcOp3->set_orbs() = indices_;
//  cdcOp3->set_initialised() = true;
//  cdcOp3->set_fermion() = true;
//  cdcOp3->set_deltaQuantum() = ( cdOp->get_deltaQuantum() + oneOp.at(0)->get_deltaQuantum() ).at(1);
//  cdcOp3->allocate( spinBlock_->get_stateInfo() );
//  operatorfunctions::Product(spinBlock_, *cdOp, *oneOp.at(0), *cdcOp3, 1.0 );
//  opReps_.push_back( cdcOp3 );

// PRINT
//pout << "CDC operator elements (built by me):\n";
//pout << "2a CDC operator elements:\n";
//pout << *(opReps_[0]);
//pout << "2b CDC operator elements:\n";
//pout << *(opReps_[1]);
//pout << "4  CDC operator elements:\n";
//pout << *(opReps_[2]);
  return false;
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
//  build_pattern_ = "((CC)D)";
  build_pattern_ = "0";
  // S={1/2,1/2,3/2}
  mults_ = { 2, 2, 4 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CCD::set_local_ops( int idx )
{
pout << "getting CCD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  opReps_ = spinBlock_->get_op_array(CRE_CRE_DES).get_local_element(idx);
  build_pattern_ = opReps_.at(0)->get_build_pattern();
  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
pout << "indices  " << ix << " " << jx << " " << kx << std::endl;
//pout << "build pattern " << opReps_.at(0)->get_build_pattern() << std::endl;
//pout << "2a CCD operator elements:\n";
//pout << *(opReps_[0]);
//pout << "2b CCD operator elements:\n";
//pout << *(opReps_[1]);
//pout << "4  CCD operator elements:\n";
//pout << *(opReps_[2]);
//  assert( build_pattern_ == opReps_.at(0)->get_build_pattern() );

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_CDD::Npdm_op_wrapper_CDD( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_DES_DES).get_size();
  factor_ = 0.0;
  transpose_ = false;
//  build_pattern_ = "((CD)D)";
  build_pattern_ = "0";
  // S={1/2,1/2,3/2}
  mults_ = { 2, 2, 4 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CDD::set_local_ops( int idx )
{
pout << "getting CDD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  opReps_ = spinBlock_->get_op_array(CRE_DES_DES).get_local_element(idx);
  build_pattern_ = opReps_.at(0)->get_build_pattern();
  if ( build_pattern_ == "(C(DD))" ) factor_ = -1.0;
  else if ( build_pattern_ == "((CD)D)" ) factor_ = 1.0;
  else assert (false);

  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
pout << "indices  " << ix << " " << jx << " " << kx << std::endl;
//pout << "build pattern " << opReps_.at(0)->get_build_pattern() << std::endl;
//pout << "2a CDD operator elements:\n";
//pout << *(opReps_[0]);
//pout << "2b CDD operator elements:\n";
//pout << *(opReps_[1]);
//pout << "4  CDD operator elements:\n";
//pout << *(opReps_[2]);
//  assert( build_pattern_ == opReps_.at(0)->get_build_pattern() );

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
//if ( (ix == jx) && (jx == kx) ) return false;
//pout << "WARNING: skipping this operator\n";
//return true;
  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_CDC::Npdm_op_wrapper_CDC( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_DES_CRE).get_size();
  factor_ = 1.0;
  transpose_ = false;
//  build_pattern_ = "((CD)C)";
  build_pattern_ = "0";
  // S={1/2,1/2,3/2}
  mults_ = { 2, 2, 4 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CDC::set_local_ops( int idx )
{
pout << "getting CDC operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  opReps_ = spinBlock_->get_op_array(CRE_DES_CRE).get_local_element(idx);
  build_pattern_ = opReps_.at(0)->get_build_pattern();
  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
pout << "indices  " << ix << " " << jx << " " << kx << std::endl;
//pout << "build pattern " << opReps_.at(0)->get_build_pattern() << std::endl;
//pout << "2a CDC operator elements:\n";
//pout << *(opReps_[0]);
//pout << "2b CDC operator elements:\n";
//pout << *(opReps_[1]);
//pout << "4  CDC operator elements:\n";
//pout << *(opReps_[2]);
//  assert( build_pattern_ == opReps_.at(0)->get_build_pattern() );

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
//FIXME fails!!!
if ( (ix == jx) && (jx == kx) ) {
//if ( jx == kx ) {
  pout << "WARNING: skipping this operator\n";
  return true;
}
  return false;
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

bool Npdm_op_wrapper_CC::set_local_ops( int idx )
{
pout << "getting CC operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx;

  opReps_ = spinBlock_->get_op_array(CRE_CRE).get_local_element(idx);
  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
pout << "indices  " << ix << " " << jx << std::endl;
//pout << "0 CC operator elements:\n";
//pout << *(opReps_[0]);
//pout << "1 CC operator elements:\n";
//pout << *(opReps_[1]);

  indices_.push_back( ix );
  indices_.push_back( jx );
  return false;
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

bool Npdm_op_wrapper_CD::set_local_ops( int idx )
{
pout << "getting CD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx;

  opReps_ = spinBlock_->get_op_array(CRE_DES).get_local_element(idx);
  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
pout << "indices  " << ix << " " << jx << std::endl;
//pout << "0 CD operator elements:\n";
//pout << *(opReps_[0]);
//pout << "1 CD operator elements:\n";
//pout << *(opReps_[1]);

  transpose_ = false;
  indices_.push_back( ix );
  indices_.push_back( jx );
  return false;
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

bool Npdm_op_wrapper_DD::set_local_ops( int idx )
{
pout << "getting DD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx;

  opReps_ = spinBlock_->get_op_array(CRE_CRE).get_local_element(idx);
  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
pout << "indices  " << ix << " " << jx << std::endl;

  // Note use of transpose means we store this as (j,i) not (i,j)
  indices_.push_back( jx );
  indices_.push_back( ix );
  return false;
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

bool Npdm_op_wrapper_C::set_local_ops( int idx )
{
pout << "getting C operator...\n";
  indices_.clear();
  opReps_ = spinBlock_->get_op_array(CRE).get_local_element(idx);
  int ix = opReps_.at(0)->get_orbs(0);
  indices_.push_back(ix);
  return false;
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

bool Npdm_op_wrapper_D::set_local_ops( int idx )
{
pout << "getting D operator...\n";
  indices_.clear();
  opReps_ = spinBlock_->get_op_array(CRE).get_local_element(idx);
  int ix = opReps_.at(0)->get_orbs(0);
  indices_.push_back(ix);
pout << "indices  " << ix << std::endl;
//pout << "D operator elements:\n";
//pout << *(opReps_[0]);
  return false;
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

bool Npdm_op_wrapper_NULL::set_local_ops( int idx )
{
  // Null operator
  indices_.clear();
  opReps_.clear();
  return false;
}


//===========================================================================================================================================================
// NOTE transpose applied to RHS operator here

boost::shared_ptr<SparseMatrix> NpdmSpinOps::build_compound_operator( bool is_fermion, int sign,
                                                                      boost::shared_ptr<SparseMatrix> lhsOp,
                                                                      boost::shared_ptr<SparseMatrix> rhsOp,
                                                                      int ispin, std::vector<int> indices, bool transpose )
{
  // Initialize new operator
//FIXME new BROKEN!
  boost::shared_ptr<SparseMatrix> newOp (new Cre);
  assert( lhsOp->get_orbs().size() + rhsOp->get_orbs().size() == indices.size() );
  newOp->set_orbs() = indices;
  newOp->set_initialised() = true;
  newOp->set_fermion() = is_fermion;
  if (sign == 1) newOp->set_deltaQuantum() = ( lhsOp->get_deltaQuantum() + rhsOp->get_deltaQuantum() ).at(ispin);
  else newOp->set_deltaQuantum() = ( lhsOp->get_deltaQuantum() - rhsOp->get_deltaQuantum() ).at(ispin);
//pout << "Building compound operator.......................................................\n";
//pout << "2*lhs spin = " << lhsOp->get_deltaQuantum().get_s() << std::endl;
//pout << "2*rhs spin = " << rhsOp->get_deltaQuantum().get_s() << std::endl;
//pout << "2*total spin = " << newOp->get_deltaQuantum().get_s() << std::endl;
  newOp->allocate( spinBlock_->get_stateInfo() );
//pout << *newOp;

  if (transpose) {
    // Build compound operator as product of LHS and TRANSPOSE( RHS )
    operatorfunctions::Product(spinBlock_, *lhsOp, Transposeview(*rhsOp), *newOp, 1.0 );
  }
  else {
    // Build compound operator as product of LHS and RHS
    operatorfunctions::Product(spinBlock_, *lhsOp, *rhsOp, *newOp, 1.0 );
  }
//pout << "Done Building compound operator.......................................................\n";

  return newOp;
}

//===========================================================================================================================================================

}
