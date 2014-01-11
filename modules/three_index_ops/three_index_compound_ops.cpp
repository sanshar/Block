/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "global.h"
#include "npdm_patterns.h"
#include "npdm_operator_wrappers.h"

namespace SpinAdapted{

//===========================================================================================================================================================

Npdm_op_wrapper_compound_CCD::Npdm_op_wrapper_compound_CCD( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
cout << "init compound CCD operator...\n";
  size_ = spinBlock_->get_op_array(RI_3INDEX).get_size();
  is_local_ = true;
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((CC)D)";
  // S={1/2,1/2,3/2}
  mults_ = { 2, 2, 4 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_compound_CCD::set_local_ops( int idx )
{
cout << "getting compound CCD operator...\n";
  // Spatial orbital indices
//FIXME don't need to keep reconstructing whole array!
  indices_ = spinBlock_->get_op_array(RI_3INDEX).get_array().at(idx);
  int ix = indices_[0];
  int jx = indices_[1];
  int kx = indices_[2];

  // Get 2-index and 1-index ops as RI building blocks
  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_CRE).get_element(ix,jx);
  assert( ix == twoOps.at(0)->get_orbs(0) );
  assert( jx == twoOps.at(0)->get_orbs(1) );
  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(CRE).get_element(kx);
  assert( kx == oneOp.at(0)->get_orbs(0) );

  // Allocate and build operator representation on the fly as RI tensor product for each spin component
  opReps_.clear();
  // S=0 (+) S=1/2  =>  S=1/2
  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(0), oneOp.at(0), 0, indices_, true ) );
  // S=1 (+) S=1/2  =>  S=1/2
  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 0, indices_, true ) );
  // S=1 (+) S=1/2  =>  S=3/2
  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 1, indices_, true ) );
  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_compound_CDD::Npdm_op_wrapper_compound_CDD( SpinBlock * spinBlock )
{
  // Assume single site block
//FIXME update to allow RI with any orbitals as above for CCD
  assert( spinBlock->size() == 1 );
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = 1;
  is_local_ = true;
  factor_ = 1.0;
  transpose_ = false;
//  build_pattern_ = "((CD)D)";
  build_pattern_ = "0";
  // S={1/2,1/2,3/2}
  mults_ = { 2, 2, 4 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_compound_CDD::set_local_ops( int idx )
{
//cout << "getting compound CDD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

//  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_CRE).get_local_element(idx);
  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_DES).get_local_element(idx);
  ix = twoOps.at(0)->get_orbs(0);
  jx = twoOps.at(0)->get_orbs(1);
  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(CRE).get_local_element(idx);
  kx = oneOp.at(0)->get_orbs(0);

  // Assumed single site (i=j=k)
  assert ( ix == jx );
  assert ( jx == kx );
  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );

//-----------------------------
// 3 ways tested
//-----------------------------

//  "((CD)D)";
//  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_DES).get_local_element(idx);
  build_pattern_ = "((CD)D)";
  opReps_.clear();
  // S=0 (+) S=1/2  =>  S=1/2
  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(0), oneOp.at(0), 0, indices_, true ) );
  // S=1 (+) S=1/2  =>  S=1/2
  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 0, indices_, true ) );
  // S=1 (+) S=1/2  =>  S=3/2
  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 1, indices_, true ) );

//  "(C(DD))";
//  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_CRE).get_local_element(idx);
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
//  transpose_ = true;
//  build_pattern_ = "(C(DD))"; //FIXME note build pattern reflects structure *after* transpose is taken
//  opReps_.clear();
//  // S=0 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(0), oneOp.at(0), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=3/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 1, indices_, true ) );

// PRINT
//cout << "2a CDD operator elements:\n";
//cout << *(opReps_[0]);
//cout << "2b CDD operator elements:\n";
//cout << *(opReps_[1]);
//cout << "4  CDD operator elements:\n";
//cout << *(opReps_[2]);
  return false;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

Npdm_op_wrapper_compound_CDC::Npdm_op_wrapper_compound_CDC( SpinBlock * spinBlock )
{
assert(false);
  // Assume single site block
  assert( spinBlock->size() == 1 );
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = 1;
  is_local_ = true;
  factor_ = 1.0;
  transpose_ = false;
//  build_pattern_ = "((CD)C)";
  // S={1/2,1/2,3/2}
  mults_ = { 2, 2, 4 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_compound_CDC::set_local_ops( int idx )
{
//cout << "getting compound CDC operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_DES).get_local_element(idx);
  ix = twoOps.at(0)->get_orbs(0);
  jx = twoOps.at(0)->get_orbs(1);
  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(CRE).get_local_element(idx);
  kx = oneOp.at(0)->get_orbs(0);
//cout << "singlet CD operator elements:\n";
//cout << *(twoOps[0]);
//cout << "triplet CD operator elements:\n";
//cout << *(twoOps[1]);
//cout << "half C operator elements:\n";
//cout << *(oneOp[0]);

  // Assumed single site (i=j=k)
  assert ( ix == jx );
  assert ( jx == kx );
  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
//cout << "indices  " << ix << " " << jx << " " << kx << std::endl;

//----------
// 1st way
//----------

  build_pattern_ = "((CD)C)";
  opReps_.clear();
  // S=0 (+) S=1/2  =>  S=1/2
  opReps_.push_back( build_compound_operator( true, 1, twoOps.at(0), oneOp.at(0), 0, indices_, false ) );
  // S=1 (+) S=1/2  =>  S=1/2
  opReps_.push_back( build_compound_operator( true, 1, twoOps.at(1), oneOp.at(0), 0, indices_, false ) );
  // S=1 (+) S=1/2  =>  S=3/2
  opReps_.push_back( build_compound_operator( true, 1, twoOps.at(1), oneOp.at(0), 1, indices_, false ) );

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

//  // "((CD)D)" + transpose
//  build_pattern_ = "(C(DC))";
//  transpose_ = true;
//  factor_ = 1.0;
//  opReps_.clear();
//  // S=0 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(0), oneOp.at(0), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=3/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 1, indices_, true ) );

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
//cout << "singlet CD operator elements (built by me):\n";
//cout << *cdOp;
//  // S=0 (+) S=1/2  =>  S=1/2
//  boost::shared_ptr<SparseMatrix> cdcOp (new Cre);
//  cdcOp->set_orbs() = indices_;
//  cdcOp->set_initialised() = true;
//  cdcOp->set_fermion() = true;
//  cdcOp->set_deltaQuantum() = ( cdOp->get_deltaQuantum() + oneOp.at(0)->get_deltaQuantum() ).at(0);
//  cdcOp->allocate( spinBlock_->get_stateInfo() );
//  operatorfunctions::Product(spinBlock_, *cdOp, *oneOp.at(0), *cdcOp, 1.0 );
//  opReps_.push_back( cdcOp );
////cout << *cdcOp;
//
//  // CD s=1
//  cdOp->set_deltaQuantum() = ( oneOp.at(0)->get_deltaQuantum() - oneOp.at(0)->get_deltaQuantum() ).at(1);
//  cdOp->allocate( spinBlock_->get_stateInfo() );
//  operatorfunctions::Product(spinBlock_, *oneOp.at(0), Transposeview(*oneOp.at(0)), *cdOp, 1.0 );
//cout << "triplet CD operator elements (built by me):\n";
//cout << *cdOp;
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
//cout << "CDC operator elements (built by me):\n";
//cout << "2a CDC operator elements:\n";
//cout << *(opReps_[0]);
//cout << "2b CDC operator elements:\n";
//cout << *(opReps_[1]);
//cout << "4  CDC operator elements:\n";
//cout << *(opReps_[2]);

  if ( jx == kx ) {
    //FIXME I think this fails because of potential problems commuting operators with same indices in spin-transformation
    //cout << "WARNING: skipping this operator\n";
    return true;
  }
  return false;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

Npdm_op_wrapper_compound_CCC::Npdm_op_wrapper_compound_CCC( SpinBlock * spinBlock )
{
  // Assume single site block
assert(false); // << this operator should always be zero on one site!
  assert( spinBlock->size() == 1 );
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = 1;
  is_local_ = true;
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((CC)C)";
  // S={1/2,1/2,3/2}
  mults_ = { 2, 2, 4 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_compound_CCC::set_local_ops( int idx )
{
//cout << "getting compound CCC operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_CRE).get_local_element(idx);
  ix = twoOps.at(0)->get_orbs(0);
  jx = twoOps.at(0)->get_orbs(1);
  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(CRE).get_local_element(idx);
  kx = oneOp.at(0)->get_orbs(0);

  // Assumed single site (i=j=k)
  assert ( ix == jx );
  assert ( jx == kx );
  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );

  opReps_.clear();
  // S=0 (+) S=1/2  =>  S=1/2
  opReps_.push_back( build_compound_operator( true, 1, twoOps.at(0), oneOp.at(0), 0, indices_, false ) );
  // S=1 (+) S=1/2  =>  S=1/2
  opReps_.push_back( build_compound_operator( true, 1, twoOps.at(1), oneOp.at(0), 0, indices_, false ) );
  // S=1 (+) S=1/2  =>  S=3/2
  opReps_.push_back( build_compound_operator( true, 1, twoOps.at(1), oneOp.at(0), 1, indices_, false ) );

  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_compound_DCD::Npdm_op_wrapper_compound_DCD( SpinBlock * spinBlock )
{
assert(false);
  // Assume single site block
  assert( spinBlock->size() == 1 );
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = 1;
  is_local_ = true;
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((DC)D)";
  // S={1/2,1/2,3/2}
  mults_ = { 2, 2, 4 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_compound_DCD::set_local_ops( int idx )
{
//cout << "getting compound DCD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(DES_CRE).get_local_element(idx);
  ix = twoOps.at(0)->get_orbs(0);
  jx = twoOps.at(0)->get_orbs(1);
  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(CRE).get_local_element(idx);
  kx = oneOp.at(0)->get_orbs(0);

  // Assumed single site (i=j=k)
  assert ( ix == jx );
  assert ( jx == kx );
  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );

  opReps_.clear();
  // S=0 (+) S=1/2  =>  S=1/2
  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(0), oneOp.at(0), 0, indices_, true ) );
  // S=1 (+) S=1/2  =>  S=1/2
  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 0, indices_, true ) );
  // S=1 (+) S=1/2  =>  S=3/2
  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 1, indices_, true ) );

  if ( ix == jx ) {
    //FIXME I think this fails because of potential problems commuting operators with same indices in spin-transformation
    //cout << "WARNING: skipping this operator\n";
    return true;
  }
  return false;

}

//===========================================================================================================================================================

Npdm_op_wrapper_compound_DDC::Npdm_op_wrapper_compound_DDC( SpinBlock * spinBlock )
{
  // Assume single site block
assert(false);
//  assert( spinBlock->size() == 1 );
//  opReps_.clear();
//  indices_.clear();
//  spinBlock_ = spinBlock;
//  size_ = 1;
//  is_local_ = true;
//  factor_ = 1.0;
//  transpose_ = false;
//  build_pattern_ = "((DC)D)";
//  // S={1/2,1/2,3/2}
//  mults_ = { 2, 2, 4 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_compound_DDC::set_local_ops( int idx )
{
//cout << "getting compound DDC operator...\n";
  // Spatial orbital indices
assert(false);
//  indices_.clear();
//  int ix, jx, kx;
//
//  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(DES_CRE).get_local_element(idx);
//  ix = twoOps.at(0)->get_orbs(0);
//  jx = twoOps.at(0)->get_orbs(1);
//  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(CRE).get_local_element(idx);
//  kx = oneOp.at(0)->get_orbs(0);
//
//  // Assumed single site (i=j=k)
//  assert ( ix == jx );
//  assert ( jx == kx );
//  indices_.push_back( ix );
//  indices_.push_back( jx );
//  indices_.push_back( kx );
//
//  opReps_.clear();
//  // S=0 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(0), oneOp.at(0), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=3/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 1, indices_, true ) );
//
//  if ( ix == jx ) {
//    //FIXME I think this fails because of potential problems commuting operators with same indices in spin-transformation
//    //cout << "WARNING: skipping this operator\n";
//    return true;
//  }
//  return false;
//
}

//===========================================================================================================================================================

}

