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
// 3-INDEX COMPOUND OPERATORS (for one site only)
//===========================================================================================================================================================

Npdm_op_wrapper_compound_CCD::Npdm_op_wrapper_compound_CCD( SpinBlock * spinBlock )
{
  // Assume single site block
  assert( spinBlock->size() == 1 );
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = 1;
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
//cout << "getting compound CCD operator...\n";
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
//cout << "indices  " << ix << " " << jx << " " << kx << std::endl;

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

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

Npdm_op_wrapper_compound_DCD::Npdm_op_wrapper_compound_DCD( SpinBlock * spinBlock )
{
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
// 3-INDEX OPERATORS
//===========================================================================================================================================================

Npdm_op_wrapper_CCC::Npdm_op_wrapper_CCC( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_CRE_CRE).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_CRE_CRE).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
  // S={1/2,1/2,3/2}
  mults_ = { 2, 2, 4 };
  // For disk-based storage
  ifile_ = spinBlock_->get_op_array(CRE_CRE_CRE).get_filename();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CCC::set_local_ops( int idx )
{
//cout << "getting CCC operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core() )
    opReps_ = spinBlock_->get_op_array(CRE_CRE_CRE).get_local_element(idx);
  else {
    assert( check_file_open( idx ) );
    std::vector< boost::shared_ptr<SparseMatrix> > opReps_tmp;
    opReps_tmp = spinBlock_->get_op_array(CRE_CRE_CRE).get_local_element(idx);
    assert( opReps_tmp.at(0)->get_built_on_disk() );
    opReps_.clear();
    // Read in full spin-set from disk
    for (int i = 0; i < opReps_tmp.size(); i++) {
       boost::archive::binary_iarchive load_op(ifs_);
       boost::shared_ptr<SparseMatrix> op (new Cre);
       load_op >> *op;
       opReps_.push_back(op);
    }
    assert( check_file_close( idx ) );
  }

  build_pattern_ = opReps_.at(0)->get_build_pattern();
  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
//cout << "indices  " << ix << " " << jx << " " << kx << std::endl;
//cout << "build pattern " << opReps_.at(0)->get_build_pattern() << std::endl;
//cout << "2a CCC operator elements:\n";
//cout << *(opReps_[0]);
//cout << "2b CCC operator elements:\n";
//cout << *(opReps_[1]);
//cout << "4  CCC operator elements:\n";
//cout << *(opReps_[2]);
//  assert( build_pattern_ == opReps_.at(0)->get_build_pattern() );

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  if ( (ix == jx) && (jx == kx) ) {
    //FIXME This operator should not be built; it's zero as we cannot create 3 spin-1/2 particles with different spins
    //cout << "WARNING: skipping this operator\n";
    return true;
  }
  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_CCD::Npdm_op_wrapper_CCD( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_CRE_DES).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_CRE_DES).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
  // S={1/2,1/2,3/2}
  mults_ = { 2, 2, 4 };
  // For disk-based storage
  ifile_ = spinBlock_->get_op_array(CRE_CRE_DES).get_filename();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CCD::set_local_ops( int idx )
{
//cout << "getting CCD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core() )
    opReps_ = spinBlock_->get_op_array(CRE_CRE_DES).get_local_element(idx);
  else {
    assert( check_file_open( idx ) );
    std::vector< boost::shared_ptr<SparseMatrix> > opReps_tmp;
    opReps_tmp = spinBlock_->get_op_array(CRE_CRE_DES).get_local_element(idx);
    assert( opReps_tmp.at(0)->get_built_on_disk() );
    opReps_.clear();
    // Read in full spin-set from disk
    for (int i = 0; i < opReps_tmp.size(); i++) {
       boost::archive::binary_iarchive load_op(ifs_);
       boost::shared_ptr<SparseMatrix> op (new Cre);
       load_op >> *op;
       opReps_.push_back(op);
    }
    assert( check_file_close( idx ) );
  }

  build_pattern_ = opReps_.at(0)->get_build_pattern();
  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
//cout << "indices  " << ix << " " << jx << " " << kx << std::endl;
//cout << "build pattern " << opReps_.at(0)->get_build_pattern() << std::endl;
//cout << "2a CCD operator elements:\n";
//cout << *(opReps_[0]);
//cout << "2b CCD operator elements:\n";
//cout << *(opReps_[1]);
//cout << "4  CCD operator elements:\n";
//cout << *(opReps_[2]);
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
  is_local_ = spinBlock_->get_op_array(CRE_DES_DES).is_local();
  factor_ = 0.0;
  transpose_ = false;
  build_pattern_ = "0";
  // S={1/2,1/2,3/2}
  mults_ = { 2, 2, 4 };
  // For disk-based storage
  ifile_ = spinBlock_->get_op_array(CRE_DES_DES).get_filename();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CDD::set_local_ops( int idx )
{
//cout << "getting CDD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core() )
    opReps_ = spinBlock_->get_op_array(CRE_DES_DES).get_local_element(idx);
  else {
    assert( check_file_open( idx ) );
    std::vector< boost::shared_ptr<SparseMatrix> > opReps_tmp;
    opReps_tmp = spinBlock_->get_op_array(CRE_DES_DES).get_local_element(idx);
    assert( opReps_tmp.at(0)->get_built_on_disk() );
    opReps_.clear();
    // Read in full spin-set from disk
    for (int i = 0; i < opReps_tmp.size(); i++) {
       boost::archive::binary_iarchive load_op(ifs_);
       boost::shared_ptr<SparseMatrix> op (new Cre);
       load_op >> *op;
       opReps_.push_back(op);
    }
    assert( check_file_close( idx ) );
  }

  build_pattern_ = opReps_.at(0)->get_build_pattern();
  if ( build_pattern_ == "(C(DD))" ) factor_ = -1.0;
  else if ( build_pattern_ == "((CD)D)" ) factor_ = 1.0;
  else assert (false);

  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
//cout << "indices  " << ix << " " << jx << " " << kx << std::endl;
//cout << "build pattern " << opReps_.at(0)->get_build_pattern() << std::endl;
//cout << "2a CDD operator elements:\n";
//cout << *(opReps_[0]);
//cout << "2b CDD operator elements:\n";
//cout << *(opReps_[1]);
//cout << "4  CDD operator elements:\n";
//cout << *(opReps_[2]);
//  assert( build_pattern_ == opReps_.at(0)->get_build_pattern() );

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
//if ( (ix == jx) && (jx == kx) ) return false;
//cout << "WARNING: skipping this operator\n";
//return true;
  return false;

//FIXME GET AS TRANSPOSE of CCD and compare!!!
}

//===========================================================================================================================================================

Npdm_op_wrapper_CDC::Npdm_op_wrapper_CDC( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_DES_CRE).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_DES_CRE).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
  // S={1/2,1/2,3/2}
  mults_ = { 2, 2, 4 };
  // For disk-based storage
  ifile_ = spinBlock_->get_op_array(CRE_DES_CRE).get_filename();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CDC::set_local_ops( int idx )
{
//cout << "getting CDC operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core() )
    opReps_ = spinBlock_->get_op_array(CRE_DES_CRE).get_local_element(idx);
  else {
    assert( check_file_open( idx ) );
    std::vector< boost::shared_ptr<SparseMatrix> > opReps_tmp;
    opReps_tmp = spinBlock_->get_op_array(CRE_DES_CRE).get_local_element(idx);
    assert( opReps_tmp.at(0)->get_built_on_disk() );
    opReps_.clear();
    // Read in full spin-set from disk
    for (int i = 0; i < opReps_tmp.size(); i++) {
       boost::archive::binary_iarchive load_op(ifs_);
       boost::shared_ptr<SparseMatrix> op (new Cre);
       load_op >> *op;
       opReps_.push_back(op);
    }
    assert( check_file_close( idx ) );
  }

  build_pattern_ = opReps_.at(0)->get_build_pattern();
  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
//cout << "indices  " << ix << " " << jx << " " << kx << std::endl;

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  if ( jx == kx ) {
    //FIXME I think this fails because of potential problems commuting operators with same indices in spin-transformation
    //cout << "WARNING: skipping this operator\n";
    return true;
  }
  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_DCD::Npdm_op_wrapper_DCD( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(DES_CRE_DES).get_size();
  is_local_ = spinBlock_->get_op_array(DES_CRE_DES).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
  // S={1/2,1/2,3/2}
  mults_ = { 2, 2, 4 };
  // For disk-based storage
  ifile_ = spinBlock_->get_op_array(DES_CRE_DES).get_filename();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_DCD::set_local_ops( int idx )
{
//cout << "getting DCD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core() )
    opReps_ = spinBlock_->get_op_array(DES_CRE_DES).get_local_element(idx);
  else {
    assert( check_file_open( idx ) );
    std::vector< boost::shared_ptr<SparseMatrix> > opReps_tmp;
    opReps_tmp = spinBlock_->get_op_array(DES_CRE_DES).get_local_element(idx);
    assert( opReps_tmp.at(0)->get_built_on_disk() );
    opReps_.clear();
    // Read in full spin-set from disk
    for (int i = 0; i < opReps_tmp.size(); i++) {
       boost::archive::binary_iarchive load_op(ifs_);
       boost::shared_ptr<SparseMatrix> op (new Cre);
       load_op >> *op;
       opReps_.push_back(op);
    }
    assert( check_file_close( idx ) );
  }

  build_pattern_ = opReps_.at(0)->get_build_pattern();
  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
//cout << "indices  " << ix << " " << jx << " " << kx << std::endl;

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  if ( ix == jx ) {
    //FIXME I think this fails because of potential problems commuting operators with same indices in spin-transformation
    //cout << "WARNING: skipping this operator\n";
    return true;
  }
  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_DDC::Npdm_op_wrapper_DDC( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(DES_DES_CRE).get_size();
  is_local_ = spinBlock_->get_op_array(DES_DES_CRE).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
  // S={1/2,1/2,3/2}
  mults_ = { 2, 2, 4 };
  // For disk-based storage
  ifile_ = spinBlock_->get_op_array(DES_DES_CRE).get_filename();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_DDC::set_local_ops( int idx )
{
//cout << "getting DDC operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core() )
    opReps_ = spinBlock_->get_op_array(DES_DES_CRE).get_local_element(idx);
  else {
    assert( check_file_open( idx ) );
    std::vector< boost::shared_ptr<SparseMatrix> > opReps_tmp;
    opReps_tmp = spinBlock_->get_op_array(DES_DES_CRE).get_local_element(idx);
    assert( opReps_tmp.at(0)->get_built_on_disk() );
    opReps_.clear();
    // Read in full spin-set from disk
    for (int i = 0; i < opReps_tmp.size(); i++) {
       boost::archive::binary_iarchive load_op(ifs_);
       boost::shared_ptr<SparseMatrix> op (new Cre);
       load_op >> *op;
       opReps_.push_back(op);
    }
    assert( check_file_close( idx ) );
  }

  build_pattern_ = opReps_.at(0)->get_build_pattern();
  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
//cout << "indices  " << ix << " " << jx << " " << kx << std::endl;

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  if ( (jx == kx) || (ix == kx) ) {
    //FIXME I think this fails because of potential problems commuting operators with same indices in spin-transformation
    //cout << "WARNING: skipping this operator\n";
    return true;
  }
  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_DCC::Npdm_op_wrapper_DCC( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(DES_DES_CRE).get_size();
  is_local_ = spinBlock_->get_op_array(DES_DES_CRE).is_local();
//FIXME minus sign?
  factor_ = 1.0;
  transpose_ = true;
  build_pattern_ = "0";
  // S={1/2,1/2,3/2}
  mults_ = { 2, 2, 4 };
  // For disk-based storage
  ifile_ = spinBlock_->get_op_array(DES_DES_CRE).get_filename();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_DCC::set_local_ops( int idx )
{
//cout << "getting DCC operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core() )
    opReps_ = spinBlock_->get_op_array(DES_DES_CRE).get_local_element(idx);
  else {
    assert( check_file_open( idx ) );
    std::vector< boost::shared_ptr<SparseMatrix> > opReps_tmp;
    opReps_tmp = spinBlock_->get_op_array(DES_DES_CRE).get_local_element(idx);
    assert( opReps_tmp.at(0)->get_built_on_disk() );
    opReps_.clear();
    // Read in full spin-set from disk
    for (int i = 0; i < opReps_tmp.size(); i++) {
       boost::archive::binary_iarchive load_op(ifs_);
       boost::shared_ptr<SparseMatrix> op (new Cre);
       load_op >> *op;
       opReps_.push_back(op);
    }
    assert( check_file_close( idx ) );
  }

  std::string tmp = opReps_.at(0)->get_build_pattern();
  if ( tmp == "((DD)C)" ) build_pattern_ = "(D(CC))";
  else if ( tmp == "(D(DC))" ) build_pattern_ = "((DC)C)";
  else assert( false );

  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
//cout << "indices  " << ix << " " << jx << " " << kx << std::endl;

  // Note use of transpose means we store this as (k,j,i) not (i,j,k)
  indices_.push_back( kx );
  indices_.push_back( jx );
  indices_.push_back( ix );
  if ( (jx == kx) || (ix == kx) ) {
    //FIXME I think this fails because of potential problems commuting operators with same indices in spin-transformation
    //cout << "WARNING: skipping this operator\n";
    return true;
  }
  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_DDD::Npdm_op_wrapper_DDD( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_CRE_CRE).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_CRE_CRE).is_local();
//FIXME minus sign?
  factor_ = 1.0;
  transpose_ = true;
  build_pattern_ = "0";
  // S={1/2,1/2,3/2}
  mults_ = { 2, 2, 4 };
  // For disk-based storage
  ifile_ = spinBlock_->get_op_array(CRE_CRE_CRE).get_filename();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_DDD::set_local_ops( int idx )
{
//cout << "getting DDD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core() )
    opReps_ = spinBlock_->get_op_array(CRE_CRE_CRE).get_local_element(idx);
  else {
    assert( check_file_open( idx ) );
    std::vector< boost::shared_ptr<SparseMatrix> > opReps_tmp;
    opReps_tmp = spinBlock_->get_op_array(CRE_CRE_CRE).get_local_element(idx);
    assert( opReps_tmp.at(0)->get_built_on_disk() );
    opReps_.clear();
    // Read in full spin-set from disk
    for (int i = 0; i < opReps_tmp.size(); i++) {
       boost::archive::binary_iarchive load_op(ifs_);
       boost::shared_ptr<SparseMatrix> op (new Cre);
       load_op >> *op;
       opReps_.push_back(op);
    }
    assert( check_file_close( idx ) );
  }

  std::string tmp = opReps_.at(0)->get_build_pattern();
  if ( tmp == "((CC)C)" ) build_pattern_ = "(D(DD))";
  else if ( tmp == "(C(CC))" ) build_pattern_ = "((DD)D)";
  else assert( false );

  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
//cout << "indices  " << ix << " " << jx << " " << kx << std::endl;

  // Note use of transpose means we store this as (k,j,i) not (i,j,k)
  indices_.push_back( kx );
  indices_.push_back( jx );
  indices_.push_back( ix );
  if ( (ix == jx) && (jx == kx) ) {
    //FIXME This operator should not be built; it's zero as we cannot destroy 3 spin-1/2 particles with different spins
    return true;
  }
  return false;
}

//===========================================================================================================================================================

}
