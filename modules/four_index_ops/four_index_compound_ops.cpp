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

//////===========================================================================================================================================================
////
////Npdm_op_wrapper_compound_CCDD::Npdm_op_wrapper_compound_CCDD( SpinBlock * spinBlock )
////{
////  // Assume single site block
////  assert( spinBlock->size() == 1 );
////  opReps_.clear();
////  indices_.clear();
////  spinBlock_ = spinBlock;
////  size_ = 1;
////  is_local_ = true;
//////FIXME why do we need -1 here -- similar to (C(DD)) case ? (use of transpose?)
////  factor_ = -1.0;
////  transpose_ = false;
////  build_pattern_ = "((CC)(DD))";
////  mults_ = { 1, 3, 3, 1, 3, 5 };
////}
////
//////-----------------------------------------------------------------------------------------------------------------------------------------------------------
////
////bool Npdm_op_wrapper_compound_CCDD::set_local_ops( int idx )
////{
//////cout << "getting CCDD operator...\n";
////  // Spatial orbital indices
////  indices_.clear();
////  int ix, jx, kx, lx;
////  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_CRE).get_local_element(idx);
////  ix = twoOps.at(0)->get_orbs(0);
////  jx = twoOps.at(0)->get_orbs(1);
//////cout << "indices  " << ix << " " << ix << " " << ix << " " << ix << std::endl;
////
////  // Assumed single site
////  assert ( ix == jx );
////  indices_.push_back( ix );
////  indices_.push_back( ix );
////  indices_.push_back( ix );
////  indices_.push_back( ix );
////
////  opReps_.clear();
////  // S=0 (+) S=0  =>  S=0
////  opReps_.push_back( build_compound_operator( false, -1, twoOps.at(0), twoOps.at(0), 0, indices_, true ) );
////  // S=1 (+) S=0  =>  S=1
////  opReps_.push_back( build_compound_operator( false, -1, twoOps.at(1), twoOps.at(0), 0, indices_, true ) );
////  // S=0 (+) S=1  =>  S=1
////  opReps_.push_back( build_compound_operator( false, -1, twoOps.at(0), twoOps.at(1), 0, indices_, true ) );
////
////  // S=1 (+) S=1  =>  S=0
////  opReps_.push_back( build_compound_operator( false, -1, twoOps.at(1), twoOps.at(1), 0, indices_, true ) );
////  // S=1 (+) S=1  =>  S=1
////  opReps_.push_back( build_compound_operator( false, -1, twoOps.at(1), twoOps.at(1), 1, indices_, true ) );
////  // S=1 (+) S=1  =>  S=2
////  opReps_.push_back( build_compound_operator( false, -1, twoOps.at(1), twoOps.at(1), 2, indices_, true ) );
////
////  return false;
////}
////
//===========================================================================================================================================================

Npdm_op_wrapper_compound_CCDD::Npdm_op_wrapper_compound_CCDD( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(RI_4INDEX).get_size();
  is_local_ = true;
//FIXME note minus sign
  factor_ = -1.0;
  transpose_ = false;
  build_pattern_ = "((CC)(DD))";
  mults_ = { 1, 3, 3, 1, 3, 5 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_compound_CCDD::set_local_ops( int idx )
{
cout << "getting RI CCDD operator...\n";
  // Spatial orbital indices
//FIXME don't need to keep reconstructing whole array!
  indices_ = spinBlock_->get_op_array(RI_4INDEX).get_array().at(idx);
  int ix = indices_[0];
  int jx = indices_[1];
  int kx = indices_[2];
  int lx = indices_[3];

  // Get 2-index ops as RI building blocks
  std::vector< boost::shared_ptr<SparseMatrix> > ijOps = spinBlock_->get_op_array(CRE_CRE).get_element(ix,jx);
  // Need (lx,kx) after transpose, but not available, so introduce minus sign for commutation
  std::vector< boost::shared_ptr<SparseMatrix> > klOps = spinBlock_->get_op_array(CRE_CRE).get_element(kx,lx);

  // Allocate and build operator representation on the fly as RI tensor product for each spin component
  opReps_.clear();
  // S=0 (+) S=0  =>  S=0
  opReps_.push_back( build_compound_operator( false, -1, ijOps.at(0), klOps.at(0), 0, indices_, true ) );
  // S=1 (+) S=0  =>  S=1
  opReps_.push_back( build_compound_operator( false, -1, ijOps.at(1), klOps.at(0), 0, indices_, true ) );
  // S=0 (+) S=1  =>  S=1
  opReps_.push_back( build_compound_operator( false, -1, ijOps.at(0), klOps.at(1), 0, indices_, true ) );

  // S=1 (+) S=1  =>  S=0
  opReps_.push_back( build_compound_operator( false, -1, ijOps.at(1), klOps.at(1), 0, indices_, true ) );
  // S=1 (+) S=1  =>  S=1
  opReps_.push_back( build_compound_operator( false, -1, ijOps.at(1), klOps.at(1), 1, indices_, true ) );
  // S=1 (+) S=1  =>  S=2
  opReps_.push_back( build_compound_operator( false, -1, ijOps.at(1), klOps.at(1), 2, indices_, true ) );

  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_compound_CCCD::Npdm_op_wrapper_compound_CCCD( SpinBlock * spinBlock )
{
  // Assume single site block
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(RI_4INDEX).get_size();
  is_local_ = true;
//FIXME sign?
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((CC)(CD))";
  mults_ = { 1, 3, 3, 1, 3, 5 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_compound_CCCD::set_local_ops( int idx )
{
cout << "getting RI CCCD operator...\n";
  // Spatial orbital indices
//FIXME don't need to keep reconstructing whole array!
  indices_ = spinBlock_->get_op_array(RI_4INDEX).get_array().at(idx);
  int ix = indices_[0];
  int jx = indices_[1];
  int kx = indices_[2];
  int lx = indices_[3];

  // Get 2-index ops as RI building blocks
  std::vector< boost::shared_ptr<SparseMatrix> > ijOps = spinBlock_->get_op_array(CRE_CRE).get_element(ix,jx);
  std::vector< boost::shared_ptr<SparseMatrix> > klOps = spinBlock_->get_op_array(CRE_DES).get_element(kx,lx);

  // Allocate and build operator representation on the fly as RI tensor product for each spin component
  opReps_.clear();
  // S=0 (+) S=0  =>  S=0
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(0), klOps.at(0), 0, indices_, false ) );
  // S=1 (+) S=0  =>  S=1
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(1), klOps.at(0), 0, indices_, false ) );
  // S=0 (+) S=1  =>  S=1
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(0), klOps.at(1), 0, indices_, false ) );

  // S=1 (+) S=1  =>  S=0
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(1), klOps.at(1), 0, indices_, false ) );
  // S=1 (+) S=1  =>  S=1
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(1), klOps.at(1), 1, indices_, false ) );
  // S=1 (+) S=1  =>  S=2
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(1), klOps.at(1), 2, indices_, false ) );

  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_compound_CCDC::Npdm_op_wrapper_compound_CCDC( SpinBlock * spinBlock )
{
  // Assume single site block
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(RI_4INDEX).get_size();
  is_local_ = true;
//FIXME sign?
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((CC)(DC))";
  mults_ = { 1, 3, 3, 1, 3, 5 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_compound_CCDC::set_local_ops( int idx )
{
cout << "getting RI CCDC operator...\n";
  // Spatial orbital indices
//FIXME don't need to keep reconstructing whole array!
  indices_ = spinBlock_->get_op_array(RI_4INDEX).get_array().at(idx);
  int ix = indices_[0];
  int jx = indices_[1];
  int kx = indices_[2];
  int lx = indices_[3];

  // Get 2-index ops as RI building blocks
  std::vector< boost::shared_ptr<SparseMatrix> > ijOps = spinBlock_->get_op_array(CRE_CRE).get_element(ix,jx);
  std::vector< boost::shared_ptr<SparseMatrix> > klOps = spinBlock_->get_op_array(DES_CRE).get_element(kx,lx);

  // Allocate and build operator representation on the fly as RI tensor product for each spin component
  opReps_.clear();
  // S=0 (+) S=0  =>  S=0
  opReps_.push_back( build_compound_operator( false, -1, ijOps.at(0), klOps.at(0), 0, indices_, false ) );
  // S=1 (+) S=0  =>  S=1
  opReps_.push_back( build_compound_operator( false, -1, ijOps.at(1), klOps.at(0), 0, indices_, false ) );
  // S=0 (+) S=1  =>  S=1
  opReps_.push_back( build_compound_operator( false, -1, ijOps.at(0), klOps.at(1), 0, indices_, false ) );

  // S=1 (+) S=1  =>  S=0
  opReps_.push_back( build_compound_operator( false, -1, ijOps.at(1), klOps.at(1), 0, indices_, false ) );
  // S=1 (+) S=1  =>  S=1
  opReps_.push_back( build_compound_operator( false, -1, ijOps.at(1), klOps.at(1), 1, indices_, false ) );
  // S=1 (+) S=1  =>  S=2
  opReps_.push_back( build_compound_operator( false, -1, ijOps.at(1), klOps.at(1), 2, indices_, false ) );

  if ( kx == lx ) {
     return true;
  }
  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_compound_CDCC::Npdm_op_wrapper_compound_CDCC( SpinBlock * spinBlock )
{
  // Assume single site block
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(RI_4INDEX).get_size();
  is_local_ = true;
//FIXME sign?
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((CD)(CC))";
  mults_ = { 1, 3, 3, 1, 3, 5 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_compound_CDCC::set_local_ops( int idx )
{
cout << "getting RI CDCC operator...\n";
  // Spatial orbital indices
//FIXME don't need to keep reconstructing whole array!
  indices_ = spinBlock_->get_op_array(RI_4INDEX).get_array().at(idx);
  int ix = indices_[0];
  int jx = indices_[1];
  int kx = indices_[2];
  int lx = indices_[3];

  // Get 2-index ops as RI building blocks
  std::vector< boost::shared_ptr<SparseMatrix> > ijOps = spinBlock_->get_op_array(CRE_DES).get_element(ix,jx);
  std::vector< boost::shared_ptr<SparseMatrix> > klOps = spinBlock_->get_op_array(CRE_CRE).get_element(kx,lx);

  // Allocate and build operator representation on the fly as RI tensor product for each spin component
  opReps_.clear();
  // S=0 (+) S=0  =>  S=0
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(0), klOps.at(0), 0, indices_, false ) );
  // S=1 (+) S=0  =>  S=1
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(1), klOps.at(0), 0, indices_, false ) );
  // S=0 (+) S=1  =>  S=1
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(0), klOps.at(1), 0, indices_, false ) );

  // S=1 (+) S=1  =>  S=0
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(1), klOps.at(1), 0, indices_, false ) );
  // S=1 (+) S=1  =>  S=1
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(1), klOps.at(1), 1, indices_, false ) );
  // S=1 (+) S=1  =>  S=2
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(1), klOps.at(1), 2, indices_, false ) );

  if ( jx == kx ) {
     return true;
  }
  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_compound_CDCD::Npdm_op_wrapper_compound_CDCD( SpinBlock * spinBlock )
{
  // Assume single site block
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(RI_4INDEX).get_size();
  is_local_ = true;
//FIXME sign?
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((CD)(CD))";
  mults_ = { 1, 3, 3, 1, 3, 5 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_compound_CDCD::set_local_ops( int idx )
{
cout << "getting RI CDCD operator...\n";
  // Spatial orbital indices
//FIXME don't need to keep reconstructing whole array!
  indices_ = spinBlock_->get_op_array(RI_4INDEX).get_array().at(idx);
  int ix = indices_[0];
  int jx = indices_[1];
  int kx = indices_[2];
  int lx = indices_[3];

  // Get 2-index ops as RI building blocks
  std::vector< boost::shared_ptr<SparseMatrix> > ijOps = spinBlock_->get_op_array(CRE_DES).get_element(ix,jx);
  std::vector< boost::shared_ptr<SparseMatrix> > klOps = spinBlock_->get_op_array(CRE_DES).get_element(kx,lx);

  // Allocate and build operator representation on the fly as RI tensor product for each spin component
  opReps_.clear();
  // S=0 (+) S=0  =>  S=0
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(0), klOps.at(0), 0, indices_, false ) );
  // S=1 (+) S=0  =>  S=1
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(1), klOps.at(0), 0, indices_, false ) );
  // S=0 (+) S=1  =>  S=1
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(0), klOps.at(1), 0, indices_, false ) );

  // S=1 (+) S=1  =>  S=0
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(1), klOps.at(1), 0, indices_, false ) );
  // S=1 (+) S=1  =>  S=1
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(1), klOps.at(1), 1, indices_, false ) );
  // S=1 (+) S=1  =>  S=2
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(1), klOps.at(1), 2, indices_, false ) );

  if ( jx == kx ) {
     return true;
  }
  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_compound_CDDC::Npdm_op_wrapper_compound_CDDC( SpinBlock * spinBlock )
{
  // Assume single site block
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(RI_4INDEX).get_size();
  is_local_ = true;
//FIXME sign?
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((CD)(DC))";
  mults_ = { 1, 3, 3, 1, 3, 5 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_compound_CDDC::set_local_ops( int idx )
{
cout << "getting RI CDDC operator...\n";
  // Spatial orbital indices
//FIXME don't need to keep reconstructing whole array!
  indices_ = spinBlock_->get_op_array(RI_4INDEX).get_array().at(idx);
  int ix = indices_[0];
  int jx = indices_[1];
  int kx = indices_[2];
  int lx = indices_[3];

  // Get 2-index ops as RI building blocks
  std::vector< boost::shared_ptr<SparseMatrix> > ijOps = spinBlock_->get_op_array(CRE_DES).get_element(ix,jx);
  std::vector< boost::shared_ptr<SparseMatrix> > klOps = spinBlock_->get_op_array(DES_CRE).get_element(kx,lx);

  // Allocate and build operator representation on the fly as RI tensor product for each spin component
  opReps_.clear();
  // S=0 (+) S=0  =>  S=0
  opReps_.push_back( build_compound_operator( false, -1, ijOps.at(0), klOps.at(0), 0, indices_, false ) );
  // S=1 (+) S=0  =>  S=1
  opReps_.push_back( build_compound_operator( false, -1, ijOps.at(1), klOps.at(0), 0, indices_, false ) );
  // S=0 (+) S=1  =>  S=1
  opReps_.push_back( build_compound_operator( false, -1, ijOps.at(0), klOps.at(1), 0, indices_, false ) );

  // S=1 (+) S=1  =>  S=0
  opReps_.push_back( build_compound_operator( false, -1, ijOps.at(1), klOps.at(1), 0, indices_, false ) );
  // S=1 (+) S=1  =>  S=1
  opReps_.push_back( build_compound_operator( false, -1, ijOps.at(1), klOps.at(1), 1, indices_, false ) );
  // S=1 (+) S=1  =>  S=2
  opReps_.push_back( build_compound_operator( false, -1, ijOps.at(1), klOps.at(1), 2, indices_, false ) );

  if ( kx == lx ) {
     return true;
  }
  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_compound_CDDD::Npdm_op_wrapper_compound_CDDD( SpinBlock * spinBlock )
{
  // Assume single site block
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(RI_4INDEX).get_size();
  is_local_ = true;
  factor_ = 0.0;
  transpose_ = false;
  build_pattern_ = "";
  mults_ = { };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_compound_CDDD::set_local_ops( int idx )
{
cout << "getting RI CDDD operator...\n";
  // Spatial orbital indices
//FIXME don't need to keep reconstructing whole array!
  indices_ = spinBlock_->get_op_array(RI_4INDEX).get_array().at(idx);
  int ix = indices_[0];
  int jx = indices_[1];
  int kx = indices_[2];
  int lx = indices_[3];

//FIXME  THIS WAY IS BROKEN!!!!
//----------------------------------
////  factor_ = -1.0;
////  transpose_ = false;
////  build_pattern_ = "((CD)(DD))";
////  mults_ = { 1, 3, 3, 1, 3, 5 };
////  // Get 2-index ops as RI building blocks
////  std::vector< boost::shared_ptr<SparseMatrix> > ijOps = spinBlock_->get_op_array(CRE_DES).get_element(ix,jx);
////  // Need (lx,kx) after transpose, but not available, so introduce minus sign for commutation
////  std::vector< boost::shared_ptr<SparseMatrix> > klOps = spinBlock_->get_op_array(CRE_CRE).get_element(kx,lx);
////
////  // Allocate and build operator representation on the fly as RI tensor product for each spin component
////  opReps_.clear();
////  // S=0 (+) S=0  =>  S=0
////  opReps_.push_back( build_compound_operator( false, -1, ijOps.at(0), klOps.at(0), 0, indices_, true ) );
////  // S=1 (+) S=0  =>  S=1
////  opReps_.push_back( build_compound_operator( false, -1, ijOps.at(1), klOps.at(0), 0, indices_, true ) );
////  // S=0 (+) S=1  =>  S=1
////  opReps_.push_back( build_compound_operator( false, -1, ijOps.at(0), klOps.at(1), 0, indices_, true ) );
////
////  // S=1 (+) S=1  =>  S=0
////  opReps_.push_back( build_compound_operator( false, -1, ijOps.at(1), klOps.at(1), 0, indices_, true ) );
////  // S=1 (+) S=1  =>  S=1
////  opReps_.push_back( build_compound_operator( false, -1, ijOps.at(1), klOps.at(1), 1, indices_, true ) );
////  // S=1 (+) S=1  =>  S=2
////  opReps_.push_back( build_compound_operator( false, -1, ijOps.at(1), klOps.at(1), 2, indices_, true ) );

//------------------------
//FIXME
// 2nd way seems to work
//------------------------

  factor_ = 1.0;
  build_pattern_ = "((CDD)(D))";
  mults_ = { 1, 3, 1, 3, 3, 5 };
  std::vector< boost::shared_ptr<SparseMatrix> > ijkOps = spinBlock_->get_op_array(CRE_DES_DES).get_element(ix,jx,kx);
  std::vector< boost::shared_ptr<SparseMatrix> > lOps = spinBlock_->get_op_array(CRE).get_element(lx);

  // Allocate and build operator representation on the fly as RI tensor product for each spin component
  opReps_.clear();
  // S=1/2 (+) S=1/2  =>  S=0
  opReps_.push_back( build_compound_operator( false, -1, ijkOps.at(0), lOps.at(0), 0, indices_, true ) );
  // S=1/2 (+) S=1/2  =>  S=1
  opReps_.push_back( build_compound_operator( false, -1, ijkOps.at(0), lOps.at(0), 1, indices_, true ) );

  // S=1/2 (+) S=1/2  =>  S=0
  opReps_.push_back( build_compound_operator( false, -1, ijkOps.at(1), lOps.at(0), 0, indices_, true ) );
  // S=1/2 (+) S=1/2  =>  S=1
  opReps_.push_back( build_compound_operator( false, -1, ijkOps.at(1), lOps.at(0), 1, indices_, true ) );

  // S=3/2 (+) S=1/2  =>  S=1
  opReps_.push_back( build_compound_operator( false, -1, ijkOps.at(2), lOps.at(0), 0, indices_, true ) );
  // S=3/2 (+) S=1/2  =>  S=2
  opReps_.push_back( build_compound_operator( false, -1, ijkOps.at(2), lOps.at(0), 1, indices_, true ) );

  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_compound_CCCC::Npdm_op_wrapper_compound_CCCC( SpinBlock * spinBlock )
{
  // Assume single site block
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(RI_4INDEX).get_size();
  is_local_ = true;
//FIXME sign?
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((CC)(CC))";
  mults_ = { 1, 3, 3, 1, 3, 5 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_compound_CCCC::set_local_ops( int idx )
{
cout << "getting RI CCCC operator...\n";
  // Spatial orbital indices
//FIXME don't need to keep reconstructing whole array!
  indices_ = spinBlock_->get_op_array(RI_4INDEX).get_array().at(idx);
  int ix = indices_[0];
  int jx = indices_[1];
  int kx = indices_[2];
  int lx = indices_[3];

  // Get 2-index ops as RI building blocks
  std::vector< boost::shared_ptr<SparseMatrix> > ijOps = spinBlock_->get_op_array(CRE_CRE).get_element(ix,jx);
  std::vector< boost::shared_ptr<SparseMatrix> > klOps = spinBlock_->get_op_array(CRE_CRE).get_element(kx,lx);

  // Allocate and build operator representation on the fly as RI tensor product for each spin component
  opReps_.clear();
  // S=0 (+) S=0  =>  S=0
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(0), klOps.at(0), 0, indices_, false ) );
  // S=1 (+) S=0  =>  S=1
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(1), klOps.at(0), 0, indices_, false ) );
  // S=0 (+) S=1  =>  S=1
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(0), klOps.at(1), 0, indices_, false ) );

  // S=1 (+) S=1  =>  S=0
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(1), klOps.at(1), 0, indices_, false ) );
  // S=1 (+) S=1  =>  S=1
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(1), klOps.at(1), 1, indices_, false ) );
  // S=1 (+) S=1  =>  S=2
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(1), klOps.at(1), 2, indices_, false ) );

  return false;
}

//===========================================================================================================================================================

}
