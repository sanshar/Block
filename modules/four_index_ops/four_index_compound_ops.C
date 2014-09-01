/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "npdm_operators.h"
#include "npdm_spin_ops.h"
#include "pario.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

Npdm_op_compound_CCDD::Npdm_op_compound_CCDD( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(RI_4INDEX).get_size();
  is_local_ = true;
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((CC)(DD))";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_compound_CCDD::set_local_ops( int idx )
{
//pout << "getting RI CCDD operator...\n";
  // Spatial orbital indices
//FIXME don't need to keep reconstructing whole array!
  indices_ = spinBlock_->get_op_array(RI_4INDEX).get_array().at(idx);
  int ix = indices_[0];
  int jx = indices_[1];
  int kx = indices_[2];
  int lx = indices_[3];
  if(!dmrginp.spinAdapted()){
    // in nonspinAdapted situation, there are two sites in a dot block. Only one creator and one destructor on one sites. Creator is on the left of destructor.
    if(ix==jx) return true;
    if(kx==lx) return true;
    //because of ix>=jx>=kx>=lx, no C1C0D1D0. 
    //C1C0D1D0 is calculation in C1D1C0D0 in CDCD pattern
    //build_pattern_ = "((CD)(CD))";
    std::vector< boost::shared_ptr<SparseMatrix> > ijOps = spinBlock_->get_op_array(CRE_CRE).get_element(ix,jx);
    std::vector< boost::shared_ptr<SparseMatrix> > klOps = spinBlock_->get_op_array(DES_DES).get_element(kx,lx);

  opReps_.clear();
  //non spin, only at(0) is needed.
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(0), klOps.at(0), 0, indices_, false ) );
  return false;
  }
//pout << "indices: " << ix << ", " << jx << ", " << kx << ", " << lx << endl;

  // Get 2-index ops as RI building blocks
  std::vector< boost::shared_ptr<SparseMatrix> > ijOps = spinBlock_->get_op_array(CRE_CRE).get_element(ix,jx);
  std::vector< boost::shared_ptr<SparseMatrix> > klOps = spinBlock_->get_op_array(DES_DES).get_element(kx,lx);
//  std::vector< boost::shared_ptr<SparseMatrix> > klOps = spinBlock_->get_op_array(CRE_CRE).get_element(kx,lx);
//  // Take into account transpose
//  indices_[2] = lx;
//  indices_[3] = kx;
  
  // Allocate and build operator representation on the fly as RI tensor product for each spin component
  // multiplicites = { 1, 3, 3, 1, 3, 5 };
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

////
////  //-----------------------------
////  // Use RI 3-index operator
////  //-----------------------------
////  std::vector< boost::shared_ptr<SparseMatrix> > ijkOps;
////
////  // Get 2-index and 1-index ops as RI building blocks
////  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_CRE).get_element(ix,jx);
////  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(CRE).get_element(kx);
////
////  // Allocate and build operator representation on the fly as RI tensor product for each spin component
////  std::string build_pattern_3 = "((CC)D)";
////  std::vector<int> indices3 = { ix, jx, kx };
////  ijkOps.clear();
////  // S=0 (+) S=1/2  =>  S=1/2
////  ijkOps.push_back( build_compound_operator( true, -1, twoOps.at(0), oneOp.at(0), 0, indices3, true ) );
////  // S=1 (+) S=1/2  =>  S=1/2
////  ijkOps.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 0, indices3, true ) );
////  // S=1 (+) S=1/2  =>  S=3/2
////  ijkOps.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 1, indices3, true ) );
////
////
////  //-----------------------------
////  // Use exact 3-index operator
////  //-----------------------------
//////  std::vector< boost::shared_ptr<SparseMatrix> > ijkOps = spinBlock_->get_op_array(CRE_CRE_DES).get_element(ix,jx,kx);
//////  std::string build_pattern_3 = ijkOps.at(0)->get_build_pattern();
////
////  std::vector< boost::shared_ptr<SparseMatrix> > lOps = spinBlock_->get_op_array(CRE).get_element(lx);
////
////  if ( build_pattern_3 == "((CC)D)" ) build_pattern_ = "(((CC)D)(D))";
////  else if ( build_pattern_3 == "(C(CD))" ) build_pattern_ = "(((C)(CD))(D))";
////  else abort();
////  factor_ = 1.0;
////  multiplicities = { 1, 3, 1, 3, 3, 5 };
////
////  // Allocate and build operator representation on the fly as RI tensor product for each spin component
////  opReps_.clear();
////  // S=1/2 (+) S=1/2  =>  S=0
////  opReps_.push_back( build_compound_operator( false, -1, ijkOps.at(0), lOps.at(0), 0, indices_, true ) );
////  // S=1/2 (+) S=1/2  =>  S=1
////  opReps_.push_back( build_compound_operator( false, -1, ijkOps.at(0), lOps.at(0), 1, indices_, true ) );
////
////  // S=1/2 (+) S=1/2  =>  S=0
////  opReps_.push_back( build_compound_operator( false, -1, ijkOps.at(1), lOps.at(0), 0, indices_, true ) );
////  // S=1/2 (+) S=1/2  =>  S=1
////  opReps_.push_back( build_compound_operator( false, -1, ijkOps.at(1), lOps.at(0), 1, indices_, true ) );
////
////  // S=3/2 (+) S=1/2  =>  S=1
////  opReps_.push_back( build_compound_operator( false, -1, ijkOps.at(2), lOps.at(0), 0, indices_, true ) );
////  // S=3/2 (+) S=1/2  =>  S=2
////  opReps_.push_back( build_compound_operator( false, -1, ijkOps.at(2), lOps.at(0), 1, indices_, true ) );
////
////  return false;
}

//===========================================================================================================================================================

Npdm_op_compound_CCCD::Npdm_op_compound_CCCD( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(RI_4INDEX).get_size();
  is_local_ = true;
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((CC)(CD))";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_compound_CCCD::set_local_ops( int idx )
{
//pout << "getting RI CCCD operator...\n";
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
  // multiplicities = { 1, 3, 3, 1, 3, 5 };
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

Npdm_op_compound_CCDC::Npdm_op_compound_CCDC( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(RI_4INDEX).get_size();
  is_local_ = true;
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((CC)(DC))";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_compound_CCDC::set_local_ops( int idx )
{
//pout << "getting RI CCDC operator...\n";
  // Spatial orbital indices
//FIXME don't need to keep reconstructing whole array!
  indices_ = spinBlock_->get_op_array(RI_4INDEX).get_array().at(idx);
  int ix = indices_[0];
  int jx = indices_[1];
  int kx = indices_[2];
  int lx = indices_[3];
  if ( kx == lx ) {
     return true;
  }

  // Get 2-index ops as RI building blocks
  std::vector< boost::shared_ptr<SparseMatrix> > ijOps = spinBlock_->get_op_array(CRE_CRE).get_element(ix,jx);
  std::vector< boost::shared_ptr<SparseMatrix> > klOps = spinBlock_->get_op_array(DES_CRE).get_element(kx,lx);

  // Allocate and build operator representation on the fly as RI tensor product for each spin component
  // multiplicities = { 1, 3, 3, 1, 3, 5 };
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

  return false;
}

//===========================================================================================================================================================

Npdm_op_compound_CDCC::Npdm_op_compound_CDCC( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(RI_4INDEX).get_size();
  is_local_ = true;
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((CD)(CC))";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_compound_CDCC::set_local_ops( int idx )
{
//pout << "getting RI CDCC operator...\n";
  // Spatial orbital indices
//FIXME don't need to keep reconstructing whole array!
  indices_ = spinBlock_->get_op_array(RI_4INDEX).get_array().at(idx);
  int ix = indices_[0];
  int jx = indices_[1];
  int kx = indices_[2];
  int lx = indices_[3];
  if ( jx == kx ) {
     return true;
  }

  // Get 2-index ops as RI building blocks
  std::vector< boost::shared_ptr<SparseMatrix> > ijOps = spinBlock_->get_op_array(CRE_DES).get_element(ix,jx);
  std::vector< boost::shared_ptr<SparseMatrix> > klOps = spinBlock_->get_op_array(CRE_CRE).get_element(kx,lx);

  // Allocate and build operator representation on the fly as RI tensor product for each spin component
  // multiplicities = { 1, 3, 3, 1, 3, 5 };
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

Npdm_op_compound_CDCD::Npdm_op_compound_CDCD( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(RI_4INDEX).get_size();
  is_local_ = true;
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((CD)(CD))";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_compound_CDCD::set_local_ops( int idx )
{
//pout << "getting RI CDCD operator...\n";
  // Spatial orbital indices
//FIXME don't need to keep reconstructing whole array!
  indices_ = spinBlock_->get_op_array(RI_4INDEX).get_array().at(idx);
  int ix = indices_[0];
  int jx = indices_[1];
  int kx = indices_[2];
  int lx = indices_[3];
  if ( jx == kx ) {
     return true;
  }
  if(!dmrginp.spinAdapted()){
    // only one creator on a site
    if(ix==kx) return true;
    if(jx==lx) return true;
    // creator is before the destructor on a site
    if(jx==kx) return true;

  std::vector< boost::shared_ptr<SparseMatrix> > ijOps = spinBlock_->get_op_array(CRE_DES).get_element(ix,jx);
  std::vector< boost::shared_ptr<SparseMatrix> > klOps = spinBlock_->get_op_array(CRE_DES).get_element(kx,lx);

  opReps_.clear();
  opReps_.push_back( build_compound_operator( false, 1, ijOps.at(0), klOps.at(0), 0, indices_, false ) );
  return false;
  }

  // Get 2-index ops as RI building blocks
  std::vector< boost::shared_ptr<SparseMatrix> > ijOps = spinBlock_->get_op_array(CRE_DES).get_element(ix,jx);
  std::vector< boost::shared_ptr<SparseMatrix> > klOps = spinBlock_->get_op_array(CRE_DES).get_element(kx,lx);

  // Allocate and build operator representation on the fly as RI tensor product for each spin component
  // multiplicities = { 1, 3, 3, 1, 3, 5 };
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

Npdm_op_compound_CDDC::Npdm_op_compound_CDDC( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(RI_4INDEX).get_size();
  is_local_ = true;
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((CD)(DC))";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_compound_CDDC::set_local_ops( int idx )
{
//pout << "getting RI CDDC operator...\n";
  // Spatial orbital indices
//FIXME don't need to keep reconstructing whole array!
  indices_ = spinBlock_->get_op_array(RI_4INDEX).get_array().at(idx);
  int ix = indices_[0];
  int jx = indices_[1];
  int kx = indices_[2];
  int lx = indices_[3];
  if ( kx == lx ) {
     return true;
  }

  // Get 2-index ops as RI building blocks
  std::vector< boost::shared_ptr<SparseMatrix> > ijOps = spinBlock_->get_op_array(CRE_DES).get_element(ix,jx);
  std::vector< boost::shared_ptr<SparseMatrix> > klOps = spinBlock_->get_op_array(DES_CRE).get_element(kx,lx);

  // Allocate and build operator representation on the fly as RI tensor product for each spin component
  // multiplicities = { 1, 3, 3, 1, 3, 5 };
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

  return false;
}

//===========================================================================================================================================================

Npdm_op_compound_CDDD::Npdm_op_compound_CDDD( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(RI_4INDEX).get_size();
  is_local_ = true;
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((CD)(DD))";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_compound_CDDD::set_local_ops( int idx )
{
//pout << "getting RI CDDD operator...\n";
  // Spatial orbital indices
//FIXME don't need to keep reconstructing whole array!
  indices_ = spinBlock_->get_op_array(RI_4INDEX).get_array().at(idx);
  int ix = indices_[0];
  int jx = indices_[1];
  int kx = indices_[2];
  int lx = indices_[3];

  // Get 2-index ops as RI building blocks
  std::vector< boost::shared_ptr<SparseMatrix> > ijOps = spinBlock_->get_op_array(CRE_DES).get_element(ix,jx);
  std::vector< boost::shared_ptr<SparseMatrix> > klOps = spinBlock_->get_op_array(DES_DES).get_element(kx,lx);

  // Allocate and build operator representation on the fly as RI tensor product for each spin component
  // multiplicities = { 1, 3, 3, 1, 3, 5 };
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

Npdm_op_compound_CCCC::Npdm_op_compound_CCCC( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(RI_4INDEX).get_size();
  is_local_ = true;
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((CC)(CC))";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_compound_CCCC::set_local_ops( int idx )
{
//pout << "getting RI CCCC operator...\n";
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
  // multiplicities = { 1, 3, 3, 1, 3, 5 };
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
}
