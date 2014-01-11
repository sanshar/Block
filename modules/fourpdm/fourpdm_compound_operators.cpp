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

Npdm_op_wrapper_compound_CCDD::Npdm_op_wrapper_compound_CCDD( SpinBlock * spinBlock )
{
  // Assume single site block
  assert( spinBlock->size() == 1 );
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = 1;
  is_local_ = true;
//FIXME why do we need -1 here -- similar to (C(DD)) case ? (use of transpose?)
  factor_ = -1.0;
  transpose_ = false;
  build_pattern_ = "((CC)(DD))";
  mults_ = { 1, 3, 3, 1, 3, 5 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_compound_CCDD::set_local_ops( int idx )
{
//cout << "getting CCDD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx, lx;
  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_CRE).get_local_element(idx);
  ix = twoOps.at(0)->get_orbs(0);
  jx = twoOps.at(0)->get_orbs(1);
//cout << "indices  " << ix << " " << ix << " " << ix << " " << ix << std::endl;

  // Assumed single site
  assert ( ix == jx );
  indices_.push_back( ix );
  indices_.push_back( ix );
  indices_.push_back( ix );
  indices_.push_back( ix );

  opReps_.clear();
  // S=0 (+) S=0  =>  S=0
  opReps_.push_back( build_compound_operator( false, -1, twoOps.at(0), twoOps.at(0), 0, indices_, true ) );
  // S=1 (+) S=0  =>  S=1
  opReps_.push_back( build_compound_operator( false, -1, twoOps.at(1), twoOps.at(0), 0, indices_, true ) );
  // S=0 (+) S=1  =>  S=1
  opReps_.push_back( build_compound_operator( false, -1, twoOps.at(0), twoOps.at(1), 0, indices_, true ) );

  // S=1 (+) S=1  =>  S=0
  opReps_.push_back( build_compound_operator( false, -1, twoOps.at(1), twoOps.at(1), 0, indices_, true ) );
  // S=1 (+) S=1  =>  S=1
  opReps_.push_back( build_compound_operator( false, -1, twoOps.at(1), twoOps.at(1), 1, indices_, true ) );
  // S=1 (+) S=1  =>  S=2
  opReps_.push_back( build_compound_operator( false, -1, twoOps.at(1), twoOps.at(1), 2, indices_, true ) );

  return false;
}

//===========================================================================================================================================================

}

//===========================================================================================================================================================
//
//Npdm_op_wrapper_compound_CCD::Npdm_op_wrapper_compound_CCD( SpinBlock * spinBlock )
//{
//  opReps_.clear();
//  indices_.clear();
//  spinBlock_ = spinBlock;
//  size_ = spinBlock_->get_op_array(RI_3INDEX).get_size();
//  is_local_ = true;
//  factor_ = 1.0;
//  transpose_ = false;
//  build_pattern_ = "((CC)D)";
//  // S={1/2,1/2,3/2}
//  mults_ = { 2, 2, 4 };
//}
//
////-----------------------------------------------------------------------------------------------------------------------------------------------------------
//
//bool Npdm_op_wrapper_compound_CCD::set_local_ops( int idx )
//{
//cout << "getting compound CCD operator...\n";
//  // Spatial orbital indices
////FIXME don't need to keep reconstructing whole array!
//  indices_ = spinBlock_->get_op_array(RI_3INDEX).get_array().at(idx);
//  int ix = indices_[0];
//  int jx = indices_[1];
//  int kx = indices_[2];
//
//  // Get 2-index and 1-index ops as RI building blocks
//  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_CRE).get_element(ix,jx);
//  assert( ix == twoOps.at(0)->get_orbs(0) );
//  assert( jx == twoOps.at(0)->get_orbs(1) );
//  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(CRE).get_element(kx);
//  assert( kx == oneOp.at(0)->get_orbs(0) );
//
//  // Allocate and build operator representation on the fly as RI tensor product for each spin component
//  opReps_.clear();
//  // S=0 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(0), oneOp.at(0), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=3/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 1, indices_, true ) );
//  return false;
//}
//
//===========================================================================================================================================================

