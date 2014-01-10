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
// 4-INDEX COMPOUND OPERATORS (for one site only)
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
// 4-INDEX OPERATORS
//===========================================================================================================================================================

Npdm_op_wrapper_CCDD::Npdm_op_wrapper_CCDD( SpinBlock * spinBlock )
{
  assert(false);
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_CRE_DES_DES).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_CRE_DES_DES).is_local();
//FIXME why do we need -1 here -- similar to (C(DD)) case ? (use of transpose?)
  factor_ = -1.0;
  transpose_ = false;
//  build_pattern_ = "((CC)(DD))";
  build_pattern_ = "0";
  mults_ = { 1, 3, 3, 1, 3, 5 };
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CCDD::set_local_ops( int idx )
{
//cout << "getting CCDD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx, lx;

  opReps_ = spinBlock_->get_op_array(CRE_CRE_DES_DES).get_local_element(idx);
  build_pattern_ = opReps_.at(0)->get_build_pattern();
  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
  lx = opReps_.at(0)->get_orbs(3);

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  indices_.push_back( lx );

  return false;
}

//===========================================================================================================================================================

}
