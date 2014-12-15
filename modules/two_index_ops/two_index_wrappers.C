/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "global.h"
#include "npdm_patterns.h"
#include "npdm_operators.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================
// 3-PDM OPERATORS
//===========================================================================================================================================================

Npdm_op_wrapper_DC::Npdm_op_wrapper_DC( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(DES_CRE).get_size();
  is_local_ = spinBlock_->get_op_array(DES_CRE).is_local();
  factor_ = 1.0;
  //Definition of D_iC_j is C_jD_i.
  //In npdm code, it is used as D_iC_j.
  //For spin adapted npdm, this difference can be detected by the solver of linear equation of spin coupling. 
  //However, for non spin adapted npdm, there is only a simply way to determine the sign before elements. It can not be detected. 
  //Therefore, -1 is necessary. 
  if(!dmrginp.spinAdapted()) factor_*=-1;
  transpose_ = false;
  build_pattern_ = "(DC)";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_DC::set_local_ops( int idx )
{
//pout << "getting DC operator...\n";

  // Spatial orbital indices
  indices_.clear();
  int ix, jx;

  opReps_ = spinBlock_->get_op_array(DES_CRE).get_local_element(idx);
  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
//pout << "singlet DC operator elements:\n";
//pout << *(opReps_[0]);
//pout << "triplet DC operator elements:\n";
//pout << *(opReps_[1]);

  indices_.push_back( ix );
  indices_.push_back( jx );
//pout << "indices  " << ix << " " << jx << std::endl;
if ( ix == jx ) {
  //FIXME I think this fails because of potential problems commuting operators with same indices in spin-transformation
  //pout << "WARNING: skipping this operator\n";
  return true;
}
  return false;

}

////// DEBUG build as (RI) product  --->  seems to give same result (within RI approx...)
////
////  indices_.clear();
////  int i, j;
////  opReps_ = spinBlock_->get_op_array(DES_CRE).get_local_element(idx);
////  i = opReps_.at(0)->get_orbs(0);
////  j = opReps_.at(0)->get_orbs(1);
////  indices_.push_back( i );
////  indices_.push_back( j );
////
////  const boost::shared_ptr<SparseMatrix>& opCi = spinBlock_->get_op_rep(CRE, getSpinQuantum(i), i);
////  const boost::shared_ptr<SparseMatrix>& opCj = spinBlock_->get_op_rep(CRE, getSpinQuantum(j), j);
//// 
////  opReps_.clear();
////  // singlet
////  boost::shared_ptr<SparseMatrix> singOp (new DesCre);
////  singOp->set_orbs() = indices_;
////  singOp->set_initialised() = true;
////  singOp->set_fermion() = false;;
//////NOTE TRANSPOSEVIEW breaks get_deltaQuantum??
////  singOp->set_deltaQuantum() = ( - opCi->get_deltaQuantum() + opCj->get_deltaQuantum() ).at(0);
////  singOp->allocate( spinBlock_->get_stateInfo() );
////  operatorfunctions::Product(spinBlock_, Transposeview(*opCi), *opCj, *singOp, 1.0 );
////  opReps_.push_back( singOp );
////
////  // triplet
////  boost::shared_ptr<SparseMatrix> tripOp (new DesCre);
////  tripOp->set_orbs() = indices_;
////  tripOp->set_initialised() = true;
////  tripOp->set_fermion() = false;;
////  tripOp->set_deltaQuantum() = ( - opCi->get_deltaQuantum() + opCj->get_deltaQuantum() ).at(1);
////  tripOp->allocate( spinBlock_->get_stateInfo() );
////  operatorfunctions::Product(spinBlock_, Transposeview(*opCi), *opCj, *tripOp, 1.0 );
////  opReps_.push_back( tripOp );
////
////if ( i == j ) {
////  //FIXME I think this fails because of potential problems commuting operators with same indices in spin-transformation
////  pout << "WARNING: skipping this operator\n";
////  return true;
////}
////  return false;
////
////}

//===========================================================================================================================================================
// 2-PDM OPERATORS
//===========================================================================================================================================================

Npdm_op_wrapper_CC::Npdm_op_wrapper_CC( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_CRE).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_CRE).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "(CC)";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CC::set_local_ops( int idx )
{
//pout << "getting CC operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx;

  opReps_ = spinBlock_->get_op_array(CRE_CRE).get_local_element(idx);
  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
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
  is_local_ = spinBlock_->get_op_array(CRE_DES).is_local();
  factor_ = 1.0;
  build_pattern_ = "(CD)";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CD::set_local_ops( int idx )
{
  // Spatial orbital indices
  indices_.clear();
  int ix, jx;

  opReps_ = spinBlock_->get_op_array(CRE_DES).get_local_element(idx);
  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  transpose_ = false;
  indices_.push_back( ix );
  indices_.push_back( jx );

  return false;
}

//===========================================================================================================================================================
// Build DD as transpose(CC)
//
//Npdm_op_wrapper_DD::Npdm_op_wrapper_DD( SpinBlock * spinBlock )
//{
//  opReps_.clear();
//  indices_.clear();
//  spinBlock_ = spinBlock;
//  size_ = spinBlock_->get_op_array(CRE_CRE).get_size();
//  is_local_ = spinBlock_->get_op_array(CRE_CRE).is_local();
////FIXME why do we need -1 here ??
//  factor_ = -1.0;
//  transpose_ = true;
//  build_pattern_ = "(DD)";
//}
//
////-----------------------------------------------------------------------------------------------------------------------------------------------------------
//
//bool Npdm_op_wrapper_DD::set_local_ops( int idx )
//{
//  // Spatial orbital indices
//  indices_.clear();
//  int ix, jx;
//
//  opReps_ = spinBlock_->get_op_array(CRE_CRE).get_local_element(idx);
//  ix = opReps_.at(0)->get_orbs(0);
//  jx = opReps_.at(0)->get_orbs(1);
//
//  // Note use of transpose means we store this as (j,i) not (i,j)
//  indices_.push_back( jx );
//  indices_.push_back( ix );
//  return false;
//}
//
//===========================================================================================================================================================

Npdm_op_wrapper_DD::Npdm_op_wrapper_DD( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(DES_DES).get_size();
  is_local_ = spinBlock_->get_op_array(DES_DES).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "(DD)";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_DD::set_local_ops( int idx )
{
  // Spatial orbital indices
  indices_.clear();
  int ix, jx;

  opReps_ = spinBlock_->get_op_array(DES_DES).get_local_element(idx);
  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  indices_.push_back( ix );
  indices_.push_back( jx );

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
  is_local_ = spinBlock_->get_op_array(CRE).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "(C)";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_C::set_local_ops( int idx )
{
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
  factor_ = 1.0;
  if(dmrginp.doimplicitTranspose()){
    size_ = spinBlock_->get_op_array(CRE).get_size();
    is_local_ = spinBlock_->get_op_array(CRE).is_local();
    transpose_ = true;
  }
  else {
    size_ = spinBlock_->get_op_array(DES).get_size();
    is_local_ = spinBlock_->get_op_array(DES).is_local();
    transpose_= false;
  }
  build_pattern_ = "(D)";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_D::set_local_ops( int idx )
{
  indices_.clear();
  if(dmrginp.doimplicitTranspose())
    opReps_ = spinBlock_->get_op_array(CRE).get_local_element(idx);
  else
    opReps_ = spinBlock_->get_op_array(DES).get_local_element(idx);
  int ix = opReps_.at(0)->get_orbs(0);
  indices_.push_back(ix);
  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_NULL::Npdm_op_wrapper_NULL(SpinBlock * spinBlock)
{
  // Null operator
  spinBlock_ = spinBlock;
  opReps_.clear();
  indices_.clear();
  size_ = 1; // For compatibility with rest of code
  is_local_ = true;
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_NULL::set_local_ops( int idx )
{
  // Null operator
  indices_.clear();
  opReps_ = spinBlock_->get_op_array(OVERLAP).get_local_element(idx);
  return false;
}

//===========================================================================================================================================================

}
}
