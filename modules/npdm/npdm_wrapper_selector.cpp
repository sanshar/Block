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
// Initialize 4-index operators built using RI approximation (exact on dot block)
boost::shared_ptr<NpdmSpinOps> init_RI_4_index_operators( SpinBlock * spinBlock, std::vector<Npdm::CD> & cd_type ) {

  std::vector<Npdm::CD> op;

  op = { Npdm::CREATION, Npdm::CREATION, Npdm::DESTRUCTION, Npdm::DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_compound_CCDD( spinBlock ) );
    return ret;
  } 
  assert(false);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Initialize 4-index operators
boost::shared_ptr<NpdmSpinOps> init_4_index_operators( SpinBlock * spinBlock, std::vector<Npdm::CD> & cd_type ) {

  std::vector<Npdm::CD> op;

  op = { Npdm::CREATION, Npdm::CREATION, Npdm::DESTRUCTION, Npdm::DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CCDD( spinBlock ) );
    return ret;
  } 
  assert(false);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Initialize 3-index operators built using RI approximation (exact on dot block)
boost::shared_ptr<NpdmSpinOps> init_RI_3_index_operators( SpinBlock * spinBlock, std::vector<Npdm::CD> & cd_type ) {

  std::vector<Npdm::CD> op;

  op = { Npdm::CREATION, Npdm::CREATION, Npdm::DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_compound_CCD( spinBlock ) );
    return ret;
  } 
  op = { Npdm::CREATION, Npdm::DESTRUCTION, Npdm::DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_compound_CDD( spinBlock ) );
    return ret;
  } 
  op = { Npdm::CREATION, Npdm::DESTRUCTION, Npdm::CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_compound_CDC( spinBlock ) );
    return ret;
  } 
  op = { Npdm::CREATION, Npdm::CREATION, Npdm::CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_compound_CCC( spinBlock ) );
    return ret;
  } 
  op = { Npdm::DESTRUCTION, Npdm::CREATION, Npdm::DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_compound_DCD( spinBlock ) );
    return ret;
  } 
  assert(false);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Initialize 3-index operators
boost::shared_ptr<NpdmSpinOps> init_3_index_operators( SpinBlock * spinBlock, std::vector<Npdm::CD> & cd_type ) {

  std::vector<Npdm::CD> op;

  op = { Npdm::CREATION, Npdm::CREATION, Npdm::DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CCD( spinBlock ) );
    return ret;
  } 
  op = { Npdm::CREATION, Npdm::DESTRUCTION, Npdm::DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CDD( spinBlock ) );
    return ret;
  } 
  op = { Npdm::CREATION, Npdm::DESTRUCTION, Npdm::CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CDC( spinBlock ) );
    return ret;
  } 
  op = { Npdm::CREATION, Npdm::CREATION, Npdm::CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CCC( spinBlock ) );
    return ret;
  } 
  op = { Npdm::DESTRUCTION, Npdm::CREATION, Npdm::DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_DCD( spinBlock ) );
    return ret;
  } 
  op = { Npdm::DESTRUCTION, Npdm::DESTRUCTION, Npdm::CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_DDC( spinBlock ) );
    return ret;
  } 
  op = { Npdm::DESTRUCTION, Npdm::DESTRUCTION, Npdm::DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_DDD( spinBlock ) );
    return ret;
  } 
  op = { Npdm::DESTRUCTION, Npdm::CREATION, Npdm::CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_DCC( spinBlock ) );
    return ret;
  } 
  assert(false);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Initialize 2-index operators
boost::shared_ptr<NpdmSpinOps> init_2_index_operators( SpinBlock * spinBlock, std::vector<Npdm::CD> & cd_type ) {

  std::vector<Npdm::CD> op;

  op = { Npdm::CREATION, Npdm::CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CC( spinBlock ) );
    return ret;
  } 
  op = { Npdm::CREATION, Npdm::DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CD( spinBlock ) );
    return ret;
  } 
  op = { Npdm::DESTRUCTION, Npdm::CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_DC( spinBlock ) );
    return ret;
  } 
  op = { Npdm::DESTRUCTION, Npdm::DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_DD( spinBlock ) );
    return ret;
  } 
  assert(false);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Initialize 1-index operators
boost::shared_ptr<NpdmSpinOps> init_1_index_operators( SpinBlock * spinBlock, std::vector<Npdm::CD> & cd_type ) {

  std::vector<Npdm::CD> op;

  op = { Npdm::CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_C( spinBlock ) );
    return ret;
  } 
  op = { Npdm::DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_D( spinBlock ) );
    return ret;
  } 
  assert(false);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Initialize NULL operators
boost::shared_ptr<NpdmSpinOps> init_0_index_operators() {
  boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_NULL() );
  return ret;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<NpdmSpinOps> select_op_wrapper( SpinBlock * spinBlock, std::vector<Npdm::CD> & cd_type ) {

  boost::shared_ptr<NpdmSpinOps> ret;

  if ( cd_type.size() == 0 ) { ret = init_0_index_operators(); return ret; }
  if ( cd_type.size() == 1 ) { ret = init_1_index_operators( spinBlock, cd_type ); return ret; }
  if ( cd_type.size() == 2 ) { ret = init_2_index_operators( spinBlock, cd_type ); return ret; }

  if (spinBlock->size() == 1) {
    // Many-body basis is complete, so exploit RI to build many-index operators on fly (e.g. dot block)
//    if      ( cd_type.size() == 3 ) ret = init_3_index_operators( spinBlock, cd_type );
    if      ( cd_type.size() == 3 ) ret = init_RI_3_index_operators( spinBlock, cd_type );
    else if ( cd_type.size() == 4 ) ret = init_RI_4_index_operators( spinBlock, cd_type );
    else assert(false);
  }
  else {
    // Many-body basis is incomplete, so cannot exploit RI exactly
//    if      ( cd_type.size() == 3 ) ret = init_RI_3_index_operators( spinBlock, cd_type ); //FIXME only works if FCI
    if      ( cd_type.size() == 3 ) ret = init_3_index_operators( spinBlock, cd_type );
    else if ( cd_type.size() == 4 ) ret = init_4_index_operators( spinBlock, cd_type );
    else assert(false);
  }

  return ret;
}

//===========================================================================================================================================================

}

