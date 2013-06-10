/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "npdm_operator_wrappers.h"

namespace SpinAdapted{

//===========================================================================================================================================================
// Init 4-index (dot) operators
boost::shared_ptr<NpdmSpinOps> init_4_index_operators( SpinBlock * spinBlock, std::vector<Npdm::CD> & cd_type ) {

  std::vector<Npdm::CD> op;

  op = { Npdm::CREATION, Npdm::CREATION, Npdm::DESTRUCTION, Npdm::DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_compound_CCDD( spinBlock ) );
    return ret;
  } 
  assert(false);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Init 3-index operators
boost::shared_ptr<NpdmSpinOps> init_3_index_operators( SpinBlock * spinBlock, std::vector<Npdm::CD> & cd_type ) {

  std::vector<Npdm::CD> op;

  op = { Npdm::CREATION, Npdm::CREATION, Npdm::DESTRUCTION };
  if ( cd_type == op ) {
//    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_compound_CCD( spinBlock ) );
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CCD( spinBlock ) );
    return ret;
  } 
  op = { Npdm::CREATION, Npdm::DESTRUCTION, Npdm::DESTRUCTION };
  if ( cd_type == op ) {
//    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_compound_CDD( spinBlock ) );
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CDD( spinBlock ) );
    return ret;
  } 
  op = { Npdm::CREATION, Npdm::DESTRUCTION, Npdm::CREATION };
  if ( cd_type == op ) {
//    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_compound_CDC( spinBlock ) );
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CDC( spinBlock ) );
    return ret;
  } 
  op = { Npdm::CREATION, Npdm::CREATION, Npdm::CREATION };
  if ( cd_type == op ) {
//    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_compound_CCC( spinBlock ) );
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CCC( spinBlock ) );
    return ret;
  } 
  assert(false);
  assert(false);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Init 2-index operators
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
// Init 1-index operators
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
// Init NULL operators
boost::shared_ptr<NpdmSpinOps> init_0_index_operators() {
  boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_NULL() );
  return ret;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<NpdmSpinOps> select_op_wrapper( SpinBlock * spinBlock, std::vector<Npdm::CD> & cd_type ) {

  boost::shared_ptr<NpdmSpinOps> ret;

  if ( cd_type.size() == 0 ) ret = init_0_index_operators();
  else if ( cd_type.size() == 1 ) ret = init_1_index_operators( spinBlock, cd_type );
  else if ( cd_type.size() == 2 ) ret = init_2_index_operators( spinBlock, cd_type );
  else if ( cd_type.size() == 3 ) ret = init_3_index_operators( spinBlock, cd_type );
  else if ( cd_type.size() == 4 ) ret = init_4_index_operators( spinBlock, cd_type );
  else assert(false);

  return ret;
}

//===========================================================================================================================================================

}

