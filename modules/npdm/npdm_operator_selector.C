/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "npdm_patterns.h"
#include "npdm_spin_ops.h"
#include "npdm_operators.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================
// Initialize 4-index operators built using RI approximation (exact on 1-site block)
boost::shared_ptr<NpdmSpinOps> init_RI_4_index_operators( SpinBlock * spinBlock,const std::vector<CD> & cd_type ) {

  std::vector<CD> op;

  op = { CREATION, CREATION, DESTRUCTION, DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_compound_CCDD( spinBlock ) );
    return ret;
  } 
  op = { CREATION, CREATION, CREATION, DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_compound_CCCD( spinBlock ) );
    return ret;
  } 
  op = { CREATION, CREATION, DESTRUCTION, CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_compound_CCDC( spinBlock ) );
    return ret;
  } 
  op = { CREATION, DESTRUCTION, CREATION, CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_compound_CDCC( spinBlock ) );
    return ret;
  } 
  op = { CREATION, DESTRUCTION, CREATION, DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_compound_CDCD( spinBlock ) );
    return ret;
  } 
  op = { CREATION, DESTRUCTION, DESTRUCTION, CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_compound_CDDC( spinBlock ) );
    return ret;
  } 
  op = { CREATION, DESTRUCTION, DESTRUCTION, DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_compound_CDDD( spinBlock ) );
    return ret;
  } 
  op = { CREATION, CREATION, CREATION, CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_compound_CCCC( spinBlock ) );
    return ret;
  } 
  abort();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Initialize 4-index operators
boost::shared_ptr<NpdmSpinOps> init_4_index_operators( SpinBlock * spinBlock,const std::vector<CD> & cd_type ) {

  std::vector<CD> op;

  op = { CREATION, CREATION, DESTRUCTION, DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CCDD( spinBlock ) );
    return ret;
  } 
  op = { CREATION, DESTRUCTION, CREATION, DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CDCD( spinBlock ) );
    return ret;
  } 
  op = { CREATION, DESTRUCTION, DESTRUCTION, CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CDDC( spinBlock ) );
    return ret;
  } 
  op = { CREATION, DESTRUCTION, DESTRUCTION, DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CDDD( spinBlock ) );
    return ret;
  } 
  op = { CREATION, CREATION, CREATION, DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CCCD( spinBlock ) );
    return ret;
  } 
  op = { CREATION, CREATION, DESTRUCTION, CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CCDC( spinBlock ) );
    return ret;
  } 
  op = { CREATION, DESTRUCTION, CREATION, CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CDCC( spinBlock ) );
    return ret;
  } 
  op = { CREATION, CREATION, CREATION, CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CCCC( spinBlock ) );
    return ret;
  } 
  abort();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Initialize 3-index operators built using RI approximation (exact on 1-site block)
boost::shared_ptr<NpdmSpinOps> init_RI_3_index_operators( SpinBlock * spinBlock,const std::vector<CD> & cd_type ) {

  std::vector<CD> op;

  op = { CREATION, CREATION, DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_compound_CCD( spinBlock ) );
    return ret;
  } 
  op = { CREATION, DESTRUCTION, DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_compound_CDD( spinBlock ) );
    return ret;
  } 
  op = { CREATION, DESTRUCTION, CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_compound_CDC( spinBlock ) );
    return ret;
  } 
  op = { CREATION, CREATION, CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_compound_CCC( spinBlock ) );
    return ret;
  } 
  op = { DESTRUCTION, CREATION, DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_compound_DCD( spinBlock ) );
    return ret;
  } 
  op = { DESTRUCTION, DESTRUCTION, CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_compound_DDC( spinBlock ) );
    return ret;
  } 
  op = { DESTRUCTION, CREATION, CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_compound_DCC( spinBlock ) );
    return ret;
  } 
  abort();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Initialize 3-index operators
boost::shared_ptr<NpdmSpinOps> init_3_index_operators( SpinBlock * spinBlock,const std::vector<CD> & cd_type ) {

  std::vector<CD> op;

  op = { CREATION, CREATION, DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CCD( spinBlock ) );
    return ret;
  } 
  op = { CREATION, DESTRUCTION, DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CDD( spinBlock ) );
    return ret;
  } 
  op = { CREATION, DESTRUCTION, CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CDC( spinBlock ) );
    return ret;
  } 
  op = { CREATION, CREATION, CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CCC( spinBlock ) );
    return ret;
  } 
  op = { DESTRUCTION, CREATION, DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_DCD( spinBlock ) );
    return ret;
  } 
  op = { DESTRUCTION, DESTRUCTION, CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_DDC( spinBlock ) );
    return ret;
  } 
  op = { DESTRUCTION, DESTRUCTION, DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_DDD( spinBlock ) );
    return ret;
  } 
  op = { DESTRUCTION, CREATION, CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_DCC( spinBlock ) );
//FIXME
//    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_compound_DCC( spinBlock ) );
    return ret;
  } 
  abort();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Initialize 2-index operators
boost::shared_ptr<NpdmSpinOps> init_2_index_operators( SpinBlock * spinBlock, const std::vector<CD> & cd_type ) {

  std::vector<CD> op;

  op = { CREATION, CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CC( spinBlock ) );
    return ret;
  } 
  op = { CREATION, DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_CD( spinBlock ) );
    return ret;
  } 
  op = { DESTRUCTION, CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_DC( spinBlock ) );
    return ret;
  } 
  op = { DESTRUCTION, DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_DD( spinBlock ) );
    return ret;
  } 
  abort();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Initialize 1-index operators
boost::shared_ptr<NpdmSpinOps> init_1_index_operators( SpinBlock * spinBlock,const std::vector<CD> & cd_type ) {

  std::vector<CD> op;

  op = { CREATION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_C( spinBlock ) );
    return ret;
  } 
  op = { DESTRUCTION };
  if ( cd_type == op ) {
    boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_D( spinBlock ) );
    return ret;
  } 
  abort();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Initialize NULL operators
boost::shared_ptr<NpdmSpinOps> init_0_index_operators(SpinBlock * spinBlock) {
  boost::shared_ptr<NpdmSpinOps> ret( new Npdm_op_wrapper_NULL(spinBlock) );
  return ret;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<NpdmSpinOps> select_op_wrapper( SpinBlock * spinBlock,const std::vector<CD> & cd_type ) {

  boost::shared_ptr<NpdmSpinOps> ret;

  if ( cd_type.size() == 0 ) { ret = init_0_index_operators( spinBlock); return ret; }
  if ( cd_type.size() == 1 ) { ret = init_1_index_operators( spinBlock, cd_type ); return ret; }
  if ( cd_type.size() == 2 ) { ret = init_2_index_operators( spinBlock, cd_type ); return ret; }

  if (spinBlock->size() == 1) {
    // Many-body basis is complete, so exploit RI to build many-index operators on fly (e.g. dot block)
    if      ( cd_type.size() == 3 ) ret = init_RI_3_index_operators( spinBlock, cd_type );
    else if ( cd_type.size() == 4 ) ret = init_RI_4_index_operators( spinBlock, cd_type );
    else abort();
  }
  else {
    //FIXME
    if(!dmrginp.spinAdapted() && spinBlock->size()==2 ){
      // Many-body basis is complete, so exploit RI to build many-index operators on fly (e.g. dot block)
      if      ( cd_type.size() == 3 ) ret = init_RI_3_index_operators( spinBlock, cd_type );
      else if ( cd_type.size() == 4 ) ret = init_RI_4_index_operators( spinBlock, cd_type );
      else abort();

    }
    else{
      // Many-body basis is incomplete, so cannot exploit RI exactly
      if      ( cd_type.size() == 3 ) ret = init_3_index_operators( spinBlock, cd_type );
      else if ( cd_type.size() == 4 ) ret = init_4_index_operators( spinBlock, cd_type );
      else abort();
    }
  }

  return ret;
}

//===========================================================================================================================================================

}
}

