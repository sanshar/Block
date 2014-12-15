#ifndef NPDM_OPERATORS_H
#define NPDM_OPERATORS_H

#include "npdm_spin_ops.h"

namespace SpinAdapted{
namespace Npdm{

//FIXME constructors / destructors??

//===========================================================================================================================================================
//  4-INDEX compound Ops (built using RI approximation, exact on dot block)
//===========================================================================================================================================================

class Npdm_op_compound_CCDD : public NpdmSpinOps {
  public:
    Npdm_op_compound_CCDD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_4INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_CCCD : public NpdmSpinOps {
  public:
    Npdm_op_compound_CCCD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_4INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_CCDC : public NpdmSpinOps {
  public:
    Npdm_op_compound_CCDC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_4INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_CDCC : public NpdmSpinOps {
  public:
    Npdm_op_compound_CDCC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_4INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_CDCD : public NpdmSpinOps {
  public:
    Npdm_op_compound_CDCD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_4INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_CDDC : public NpdmSpinOps {
  public:
    Npdm_op_compound_CDDC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_4INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_CDDD : public NpdmSpinOps {
  public:
    Npdm_op_compound_CDDD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_4INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_CCCC : public NpdmSpinOps {
  public:
    Npdm_op_compound_CCCC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_4INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_DCCD : public NpdmSpinOps {
  public:
    Npdm_op_compound_DCCD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_4INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_DCDC : public NpdmSpinOps {
  public:
    Npdm_op_compound_DCDC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_4INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_DDCC : public NpdmSpinOps {
  public:
    Npdm_op_compound_DDCC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_4INDEX).get_local_indices(); }
};

//===========================================================================================================================================================
//  4-INDEX Ops
//===========================================================================================================================================================

class Npdm_op_wrapper_CCDD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CCDD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_CRE_DES_DES).get_local_indices(); }
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_CRE_DES_DES).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CDCD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CDCD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_DES_CRE_DES).get_local_indices(); }
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_DES_CRE_DES).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CDDC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CDDC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_DES_DES_CRE).get_local_indices(); }
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_DES_DES_CRE).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CDDD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CDDD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_DES_DES_DES).get_local_indices(); }
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_DES_DES_DES).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CCCD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CCCD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_CRE_CRE_DES).get_local_indices(); }
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_CRE_CRE_DES).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CCDC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CCDC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_CRE_DES_CRE).get_local_indices(); }
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_CRE_DES_CRE).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CDCC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CDCC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_DES_CRE_CRE).get_local_indices(); }
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_DES_CRE_CRE).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CCCC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CCCC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_CRE_CRE_CRE).get_local_indices(); }
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_CRE_CRE_CRE).get_array(); }
};

//===========================================================================================================================================================
//  3-INDEX compound Ops (built using RI approximation, exact on dot block)
//===========================================================================================================================================================

class Npdm_op_compound_CCD : public NpdmSpinOps {
  public:
    Npdm_op_compound_CCD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_3INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_CDD : public NpdmSpinOps {
  public:
    Npdm_op_compound_CDD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_3INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_CDC : public NpdmSpinOps {
  public:
    Npdm_op_compound_CDC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_3INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_CCC : public NpdmSpinOps {
  public:
    Npdm_op_compound_CCC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_3INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_DCD : public NpdmSpinOps {
  public:
    Npdm_op_compound_DCD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_3INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_DDC : public NpdmSpinOps {
  public:
    Npdm_op_compound_DDC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_3INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_DCC : public NpdmSpinOps {
  public:
    Npdm_op_compound_DCC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_3INDEX).get_local_indices(); }
};

//===========================================================================================================================================================
//  3-INDEX Ops
//===========================================================================================================================================================

class Npdm_op_wrapper_CCC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CCC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_CRE_CRE).get_local_indices(); }
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_CRE_CRE).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CCD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CCD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_CRE_DES).get_local_indices(); }
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_CRE_DES).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CDD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CDD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_DES_DES).get_local_indices(); }
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_DES_DES).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CDC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CDC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_DES_CRE).get_local_indices(); }
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_DES_CRE).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_DCD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_DCD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(DES_CRE_DES).get_local_indices(); }
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(DES_CRE_DES).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_DDC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_DDC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(DES_DES_CRE).get_local_indices(); }
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(DES_DES_CRE).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_DCC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_DCC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(DES_DES_CRE).get_local_indices(); }
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(DES_DES_CRE).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_DDD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_DDD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_CRE_CRE).get_local_indices(); }
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_CRE_CRE).get_array(); }
};

//===========================================================================================================================================================
//  2-INDEX Ops
//===========================================================================================================================================================

class Npdm_op_wrapper_CC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_CRE).get_local_indices(); }
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_CRE).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_DES).get_local_indices(); }
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_DES).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_DC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_DC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(DES_CRE).get_local_indices(); }
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(DES_CRE).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_DD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_DD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(DES_DES).get_local_indices(); }
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(DES_DES).get_array(); }
};

//===========================================================================================================================================================
//  1-INDEX Ops
//===========================================================================================================================================================

class Npdm_op_wrapper_C : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_C( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE).get_local_indices(); }
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_D : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_D( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE).get_local_indices(); }
//    FIXME

      std::vector< std::vector<int> > get_indices() { 
        if(dmrginp.doimplicitTranspose())
          return spinBlock_->get_op_array(CRE).get_array(); 
        else
          return spinBlock_->get_op_array(DES).get_array(); 
      }
};

//===========================================================================================================================================================
//  NULL case (for empty creation-destruction patterns)
//===========================================================================================================================================================

class Npdm_op_wrapper_NULL : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_NULL(SpinBlock * spinBlock);
    bool set_local_ops( int idx );
};

//===========================================================================================================================================================

}
}

#endif

