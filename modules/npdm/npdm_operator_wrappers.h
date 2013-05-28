/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NPDM_OP_WRAPPERS_H
#define NPDM_OP_WRAPPERS_H

#include "operatorfunctions.h"
#include "npdm_patterns.h"

namespace SpinAdapted{

//===========================================================================================================================================================
//  BASE CLASS
//===========================================================================================================================================================
//FIXME constructors / destructors
class NpdmSpinOps {

  public:
    int size() { return size_; };
    virtual void set_local_ops( int idx ) { assert(false); };

    // Numerical representation of the operators for several total spins (e.g. 2-index op has two forms with spin-1/2 particles)
    std::vector< boost::shared_ptr<SparseMatrix> > opReps_;
    // Spin multiplicity of each operator (this info should be in each OpReps element, but we can use this for diagnostics)
    std::vector<int> mults_;
    // How the operator is built (e.g. 3-index from product of 2-index cre-cre and 1-index destruction, otherwise)
    std::string build_pattern_;
    // Do we need to transpose the representation before using it?
    bool transpose_;
    // Do we need to multiply by any constant factors when using it (due to implicit use of commutation relations)?
    double factor_;
    // Effective spatial orbital indices (since due to use of transposition / commutation may not match OpRep.get_orbs() etc)
    std::vector<int> indices_;

  protected:
    boost::shared_ptr<SparseMatrix> build_compound_operator( bool is_fermion,
                                                             boost::shared_ptr<SparseMatrix> lhsOp,
                                                             boost::shared_ptr<SparseMatrix> rhsOp,
                                                             int ispin, std::vector<int> indices, bool transpose );

    SpinBlock * spinBlock_;
    // Number of spatial orbital combinations
    int size_;

};

//===========================================================================================================================================================
//  4-INDEX compound Ops
//===========================================================================================================================================================

class Npdm_op_wrapper_compound_CCDD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_compound_CCDD( SpinBlock * spinBlock );
    void set_local_ops( int idx );
};

//===========================================================================================================================================================
//  3-INDEX compound Ops
//===========================================================================================================================================================

class Npdm_op_wrapper_compound_CCD : public NpdmSpinOps {
  public: 
    Npdm_op_wrapper_compound_CCD( SpinBlock * spinBlock );
    void set_local_ops( int idx );
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_compound_CDD : public NpdmSpinOps {
  public: 
    Npdm_op_wrapper_compound_CDD( SpinBlock * spinBlock );
    void set_local_ops( int idx );
};

//===========================================================================================================================================================
//  3-INDEX Ops
//===========================================================================================================================================================

class Npdm_op_wrapper_CCD : public NpdmSpinOps {
  public: 
    Npdm_op_wrapper_CCD( SpinBlock * spinBlock );
    void set_local_ops( int idx );
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CDD : public NpdmSpinOps {
  public: 
    Npdm_op_wrapper_CDD( SpinBlock * spinBlock );
    void set_local_ops( int idx );
};

//===========================================================================================================================================================
//  2-INDEX Ops
//===========================================================================================================================================================

class Npdm_op_wrapper_CC : public NpdmSpinOps {
  public: 
    Npdm_op_wrapper_CC( SpinBlock * spinBlock );
    void set_local_ops( int idx );
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CD : public NpdmSpinOps {
  public: 
    Npdm_op_wrapper_CD( SpinBlock * spinBlock );
    void set_local_ops( int idx );
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_DD : public NpdmSpinOps {
  public: 
    Npdm_op_wrapper_DD( SpinBlock * spinBlock );
    void set_local_ops( int idx );
};

//===========================================================================================================================================================
//  1-INDEX Ops
//===========================================================================================================================================================

class Npdm_op_wrapper_C : public NpdmSpinOps {
  public: 
    Npdm_op_wrapper_C( SpinBlock * spinBlock );
    void set_local_ops( int idx );
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_D : public NpdmSpinOps {
  public: 
    Npdm_op_wrapper_D( SpinBlock * spinBlock );
    void set_local_ops( int idx );
};

//===========================================================================================================================================================
//  NULL case (for empty creation-destruction patterns)
//===========================================================================================================================================================

class Npdm_op_wrapper_NULL : public NpdmSpinOps {
  public: 
    Npdm_op_wrapper_NULL();
    void set_local_ops( int idx );
};

//===========================================================================================================================================================

}

#endif
 
