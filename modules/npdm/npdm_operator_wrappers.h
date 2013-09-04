/*
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012
Copyright (c) 2012, Garnet K.-L. Chan

This program is integrated in Molpro with the permission of
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NPDM_OP_WRAPPERS_H
#define NPDM_OP_WRAPPERS_H

#include <boost/mpi.hpp>
#include "operatorfunctions.h"
#include "npdm_patterns.h"
//FIXME serialize bullshit
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

namespace SpinAdapted{

//===========================================================================================================================================================

class NpdmSpinOps_base {

  public:
//FIXME
    NpdmSpinOps_base() {}
// FIXME Shallow copy constructor
    NpdmSpinOps_base( const NpdmSpinOps_base & obj ) {
//      opReps_.clear();
//      for (int i=0; i < obj.opReps_.size(); ++i) {
//        boost::shared_ptr<SparseMatrix> op (new Cre);
//        opReps_.push_back(op);
//      }
//      for (int i=0; i < obj.opReps_.size(); ++i) {
//        opReps_.at(i) = obj.opReps_.at(i);
//      }
      opReps_ = obj.opReps_;
      mults_ = obj.mults_;
      build_pattern_ = obj.build_pattern_;
      transpose_ = obj.transpose_;
      factor_ = obj.factor_;
      indices_ = obj.indices_;
    }

    // Numerical representation of the operators for several total spins (e.g. 2-index op has two forms with spin-1/2 particles)
    std::vector< boost::shared_ptr<SparseMatrix> > opReps_;
    // Spin multiplicity of each operator (this info should be in each OpReps element, but we can use this for diagnostics)
    std::vector<int> mults_;
    // How the operator is built (e.g. 3-index from product of 2-index cre-cre and 1-index destruction)
    std::string build_pattern_;
    // Do we need to transpose the representation before using it?
    bool transpose_;
    // Do we need to multiply by any constant factors when using it (due to implicit use of commutation relations or such like)?
    double factor_;
    // Effective spatial orbital indices (since due to use of transposition / commutation may not match OpRep.get_orbs() etc)
    std::vector<int> indices_;

//FIXME put in implementation file
//FIXME  Do like this since serialization of whole object broken!!
    std::vector< boost::mpi::request > isend_mpi_obj( int rank, unsigned int tag_lo, unsigned int tag_hi )
    {
      boost::mpi::communicator world;
      std::vector< boost::mpi::request > reqs;
      int k = tag_lo;
      for ( int i = 0; i < opReps_.size(); ++i) {
        reqs.push_back( world.isend(rank, k++, *(opReps_.at(i))) );
      }
      reqs.push_back( world.isend(rank, k++, mults_) );
      reqs.push_back( world.isend(rank, k++, build_pattern_) );
      reqs.push_back( world.isend(rank, k++, transpose_) );
      reqs.push_back( world.isend(rank, k++, factor_) );
      reqs.push_back( world.isend(rank, k++, indices_) );
      assert( k < tag_hi );
      return reqs;
    }
      
    std::vector< boost::mpi::request > irecv_mpi_obj( int rank, unsigned int tag_lo, unsigned int tag_hi, int size )
    {
      boost::mpi::communicator world;
      std::vector< boost::mpi::request > reqs;
      assert( opReps_.size() == 0 );
      int k = tag_lo;
      for ( int i = 0; i < size; ++i) {
        boost::shared_ptr<SparseMatrix> op (new Cre);
        reqs.push_back( world.irecv(rank, k++, *op ) );
        opReps_.push_back(op);
      }
      reqs.push_back( world.irecv(rank, k++, mults_) );
      reqs.push_back( world.irecv(rank, k++, build_pattern_) );
      reqs.push_back( world.irecv(rank, k++, transpose_) );
      reqs.push_back( world.irecv(rank, k++, factor_) );
      reqs.push_back( world.irecv(rank, k++, indices_) );
      assert( k < tag_hi );
      return reqs;
    }

//  private:
//  friend class boost::serialization::access;
//  template<class Archive>
//    void serialize(Archive & ar, const unsigned int version)
//    {
////      boost::serialization::void_cast_register<Cre, SparseMatrix>();
//
//      ar & mults_ \
//         & build_pattern_ \
//         & transpose_ \
//         & factor_ \
//         & indices_;
////FIXME this is horrible!  See also spinblock.h
//      ar.register_type(static_cast<Cre *>(NULL));
//      ar.register_type(static_cast<CreDes *>(NULL));
//      ar.register_type(static_cast<CreCre *>(NULL));
//      ar.register_type(static_cast<CreDesComp *>(NULL));
//      ar.register_type(static_cast<DesDesComp *>(NULL));
//      ar.register_type(static_cast<CreCreDesComp *>(NULL));
//      ar.register_type(static_cast<Ham *>(NULL));
////MAW 3PDM
//      ar.register_type(static_cast<DesCre *>(NULL));
//      ar.register_type(static_cast<CreCreDes *>(NULL));
//      ar.register_type(static_cast<CreDesDes *>(NULL));
//      ar.register_type(static_cast<CreDesCre *>(NULL));
//      ar.register_type(static_cast<CreCreCre *>(NULL));
////MAW 4PDM
//      ar.register_type(static_cast<DesCreDes *>(NULL));
//      ar.register_type(static_cast<DesDesCre *>(NULL));
//      ar & opReps_;
//    }

};

//===========================================================================================================================================================
//FIXME constructors / destructors
class NpdmSpinOps : public NpdmSpinOps_base {

  public:
    int size() { return size_; }
    virtual bool set_local_ops( int idx ) { assert(false); }

    // Input file stream for disk-based operators used to build NPDM
    std::ifstream ifs_;
    // Number of spatial orbital combinations
    int size_;

  protected:
    boost::shared_ptr<SparseMatrix> build_compound_operator( bool is_fermion, int sign,
                                                             boost::shared_ptr<SparseMatrix> lhsOp,
                                                             boost::shared_ptr<SparseMatrix> rhsOp,
                                                             int ispin, std::vector<int> indices, bool transpose );

    SpinBlock* spinBlock_;
};

////===========================================================================================================================================================
////FIXME constructors / destructors
//class NpdmSpinOps {
//
//  private:
//  friend class boost::serialization::access;
//  template<class Archive>
//    void serialize(Archive & ar, const unsigned int version)
//    {
//cout << "serializing class NpdmSpinOps base\n";
//      // Note we don't archive ifs_ or spinBlock_
//      boost::serialization::void_cast_register<Cre, SparseMatrix>();
//      ar & factor_ \
//         & transpose_;
////      ar & mults_ \
////         & build_pattern_ \
////         & transpose_ \
////         & factor_ \
////         & build_pattern_ \
////         & size_;
////      ar.register_type(static_cast<Cre *>(NULL));
////      ar & opReps_;
//    }
////         & ifs_ \
////         & spinBlock_ \
//
//  public:
//
////FIXME filename for disk-based storage
//    int size() { return size_; }
//    virtual bool set_local_ops( int idx ) { assert(false); }
//    // Input file stream for disk-based operators used to build NPDM
//    std::ifstream ifs_;
//
//    // Numerical representation of the operators for several total spins (e.g. 2-index op has two forms with spin-1/2 particles)
//    std::vector< boost::shared_ptr<SparseMatrix> > opReps_;
//    // Spin multiplicity of each operator (this info should be in each OpReps element, but we can use this for diagnostics)
//    std::vector<int> mults_;
//    // How the operator is built (e.g. 3-index from product of 2-index cre-cre and 1-index destruction, otherwise)
//    std::string build_pattern_;
//    // Do we need to transpose the representation before using it?
//    bool transpose_;
//    // Do we need to multiply by any constant factors when using it (due to implicit use of commutation relations or such like)?
//    double factor_;
//    // Effective spatial orbital indices (since due to use of transposition / commutation may not match OpRep.get_orbs() etc)
//    std::vector<int> indices_;
//    // Number of spatial orbital combinations
//    int size_;
//
//  protected:
//    boost::shared_ptr<SparseMatrix> build_compound_operator( bool is_fermion, int sign,
//                                                             boost::shared_ptr<SparseMatrix> lhsOp,
//                                                             boost::shared_ptr<SparseMatrix> rhsOp,
//                                                             int ispin, std::vector<int> indices, bool transpose );
//
//    SpinBlock* spinBlock_;
//};
//
////===========================================================================================================================================================
//  4-INDEX compound Ops (build using RI approximation, exact on dot block)
//===========================================================================================================================================================

class Npdm_op_wrapper_compound_CCDD : public NpdmSpinOps {
  public:
//FIXME constructors / destructors??
    Npdm_op_wrapper_compound_CCDD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
};

//===========================================================================================================================================================
//  4-INDEX Ops
//===========================================================================================================================================================

class Npdm_op_wrapper_CCDD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CCDD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
};

//===========================================================================================================================================================
//  3-INDEX compound Ops (build using RI approximation, exact on dot block)
//===========================================================================================================================================================

class Npdm_op_wrapper_compound_CCD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_compound_CCD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_compound_CDD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_compound_CDD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_compound_CDC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_compound_CDC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_compound_CCC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_compound_CCC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_compound_DCD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_compound_DCD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
};

//===========================================================================================================================================================
//  3-INDEX Ops
//===========================================================================================================================================================

class Npdm_op_wrapper_CCC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CCC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CCD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CCD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CDD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CDD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CDC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CDC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_DCD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_DCD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_DDC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_DDC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_DCC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_DCC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_DDD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_DDD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
};

//===========================================================================================================================================================
//  2-INDEX Ops
//===========================================================================================================================================================

class Npdm_op_wrapper_CC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_DC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_DC( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_DD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_DD( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
};

//===========================================================================================================================================================
//  1-INDEX Ops
//===========================================================================================================================================================

class Npdm_op_wrapper_C : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_C( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_D : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_D( SpinBlock * spinBlock );
    bool set_local_ops( int idx );
};

//===========================================================================================================================================================
//  NULL case (for empty creation-destruction patterns)
//===========================================================================================================================================================

class Npdm_op_wrapper_NULL : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_NULL();
    bool set_local_ops( int idx );
};

//===========================================================================================================================================================

//BOOST_CLASS_EXPORT( NpdmSpinOps );
//BOOST_CLASS_EXPORT( Npdm_op_wrapper_compound_CCDD );
//BOOST_CLASS_EXPORT( Npdm_op_wrapper_CCDD );
//BOOST_CLASS_EXPORT( Npdm_op_wrapper_compound_CCD );
//BOOST_CLASS_EXPORT( Npdm_op_wrapper_compound_CDD );
//BOOST_CLASS_EXPORT( Npdm_op_wrapper_compound_CDC );
//BOOST_CLASS_EXPORT( Npdm_op_wrapper_compound_CCC );
//BOOST_CLASS_EXPORT( Npdm_op_wrapper_compound_DCD );
//BOOST_CLASS_EXPORT( Npdm_op_wrapper_CCC );
//BOOST_CLASS_EXPORT( Npdm_op_wrapper_CCD );
//BOOST_CLASS_EXPORT( Npdm_op_wrapper_CDD );
//BOOST_CLASS_EXPORT( Npdm_op_wrapper_CDC );
//BOOST_CLASS_EXPORT( Npdm_op_wrapper_DCD );
//BOOST_CLASS_EXPORT( Npdm_op_wrapper_DDC );
//BOOST_CLASS_EXPORT( Npdm_op_wrapper_DCC );
//BOOST_CLASS_EXPORT( Npdm_op_wrapper_DDD );
//BOOST_CLASS_EXPORT( Npdm_op_wrapper_CC );
//BOOST_CLASS_EXPORT( Npdm_op_wrapper_CD );
//BOOST_CLASS_EXPORT( Npdm_op_wrapper_DC );
//BOOST_CLASS_EXPORT( Npdm_op_wrapper_DD );
//BOOST_CLASS_EXPORT( Npdm_op_wrapper_C );
//BOOST_CLASS_EXPORT( Npdm_op_wrapper_D );
//BOOST_CLASS_EXPORT( Npdm_op_wrapper_NULL );
//

}

#endif

