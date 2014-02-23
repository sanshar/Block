/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "global.h"
#include "npdm_patterns.h"
#include "npdm_operator_wrappers.h"

//FIXME streamline code and reuse??  functions are almost identical here
namespace SpinAdapted{

//===========================================================================================================================================================

Npdm_op_wrapper_CCDD::Npdm_op_wrapper_CCDD( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_CRE_DES_DES).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_CRE_DES_DES).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
  // For disk-based storage
  ifile_ = spinBlock_->get_op_array(CRE_CRE_DES_DES).get_filename();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CCDD::set_local_ops( int idx )
{
cout << "getting CCDD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx, lx;

  // Read in operator representations from disk or memory
//  if ( dmrginp.do_npdm_in_core() ) {
  if ( false ) {
    opReps_ = spinBlock_->get_op_array(CRE_CRE_DES_DES).get_local_element(idx);
    build_pattern_ = opReps_.at(0)->get_build_pattern();
  }
  else {
cout << "...from disk...\n";
    // Get first operator of spin set
    assert( check_file_open( idx ) );
    boost::shared_ptr<SparseMatrix> op (new CreCreDesDes);
    boost::archive::binary_iarchive load_op(ifs_);
    load_op >> *op;
    build_pattern_ = op->get_build_pattern();
    // Get order of spin components
    ix = op->get_orbs(0); jx = op->get_orbs(1); kx = op->get_orbs(2); lx = op->get_orbs(3);
    opReps_ = spinBlock_->get_op_array(CRE_CRE_DES_DES).get_element(ix,jx,kx,lx);
    // Store spin components in correct order
    for (int j = 0; j < opReps_.size(); j++) {
      if ( opReps_[j]->get_quantum_ladder().at(build_pattern_) == op->get_quantum_ladder().at(build_pattern_) ) {
//FIXME break
        opReps_[j] = op; break;
      }
    }
    // Read in remaining spin-set from disk 
    for (int i = 1; i < opReps_.size(); i++) {
      boost::shared_ptr<SparseMatrix> op (new CreCreDesDes);
      boost::archive::binary_iarchive load_op(ifs_);
      load_op >> *op;
      for (int j = 0; j < opReps_.size(); j++) {
        // Store spin components in correct order
        if ( opReps_[j]->get_quantum_ladder().at(build_pattern_) == op->get_quantum_ladder().at(build_pattern_) ) {
          opReps_[j] = op; break;
        }
      }
    }
    assert( check_file_close( idx ) );
  }

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

Npdm_op_wrapper_CDCD::Npdm_op_wrapper_CDCD( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_DES_CRE_DES).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_DES_CRE_DES).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CDCD::set_local_ops( int idx )
{
cout << "getting CDCD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx, lx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core() ) {
    opReps_ = spinBlock_->get_op_array(CRE_DES_CRE_DES).get_local_element(idx);
    build_pattern_ = opReps_.at(0)->get_build_pattern();
  }
  else {
cout << "...from disk...\n";
    // Get first operator of spin set
    assert( check_file_open( idx ) );
    boost::shared_ptr<SparseMatrix> op (new CreDesCreDes);
    boost::archive::binary_iarchive load_op(ifs_);
    load_op >> *op;
    build_pattern_ = op->get_build_pattern();
    // Get order of spin components
    ix = op->get_orbs(0); jx = op->get_orbs(1); kx = op->get_orbs(2); lx = op->get_orbs(3);
    opReps_ = spinBlock_->get_op_array(CRE_DES_CRE_DES).get_element(ix,jx,kx,lx);
    // Store spin components in correct order
    for (int j = 0; j < opReps_.size(); j++) {
      if ( opReps_[j]->get_quantum_ladder().at(build_pattern_) == op->get_quantum_ladder().at(build_pattern_) ) {
        opReps_[j] = op;
      }
    }
    // Read in remaining spin-set from disk 
    for (int i = 1; i < opReps_.size(); i++) {
      boost::shared_ptr<SparseMatrix> op (new CreDesCreDes);
      boost::archive::binary_iarchive load_op(ifs_);
      load_op >> *op;
      for (int j = 0; j < opReps_.size(); j++) {
        // Store spin components in correct order
        if ( opReps_[j]->get_quantum_ladder().at(build_pattern_) == op->get_quantum_ladder().at(build_pattern_) ) {
          opReps_[j] = op;
        }
      }
    }
    assert( check_file_close( idx ) );
  }

  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
  lx = opReps_.at(0)->get_orbs(3);
  // Skip if commutation problematic
  if ( jx == kx ) return true;

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  indices_.push_back( lx );

  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_CDDC::Npdm_op_wrapper_CDDC( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_DES_DES_CRE).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_DES_DES_CRE).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CDDC::set_local_ops( int idx )
{
cout << "getting CDDC operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx, lx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core() ) {
    opReps_ = spinBlock_->get_op_array(CRE_DES_DES_CRE).get_local_element(idx);
    build_pattern_ = opReps_.at(0)->get_build_pattern();
  }
  else {
cout << "...from disk...\n";
    // Get first operator of spin set
    assert( check_file_open( idx ) );
    boost::shared_ptr<SparseMatrix> op (new CreDesDesCre);
    boost::archive::binary_iarchive load_op(ifs_);
    load_op >> *op;
    build_pattern_ = op->get_build_pattern();
    // Get order of spin components
    ix = op->get_orbs(0); jx = op->get_orbs(1); kx = op->get_orbs(2); lx = op->get_orbs(3);
    opReps_ = spinBlock_->get_op_array(CRE_DES_DES_CRE).get_element(ix,jx,kx,lx);
    // Store spin components in correct order
    for (int j = 0; j < opReps_.size(); j++) {
      if ( opReps_[j]->get_quantum_ladder().at(build_pattern_) == op->get_quantum_ladder().at(build_pattern_) ) {
        opReps_[j] = op;
      }
    }
    // Read in remaining spin-set from disk 
    for (int i = 1; i < opReps_.size(); i++) {
      boost::shared_ptr<SparseMatrix> op (new CreDesDesCre);
      boost::archive::binary_iarchive load_op(ifs_);
      load_op >> *op;
      for (int j = 0; j < opReps_.size(); j++) {
        // Store spin components in correct order
        if ( opReps_[j]->get_quantum_ladder().at(build_pattern_) == op->get_quantum_ladder().at(build_pattern_) ) {
          opReps_[j] = op;
        }
      }
    }
    assert( check_file_close( idx ) );
  }

  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
  lx = opReps_.at(0)->get_orbs(3);
  // Skip if commutation problematic
  if ( (kx == lx) || (jx == lx) ) return true;

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  indices_.push_back( lx );

  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_CDDD::Npdm_op_wrapper_CDDD( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_DES_DES_DES).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_DES_DES_DES).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CDDD::set_local_ops( int idx )
{
cout << "getting CDDD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx, lx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core() ) {
    opReps_ = spinBlock_->get_op_array(CRE_DES_DES_DES).get_local_element(idx);
    build_pattern_ = opReps_.at(0)->get_build_pattern();
  }
  else {
cout << "...from disk...\n";
    // Get first operator of spin set
    assert( check_file_open( idx ) );
    boost::shared_ptr<SparseMatrix> op (new CreDesDesDes);
    boost::archive::binary_iarchive load_op(ifs_);
    load_op >> *op;
    build_pattern_ = op->get_build_pattern();
    // Get order of spin components
    ix = op->get_orbs(0); jx = op->get_orbs(1); kx = op->get_orbs(2); lx = op->get_orbs(3);
    opReps_ = spinBlock_->get_op_array(CRE_DES_DES_DES).get_element(ix,jx,kx,lx);
    // Store spin components in correct order
    for (int j = 0; j < opReps_.size(); j++) {
      if ( opReps_[j]->get_quantum_ladder().at(build_pattern_) == op->get_quantum_ladder().at(build_pattern_) ) {
        opReps_[j] = op;
      }
    }
    // Read in remaining spin-set from disk 
    for (int i = 1; i < opReps_.size(); i++) {
      boost::shared_ptr<SparseMatrix> op (new CreDesDesDes);
      boost::archive::binary_iarchive load_op(ifs_);
      load_op >> *op;
      for (int j = 0; j < opReps_.size(); j++) {
        // Store spin components in correct order
        if ( opReps_[j]->get_quantum_ladder().at(build_pattern_) == op->get_quantum_ladder().at(build_pattern_) ) {
          opReps_[j] = op;
        }
      }
    }
    assert( check_file_close( idx ) );
  }

  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
  lx = opReps_.at(0)->get_orbs(3);
  // Skip if commutation problematic
  //if ( kx == lx ) return true;

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  indices_.push_back( lx );

  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_CCCD::Npdm_op_wrapper_CCCD( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_CRE_CRE_DES).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_CRE_CRE_DES).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CCCD::set_local_ops( int idx )
{
cout << "getting CCCD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx, lx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core() ) {
    opReps_ = spinBlock_->get_op_array(CRE_CRE_CRE_DES).get_local_element(idx);
    build_pattern_ = opReps_.at(0)->get_build_pattern();
  }
  else {
cout << "...from disk...\n";
    // Get first operator of spin set
    assert( check_file_open( idx ) );
    boost::shared_ptr<SparseMatrix> op (new CreCreCreDes);
    boost::archive::binary_iarchive load_op(ifs_);
    load_op >> *op;
    build_pattern_ = op->get_build_pattern();
    // Get order of spin components
    ix = op->get_orbs(0); jx = op->get_orbs(1); kx = op->get_orbs(2); lx = op->get_orbs(3);
    opReps_ = spinBlock_->get_op_array(CRE_CRE_CRE_DES).get_element(ix,jx,kx,lx);
    // Store spin components in correct order
    for (int j = 0; j < opReps_.size(); j++) {
      if ( opReps_[j]->get_quantum_ladder().at(build_pattern_) == op->get_quantum_ladder().at(build_pattern_) ) {
        opReps_[j] = op;
      }
    }
    // Read in remaining spin-set from disk 
    for (int i = 1; i < opReps_.size(); i++) {
      boost::shared_ptr<SparseMatrix> op (new CreCreCreDes);
      boost::archive::binary_iarchive load_op(ifs_);
      load_op >> *op;
      for (int j = 0; j < opReps_.size(); j++) {
        // Store spin components in correct order
        if ( opReps_[j]->get_quantum_ladder().at(build_pattern_) == op->get_quantum_ladder().at(build_pattern_) ) {
          opReps_[j] = op;
        }
      }
    }
    assert( check_file_close( idx ) );
  }

  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
  lx = opReps_.at(0)->get_orbs(3);
  // Skip if commutation problematic
  //if ( kx == lx ) return true;

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  indices_.push_back( lx );

  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_CCDC::Npdm_op_wrapper_CCDC( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_CRE_DES_CRE).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_CRE_DES_CRE).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CCDC::set_local_ops( int idx )
{
cout << "getting CCDC operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx, lx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core() ) {
    opReps_ = spinBlock_->get_op_array(CRE_CRE_DES_CRE).get_local_element(idx);
    build_pattern_ = opReps_.at(0)->get_build_pattern();
  }
  else {
cout << "...from disk...\n";
    // Get first operator of spin set
    assert( check_file_open( idx ) );
    boost::shared_ptr<SparseMatrix> op (new CreCreDesCre);
    boost::archive::binary_iarchive load_op(ifs_);
    load_op >> *op;
    build_pattern_ = op->get_build_pattern();
    // Get order of spin components
    ix = op->get_orbs(0); jx = op->get_orbs(1); kx = op->get_orbs(2); lx = op->get_orbs(3);
    opReps_ = spinBlock_->get_op_array(CRE_CRE_DES_CRE).get_element(ix,jx,kx,lx);
    // Store spin components in correct order
    for (int j = 0; j < opReps_.size(); j++) {
      if ( opReps_[j]->get_quantum_ladder().at(build_pattern_) == op->get_quantum_ladder().at(build_pattern_) ) {
        opReps_[j] = op;
      }
    }
    // Read in remaining spin-set from disk 
    for (int i = 1; i < opReps_.size(); i++) {
      boost::shared_ptr<SparseMatrix> op (new CreCreDesCre);
      boost::archive::binary_iarchive load_op(ifs_);
      load_op >> *op;
      for (int j = 0; j < opReps_.size(); j++) {
        // Store spin components in correct order
        if ( opReps_[j]->get_quantum_ladder().at(build_pattern_) == op->get_quantum_ladder().at(build_pattern_) ) {
          opReps_[j] = op;
        }
      }
    }
    assert( check_file_close( idx ) );
  }

  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
  lx = opReps_.at(0)->get_orbs(3);
  // Skip if commutation problematic
  if ( kx == lx ) return true;

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  indices_.push_back( lx );

  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_CDCC::Npdm_op_wrapper_CDCC( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_DES_CRE_CRE).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_DES_CRE_CRE).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CDCC::set_local_ops( int idx )
{
cout << "getting CDCC operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx, lx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core() ) {
    opReps_ = spinBlock_->get_op_array(CRE_DES_CRE_CRE).get_local_element(idx);
    build_pattern_ = opReps_.at(0)->get_build_pattern();
  }
  else {
cout << "...from disk...\n";
    // Get first operator of spin set
    assert( check_file_open( idx ) );
    boost::shared_ptr<SparseMatrix> op (new CreDesCreCre);
    boost::archive::binary_iarchive load_op(ifs_);
    load_op >> *op;
    build_pattern_ = op->get_build_pattern();
    // Get order of spin components
    ix = op->get_orbs(0); jx = op->get_orbs(1); kx = op->get_orbs(2); lx = op->get_orbs(3);
    opReps_ = spinBlock_->get_op_array(CRE_DES_CRE_CRE).get_element(ix,jx,kx,lx);
    // Store spin components in correct order
    for (int j = 0; j < opReps_.size(); j++) {
      if ( opReps_[j]->get_quantum_ladder().at(build_pattern_) == op->get_quantum_ladder().at(build_pattern_) ) {
        opReps_[j] = op;
      }
    }
    // Read in remaining spin-set from disk 
    for (int i = 1; i < opReps_.size(); i++) {
      boost::shared_ptr<SparseMatrix> op (new CreDesCreCre);
      boost::archive::binary_iarchive load_op(ifs_);
      load_op >> *op;
      for (int j = 0; j < opReps_.size(); j++) {
        // Store spin components in correct order
        if ( opReps_[j]->get_quantum_ladder().at(build_pattern_) == op->get_quantum_ladder().at(build_pattern_) ) {
          opReps_[j] = op;
        }
      }
    }
    assert( check_file_close( idx ) );
  }

  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
  lx = opReps_.at(0)->get_orbs(3);
  // Skip if commutation problematic
  if ( (jx==kx) || (jx==lx) ) return true;

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  indices_.push_back( lx );

  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_CCCC::Npdm_op_wrapper_CCCC( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_CRE_CRE_CRE).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_CRE_CRE_CRE).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CCCC::set_local_ops( int idx )
{
cout << "getting CCCC operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx, lx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core() ) {
    opReps_ = spinBlock_->get_op_array(CRE_CRE_CRE_CRE).get_local_element(idx);
    build_pattern_ = opReps_.at(0)->get_build_pattern();
  }
  else {
cout << "...from disk...\n";
    // Get first operator of spin set
    assert( check_file_open( idx ) );
    boost::shared_ptr<SparseMatrix> op (new CreCreCreCre);
    boost::archive::binary_iarchive load_op(ifs_);
    load_op >> *op;
    build_pattern_ = op->get_build_pattern();
    // Get order of spin components
    ix = op->get_orbs(0); jx = op->get_orbs(1); kx = op->get_orbs(2); lx = op->get_orbs(3);
    opReps_ = spinBlock_->get_op_array(CRE_CRE_CRE_CRE).get_element(ix,jx,kx,lx);
    // Store spin components in correct order
    for (int j = 0; j < opReps_.size(); j++) {
      if ( opReps_[j]->get_quantum_ladder().at(build_pattern_) == op->get_quantum_ladder().at(build_pattern_) ) {
        opReps_[j] = op;
      }
    }
    // Read in remaining spin-set from disk 
    for (int i = 1; i < opReps_.size(); i++) {
      boost::shared_ptr<SparseMatrix> op (new CreCreCreCre);
      boost::archive::binary_iarchive load_op(ifs_);
      load_op >> *op;
      for (int j = 0; j < opReps_.size(); j++) {
        // Store spin components in correct order
        if ( opReps_[j]->get_quantum_ladder().at(build_pattern_) == op->get_quantum_ladder().at(build_pattern_) ) {
          opReps_[j] = op;
        }
      }
    }
    assert( check_file_close( idx ) );
  }

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

