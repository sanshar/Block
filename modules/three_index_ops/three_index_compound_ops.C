/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "npdm_operators.h"
#include "npdm_spin_ops.h"
#include "pario.h"

//FIXME update all these to allow RI with any orbitals, not just 1-site

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

Npdm_op_compound_CCD::Npdm_op_compound_CCD( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
//FIXMEMAW
//  size_ = 1;
  size_ = spinBlock_->get_op_array(RI_3INDEX).get_size();
//  if(!dmrginp.spinAdapted()) size_=1;
  is_local_ = true;
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((CC)D)";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_compound_CCD::set_local_ops( int idx )
{
  if(dmrginp.doimplicitTranspose()){
//pout << "getting compound CCD operator...\n";
//pout << "size_ = " << size_ << endl;
//assert( idx == 0 );
  // Spatial orbital indices
//FIXME don't need to keep reconstructing whole array!
  indices_ = spinBlock_->get_op_array(RI_3INDEX).get_array().at(idx);
  int ix = indices_[0];
  int jx = indices_[1];
  int kx = indices_[2];
  //pout <<ix<<','<<jx<<','<<kx<<endl;
  if(!dmrginp.spinAdapted()){
    // only suitble for single site block. ( there are two sites in single site block in non-spinadapted dmrg.)
    if(ix==jx) return true;
   // if(ix==kx) return true;
  }

  // Get 2-index and 1-index ops as RI building blocks
  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_CRE).get_element(ix,jx);
  //std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(DES).get_element(kx);
  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(CRE).get_element(kx);

////  // Spatial orbital indices
////  indices_.clear();
////  int ix, jx, kx;
////
////  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_CRE).get_local_element(idx);
////  ix = twoOps.at(0)->get_orbs(0);
////  jx = twoOps.at(0)->get_orbs(1);
////  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(CRE).get_local_element(idx);
////  kx = oneOp.at(0)->get_orbs(0);
////  indices_.push_back( ix );
////  indices_.push_back( jx );
////  indices_.push_back( kx );
////pout << "CCD indices = " << ix << ", " << jx << ", " << kx << endl;

  // Allocate and build operator representation on the fly as RI tensor product for each spin component
  opReps_.clear();
  if(dmrginp.spinAdapted()){
    // S=0 (+) S=1/2  =>  S=1/2
    opReps_.push_back( build_compound_operator( true, -1, twoOps.at(0), oneOp.at(0), 0, indices_, true ) );
    // S=1 (+) S=1/2  =>  S=1/2
    opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 0, indices_, true ) );
    // S=1 (+) S=1/2  =>  S=3/2
    opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 1, indices_, true ) );
  }
  else{
    opReps_.push_back( build_compound_operator( true, -1, twoOps.at(0), oneOp.at(0), 0, indices_, true ) );
  }
  return false;
  }
  else{
//FIXME don't need to keep reconstructing whole array!
  indices_ = spinBlock_->get_op_array(RI_3INDEX).get_array().at(idx);
  int ix = indices_[0];
  int jx = indices_[1];
  int kx = indices_[2];
  if(!dmrginp.spinAdapted()){
    if(ix==jx) return true;
    if(ix==kx) return true;
  }

  // Get 2-index and 1-index ops as RI building blocks
  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_CRE).get_element(ix,jx);
  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(DES).get_element(kx);

  std::vector< boost::shared_ptr<SparseMatrix> > testOp = spinBlock_->get_op_array(DES_CRE).get_element(ix,jx);
  


  // Allocate and build operator representation on the fly as RI tensor product for each spin component
  opReps_.clear();
  if(dmrginp.spinAdapted()){
    // S=0 (+) S=1/2  =>  S=1/2
    opReps_.push_back( build_compound_operator( true, 1, twoOps.at(0), oneOp.at(0), 0, indices_, false) );
    // S=1 (+) S=1/2  =>  S=1/2
    opReps_.push_back( build_compound_operator( true, 1, twoOps.at(1), oneOp.at(0), 0, indices_, false) );
    // S=1 (+) S=1/2  =>  S=3/2
    opReps_.push_back( build_compound_operator( true, 1, twoOps.at(1), oneOp.at(0), 1, indices_, false) );
  }
  else{
    opReps_.push_back( build_compound_operator( true, 1, twoOps.at(0), oneOp.at(0), 0, indices_, false) );
  }
  return false;
  }
}

//===========================================================================================================================================================

Npdm_op_compound_CDD::Npdm_op_compound_CDD( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
//FIXMEMAW
//  size_ = 1;
  size_ = spinBlock_->get_op_array(RI_3INDEX).get_size();
//  if(!dmrginp.spinAdapted()) size_=1;
  is_local_ = true;
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((CD)D)";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_compound_CDD::set_local_ops( int idx )
{
  if(dmrginp.doimplicitTranspose()){
//pout << "getting compound CDD operator...\n";
//assert( idx == 0 );
  // Spatial orbital indices
//FIXME don't need to keep reconstructing whole array!
  indices_ = spinBlock_->get_op_array(RI_3INDEX).get_array().at(idx);
  int ix = indices_[0];
  int jx = indices_[1];
  int kx = indices_[2];
  if(!dmrginp.spinAdapted()) {
    if(kx==jx) return true;
    if(kx==ix) return true;
  }

  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_DES).get_element(ix,jx);
  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(CRE).get_element(kx);

////  // Spatial orbital indices
////  indices_.clear();
////  int ix, jx, kx;
////
////  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_DES).get_local_element(idx);
////  ix = twoOps.at(0)->get_orbs(0);
////  jx = twoOps.at(0)->get_orbs(1);
////  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(CRE).get_local_element(idx);
////  kx = oneOp.at(0)->get_orbs(0);
////  indices_.push_back( ix );
////  indices_.push_back( jx );
////  indices_.push_back( kx );
////pout << "CDD indices = " << ix << ", " << jx << ", " << kx << endl;

  opReps_.clear();
  if(dmrginp.spinAdapted()){
    // S=0 (+) S=1/2  =>  S=1/2
    opReps_.push_back( build_compound_operator( true, -1, twoOps.at(0), oneOp.at(0), 0, indices_, true ) );
    // S=1 (+) S=1/2  =>  S=1/2
    opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 0, indices_, true ) );
    // S=1 (+) S=1/2  =>  S=3/2
    opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 1, indices_, true ) );
  }
  else{
    opReps_.push_back( build_compound_operator( true, -1, twoOps.at(0), oneOp.at(0), 0, indices_, true ) );
    //opReps_.push_back( build_compound_operator( true, 1, twoOps.at(0), oneOp.at(0), 0, indices_, false ) );
  }


//-----------------------------
// 2 other ways
//-----------------------------

//  "(C(DD))";
//  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_CRE).get_local_element(idx);
//  build_pattern_ = "(C(DD))";
//  factor_ = -1.0;  //FIXME note we need -1 here!
//  opReps_.clear();
//  // S=0 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, oneOp.at(0), twoOps.at(0), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, oneOp.at(0), twoOps.at(1), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=3/2
//  opReps_.push_back( build_compound_operator( true, -1, oneOp.at(0), twoOps.at(1), 1, indices_, true ) );

//  "((CC)D)" + tranpose;
//  "(C(DD))";
//  transpose_ = true;
//  build_pattern_ = "(C(DD))"; //FIXME note build pattern reflects structure *after* transpose is taken
//  opReps_.clear();
//  // S=0 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(0), oneOp.at(0), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=3/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 1, indices_, true ) );

// PRINT
//pout << "2a CDD operator elements:\n";
//pout << *(opReps_[0]);
//pout << "2b CDD operator elements:\n";
//pout << *(opReps_[1]);
//pout << "4  CDD operator elements:\n";
//pout << *(opReps_[2]);
  return false;
  }
  else{
//FIXME don't need to keep reconstructing whole array!
  indices_ = spinBlock_->get_op_array(RI_3INDEX).get_array().at(idx);
  int ix = indices_[0];
  int jx = indices_[1];
  int kx = indices_[2];
  if(!dmrginp.spinAdapted()) {
    if(kx==jx) return true;
    if(kx==ix) return true;
  }

  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_DES).get_element(ix,jx);
  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(DES).get_element(kx);


  opReps_.clear();
  if(dmrginp.spinAdapted()){
    // S=0 (+) S=1/2  =>  S=1/2
    opReps_.push_back( build_compound_operator( true, 1, twoOps.at(0), oneOp.at(0), 0, indices_, false ) );
    // S=1 (+) S=1/2  =>  S=1/2
    opReps_.push_back( build_compound_operator( true, 1, twoOps.at(1), oneOp.at(0), 0, indices_, false ) );
    // S=1 (+) S=1/2  =>  S=3/2
    opReps_.push_back( build_compound_operator( true, 1, twoOps.at(1), oneOp.at(0), 1, indices_, false ) );
  }
  else{
    opReps_.push_back( build_compound_operator( true, 1, twoOps.at(0), oneOp.at(0), 0, indices_, false ) );
  }


//-----------------------------
// 2 other ways
//-----------------------------

//  "(C(DD))";
//  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_CRE).get_local_element(idx);
//  build_pattern_ = "(C(DD))";
//  factor_ = -1.0;  //FIXME note we need -1 here!
//  opReps_.clear();
//  // S=0 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, oneOp.at(0), twoOps.at(0), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, oneOp.at(0), twoOps.at(1), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=3/2
//  opReps_.push_back( build_compound_operator( true, -1, oneOp.at(0), twoOps.at(1), 1, indices_, true ) );

//  "((CC)D)" + tranpose;
//  "(C(DD))";
//  transpose_ = true;
//  build_pattern_ = "(C(DD))"; //FIXME note build pattern reflects structure *after* transpose is taken
//  opReps_.clear();
//  // S=0 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(0), oneOp.at(0), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=3/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 1, indices_, true ) );

// PRINT
//pout << "2a CDD operator elements:\n";
//pout << *(opReps_[0]);
//pout << "2b CDD operator elements:\n";
//pout << *(opReps_[1]);
//pout << "4  CDD operator elements:\n";
//pout << *(opReps_[2]);
  return false;
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

Npdm_op_compound_CDC::Npdm_op_compound_CDC( SpinBlock * spinBlock )
{
//abort();
  // Assume single site block
  // for spinAdpated dot block, creator must be on the left of destructor
  if(dmrginp.spinAdapted()) assert( spinBlock->size() == 1 && 1==0 );
  else assert(spinBlock->size() == 2);
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = 1;
  if(!dmrginp.spinAdapted()) size_ = spinBlock_->get_op_array(RI_3INDEX).get_size();
  is_local_ = true;
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((CD)C)";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_compound_CDC::set_local_ops( int idx )
{
//pout << "getting compound CDC operator...\n";
  // Spatial orbital indices
  if(!dmrginp.spinAdapted()){
    //pout << spinBlock_->get_op_array(RI_3INDEX).get_array().size();
    indices_ = spinBlock_->get_op_array(RI_3INDEX).get_array().at(idx);
    int ix = indices_[0];
    int jx = indices_[1];
    int kx = indices_[2];
    if(ix==kx) return true;
    if(jx==kx) return true;

    std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_DES).get_element(ix,jx);
    std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(CRE).get_element(kx);


    opReps_.clear();
    opReps_.push_back( build_compound_operator( true, 1, twoOps.at(0), oneOp.at(0), 0, indices_, false ) );
    return false;
  }
  else{

  indices_.clear();
  int ix, jx, kx;

  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_DES).get_local_element(idx);
  ix = twoOps.at(0)->get_orbs(0);
  jx = twoOps.at(0)->get_orbs(1);
  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(CRE).get_local_element(idx);
  kx = oneOp.at(0)->get_orbs(0);
//pout << "singlet CD operator elements:\n";
//pout << *(twoOps[0]);
//pout << "triplet CD operator elements:\n";
//pout << *(twoOps[1]);
//pout << "half C operator elements:\n";
//pout << *(oneOp[0]);

    // Assumed single site (i=j=k)
    assert ( ix == jx );
    assert ( jx == kx );
  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  if ( jx == kx ) {
    //FIXME I think this fails because of potential problems commuting operators with same indices in spin-transformation
    //pout << "WARNING: skipping this operator\n";
    return true;
  }
//pout << "indices  " << ix << " " << jx << " " << kx << std::endl;

//----------
// 1st way
//----------

  opReps_.clear();
  if(dmrginp.spinAdapted()){
    // S=0 (+) S=1/2  =>  S=1/2
    opReps_.push_back( build_compound_operator( true, 1, twoOps.at(0), oneOp.at(0), 0, indices_, false ) );
    // S=1 (+) S=1/2  =>  S=1/2
    opReps_.push_back( build_compound_operator( true, 1, twoOps.at(1), oneOp.at(0), 0, indices_, false ) );
    // S=1 (+) S=1/2  =>  S=3/2
    opReps_.push_back( build_compound_operator( true, 1, twoOps.at(1), oneOp.at(0), 1, indices_, false ) );
  }
  else{
    opReps_.push_back( build_compound_operator( true, 1, twoOps.at(0), oneOp.at(0), 0, indices_, false ) );
  }

//----------
// 2nd way
//----------
//
//  build_pattern_ = "(C(DC))";
//  factor_ = 1.0;
//  opReps_.clear();
//  //FIXME what should second argument be?
//  // S=0 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, oneOp.at(0), twoOps.at(0), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, oneOp.at(0), twoOps.at(1), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=3/2
//  opReps_.push_back( build_compound_operator( true, -1, oneOp.at(0), twoOps.at(1), 1, indices_, true ) );

//----------
// 3rd way
//----------

//  // "((CD)D)" + transpose
//  build_pattern_ = "(C(DC))";
//  transpose_ = true;
//  factor_ = 1.0;
//  opReps_.clear();
//  // S=0 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(0), oneOp.at(0), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=3/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 1, indices_, true ) );

//---------------------
// 1st way re-written
//---------------------
//
//  build_pattern_ = "((CD)C)";
//  factor_ = 1.0;
//  opReps_.clear();
//
//  // Singlet or triplet CD operator
//  boost::shared_ptr<SparseMatrix> cdOp (new CreDes);
//  cdOp->set_orbs() = indices_;
//  cdOp->set_orbs().pop_back();
//  cdOp->set_initialised() = true;
//  cdOp->set_fermion() = false;
//
//  // CD s=0
//  cdOp->set_deltaQuantum() = ( oneOp.at(0)->get_deltaQuantum() - oneOp.at(0)->get_deltaQuantum() ).at(0);
//  cdOp->allocate( spinBlock_->get_stateInfo() );
//  operatorfunctions::Product(spinBlock_, *oneOp.at(0), Transposeview(*oneOp.at(0)), *cdOp, 1.0 );
//pout << "singlet CD operator elements (built by me):\n";
//pout << *cdOp;
//  // S=0 (+) S=1/2  =>  S=1/2
//  boost::shared_ptr<SparseMatrix> cdcOp (new Cre);
//  cdcOp->set_orbs() = indices_;
//  cdcOp->set_initialised() = true;
//  cdcOp->set_fermion() = true;
//  cdcOp->set_deltaQuantum() = ( cdOp->get_deltaQuantum() + oneOp.at(0)->get_deltaQuantum() ).at(0);
//  cdcOp->allocate( spinBlock_->get_stateInfo() );
//  operatorfunctions::Product(spinBlock_, *cdOp, *oneOp.at(0), *cdcOp, 1.0 );
//  opReps_.push_back( cdcOp );
////pout << *cdcOp;
//
//  // CD s=1
//  cdOp->set_deltaQuantum() = ( oneOp.at(0)->get_deltaQuantum() - oneOp.at(0)->get_deltaQuantum() ).at(1);
//  cdOp->allocate( spinBlock_->get_stateInfo() );
//  operatorfunctions::Product(spinBlock_, *oneOp.at(0), Transposeview(*oneOp.at(0)), *cdOp, 1.0 );
//pout << "triplet CD operator elements (built by me):\n";
//pout << *cdOp;
//  // S=1 (+) S=1/2  =>  S=1/2
//  boost::shared_ptr<SparseMatrix> cdcOp2 (new Cre);
//  cdcOp2->set_orbs() = indices_;
//  cdcOp2->set_initialised() = true;
//  cdcOp2->set_fermion() = true;
//  cdcOp2->set_deltaQuantum() = ( cdOp->get_deltaQuantum() + oneOp.at(0)->get_deltaQuantum() ).at(0);
//  cdcOp2->allocate( spinBlock_->get_stateInfo() );
//  operatorfunctions::Product(spinBlock_, *cdOp, *oneOp.at(0), *cdcOp2, 1.0 );
//  opReps_.push_back( cdcOp2 );
//  // S=1 (+) S=1/2  =>  S=3/2
//  boost::shared_ptr<SparseMatrix> cdcOp3 (new Cre);
//  cdcOp3->set_orbs() = indices_;
//  cdcOp3->set_initialised() = true;
//  cdcOp3->set_fermion() = true;
//  cdcOp3->set_deltaQuantum() = ( cdOp->get_deltaQuantum() + oneOp.at(0)->get_deltaQuantum() ).at(1);
//  cdcOp3->allocate( spinBlock_->get_stateInfo() );
//  operatorfunctions::Product(spinBlock_, *cdOp, *oneOp.at(0), *cdcOp3, 1.0 );
//  opReps_.push_back( cdcOp3 );

// PRINT
//pout << "CDC operator elements (built by me):\n";
//pout << "2a CDC operator elements:\n";
//pout << *(opReps_[0]);
//pout << "2b CDC operator elements:\n";
//pout << *(opReps_[1]);
//pout << "4  CDC operator elements:\n";
//pout << *(opReps_[2]);

  return false;
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

Npdm_op_compound_CCC::Npdm_op_compound_CCC( SpinBlock * spinBlock )
{
  // Assume single site block
abort(); // << this operator should always be zero on one site!
  assert( spinBlock->size() == 1 );
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = 1;
  is_local_ = true;
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((CC)C)";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_compound_CCC::set_local_ops( int idx )
{
//pout << "getting compound CCC operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_CRE).get_local_element(idx);
  ix = twoOps.at(0)->get_orbs(0);
  jx = twoOps.at(0)->get_orbs(1);
  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(CRE).get_local_element(idx);
  kx = oneOp.at(0)->get_orbs(0);

  if(dmrginp.spinAdapted()) return true;
  // Assumed single site (i=j=k)
  assert ( ix == jx );
  assert ( jx == kx );
  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );

  opReps_.clear();
  // S=0 (+) S=1/2  =>  S=1/2
  opReps_.push_back( build_compound_operator( true, 1, twoOps.at(0), oneOp.at(0), 0, indices_, false ) );
  // S=1 (+) S=1/2  =>  S=1/2
  opReps_.push_back( build_compound_operator( true, 1, twoOps.at(1), oneOp.at(0), 0, indices_, false ) );
  // S=1 (+) S=1/2  =>  S=3/2
  opReps_.push_back( build_compound_operator( true, 1, twoOps.at(1), oneOp.at(0), 1, indices_, false ) );

  return false;
}

//===========================================================================================================================================================

Npdm_op_compound_DCD::Npdm_op_compound_DCD( SpinBlock * spinBlock )
{
//abort();
  // Assume single site block
  // for spinAdpated dot block, creator must be on the left of destructor
  if(dmrginp.spinAdapted()) assert( spinBlock->size() == 1 && 1==0 );
  else assert(spinBlock->size()==2);
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  if(dmrginp.spinAdapted()) size_ = 1;
  else size_= spinBlock_->get_op_array(RI_3INDEX).get_size();
  is_local_ = true;
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((DC)D)";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_compound_DCD::set_local_ops( int idx )
{
  if(dmrginp.doimplicitTranspose()){
  if(!dmrginp.spinAdapted()){
    indices_ = spinBlock_->get_op_array(RI_3INDEX).get_array().at(idx);
    int ix = indices_[0];
    int jx = indices_[1];
    int kx = indices_[2];
    if(ix==jx) return true;
    if(ix==kx) return true;

    std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(CRE).get_element(ix);
    boost::shared_ptr<SparseMatrix>  t_oneOp(new Transposeview(*oneOp.at(0)));
    //pout << Transposeview(*oneOp.at(0));
    //*t_oneOp = Transposeview(*oneOp.at(0));
    std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(CRE_DES).get_element(jx,kx);



    opReps_.clear();
    opReps_.push_back( build_compound_operator( true, 1, t_oneOp, twoOps.at(0), 0, indices_, false ) );
    return false;
  }
  pout << " no need to build DCD in a dot site for spin Adapted\n";
  abort();
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(DES_CRE).get_local_element(idx);
  ix = twoOps.at(0)->get_orbs(0);
  jx = twoOps.at(0)->get_orbs(1);
  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(CRE).get_local_element(idx);
  kx = oneOp.at(0)->get_orbs(0);

  if(ix==jx) return true;
  if(ix==kx) return true;
  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  if ( ix == jx ) {
    //FIXME I think this fails because of potential problems commuting operators with same indices in spin-transformation
    //pout << "WARNING: skipping this operator\n";
    return true;
  }

  opReps_.clear();
  if(dmrginp.spinAdapted()){
    // S=0 (+) S=1/2  =>  S=1/2
    opReps_.push_back( build_compound_operator( true, -1, twoOps.at(0), oneOp.at(0), 0, indices_, true ) );
    // S=1 (+) S=1/2  =>  S=1/2
    opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 0, indices_, true ) );
    // S=1 (+) S=1/2  =>  S=3/2
    opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 1, indices_, true ) );
  }
  else{
    opReps_.push_back( build_compound_operator( true, -1, twoOps.at(0), oneOp.at(0), 0, indices_, true ) );
  }

  return false;
  }
  else{
  if(!dmrginp.spinAdapted()){
    indices_ = spinBlock_->get_op_array(RI_3INDEX).get_array().at(idx);
    int ix = indices_[0];
    int jx = indices_[1];
    int kx = indices_[2];
    if(ix==jx) return true;
    if(ix==kx) return true;

    std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(DES_CRE).get_element(ix,jx);
    std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(DES).get_element(kx);


    opReps_.clear();
    opReps_.push_back( build_compound_operator( true, 1, twoOps.at(0), oneOp.at(0), 0, indices_, false ) );
    return false;
  }
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(DES_CRE).get_local_element(idx);
  ix = twoOps.at(0)->get_orbs(0);
  jx = twoOps.at(0)->get_orbs(1);
  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(DES).get_local_element(idx);
  kx = oneOp.at(0)->get_orbs(0);

  if(jx==ix) return true;
  if(kx==ix) return true;
  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  if ( ix == jx ) {
    //FIXME I think this fails because of potential problems commuting operators with same indices in spin-transformation
    //pout << "WARNING: skipping this operator\n";
    return true;
  }

  opReps_.clear();
  // S=0 (+) S=1/2  =>  S=1/2
  opReps_.push_back( build_compound_operator( true, 1, twoOps.at(0), oneOp.at(0), 0, indices_, false ) );
  // S=1 (+) S=1/2  =>  S=1/2
  opReps_.push_back( build_compound_operator( true, 1, twoOps.at(1), oneOp.at(0), 0, indices_, false ) );
  // S=1 (+) S=1/2  =>  S=3/2
  opReps_.push_back( build_compound_operator( true, 1, twoOps.at(1), oneOp.at(0), 1, indices_, false ) );
  return false;
  }

}

//===========================================================================================================================================================

Npdm_op_compound_DDC::Npdm_op_compound_DDC( SpinBlock * spinBlock )
{
  // Assume single site block
abort();
//  assert( spinBlock->size() == 1 );
//  opReps_.clear();
//  indices_.clear();
//  spinBlock_ = spinBlock;
//  size_ = 1;
//  is_local_ = true;
//  factor_ = 1.0;
//  transpose_ = false;
//  build_pattern_ = "((DC)D)";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_compound_DDC::set_local_ops( int idx )
{
//pout << "getting compound DDC operator...\n";
  // Spatial orbital indices
  // FIXME, really no DDC ? 
abort();
//  indices_.clear();
//  int ix, jx, kx;
//
//  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(DES_CRE).get_local_element(idx);
//  ix = twoOps.at(0)->get_orbs(0);
//  jx = twoOps.at(0)->get_orbs(1);
//  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(CRE).get_local_element(idx);
//  kx = oneOp.at(0)->get_orbs(0);
//
//  // Assumed single site (i=j=k)
//  assert ( ix == jx );
//  assert ( jx == kx );
//  indices_.push_back( ix );
//  indices_.push_back( jx );
//  indices_.push_back( kx );
//  if ( ix == jx ) {
//    //FIXME I think this fails because of potential problems commuting operators with same indices in spin-transformation
//    //pout << "WARNING: skipping this operator\n";
//    return true;
//  }
//
//  opReps_.clear();
//  // S=0 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(0), oneOp.at(0), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=1/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 0, indices_, true ) );
//  // S=1 (+) S=1/2  =>  S=3/2
//  opReps_.push_back( build_compound_operator( true, -1, twoOps.at(1), oneOp.at(0), 1, indices_, true ) );
//
//  return false;
//
}

//===========================================================================================================================================================

Npdm_op_compound_DCC::Npdm_op_compound_DCC( SpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(RI_3INDEX).get_size();
  if(!dmrginp.spinAdapted()) size_=1;
  is_local_ = true;
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "((DC)C)";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_compound_DCC::set_local_ops( int idx )
{
//pout << "getting compound DCC operator...\n";
  // Spatial orbital indices
//FIXME don't need to keep reconstructing whole array!
//
  pout << " no need to build DCD in a dot site for spin Adapted\n";
  abort();
  indices_ = spinBlock_->get_op_array(RI_3INDEX).get_array().at(idx);
  int ix = indices_[0];
  int jx = indices_[1];
  int kx = indices_[2];

  if(!dmrginp.spinAdapted()){
    if(ix==jx) return true;
    if(kx==jx) return true;
    if(kx==ix) return true;
  }

  // Skip operator if i,j equal
  if ( ix == jx ) return true;

  // Get 2-index and 1-index ops as RI building blocks
  std::vector< boost::shared_ptr<SparseMatrix> > twoOps = spinBlock_->get_op_array(DES_CRE).get_element(ix,jx);
  std::vector< boost::shared_ptr<SparseMatrix> > oneOp = spinBlock_->get_op_array(CRE).get_element(kx);

  // Allocate and build operator representation on the fly as RI tensor product for each spin component
  opReps_.clear();
  if(dmrginp.spinAdapted()){
    // S=0 (+) S=1/2  =>  S=1/2
    opReps_.push_back( build_compound_operator( true, 1, twoOps.at(0), oneOp.at(0), 0, indices_, false ) );
    // S=1 (+) S=1/2  =>  S=1/2
    opReps_.push_back( build_compound_operator( true, 1, twoOps.at(1), oneOp.at(0), 0, indices_, false ) );
    // S=1 (+) S=1/2  =>  S=3/2
    opReps_.push_back( build_compound_operator( true, 1, twoOps.at(1), oneOp.at(0), 1, indices_, false ) );
  }
  else{
    opReps_.push_back( build_compound_operator( true, 1, twoOps.at(0), oneOp.at(0), 0, indices_, false ) );

  }

  return false;
}

//===========================================================================================================================================================

}
}

