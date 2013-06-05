/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//
// Extension of Operators.C for 3-index operators
//
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "op_components.h"
#include "BaseOperator.h"
#include "npdm_operators.h"
#include "spinblock.h"
#include "operatorfunctions.h"
#include "tensor_operator.h"

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Cre,Des)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinAdapted::CreCreDes::build(const SpinBlock& b)
{
pout << "building CreCreDes renormalized operator...\n";
//build_in_csf_space(b); //fails

  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_stateInfo());
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
pout << "indices  " << i << " " << j << " " << k << std::endl;

  SpinBlock* sysBlock = b.get_leftBlock();
  SpinBlock* dotBlock = b.get_rightBlock();
pout << "sysBlock:\n";
pout << *sysBlock;
pout << "dotBlock:\n";
pout << *dotBlock;

//  SpinQuantum deltaQuantum12 = get_quantum_ladder().at("0").at(0);
//  assert( get_quantum_ladder().size() == 2 );

  // Sys has i,j,k
  if (sysBlock->get_op_array(CRE_CRE_DES).has(i,j,k)) {
pout << "maw sys(i,j,k)\n";
    std::string build_pattern_old = sysBlock->get_op_array(CRE_CRE_DES).get_element(i,j,k).at(0)->get_build_pattern();
    assert( build_pattern_old == sysBlock->get_op_array(CRE_CRE_DES).get_element(i,j,k).at(1)->get_build_pattern() );
    assert( build_pattern_old == sysBlock->get_op_array(CRE_CRE_DES).get_element(i,j,k).at(2)->get_build_pattern() );
    const boost::shared_ptr<SparseMatrix>& op = sysBlock->get_op_rep(CRE_CRE_DES, quantum_ladder.at(build_pattern_old), i,j,k);
    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *op, &b, &(b.get_stateInfo()), *this);
    build_pattern = build_pattern_old;
  }
  // Dot has i,j,k
  else if (dotBlock->get_op_array(CRE_CRE_DES).has(i,j,k)) {
pout << "maw dot(i,j,k)\n";
    std::string build_pattern_old = dotBlock->get_op_array(CRE_CRE_DES).get_element(i,j,k).at(0)->get_build_pattern();
    assert( build_pattern_old == dotBlock->get_op_array(CRE_CRE_DES).get_element(i,j,k).at(1)->get_build_pattern() );
pout <<  dotBlock->get_op_array(CRE_CRE_DES).get_element(i,j,k).at(2)->get_build_pattern()  << std::endl;
    assert( build_pattern_old == dotBlock->get_op_array(CRE_CRE_DES).get_element(i,j,k).at(2)->get_build_pattern() );
    const boost::shared_ptr<SparseMatrix>& op = dotBlock->get_op_rep(CRE_CRE_DES, quantum_ladder.at(build_pattern_old), i,j,k);
    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *op, &b, &(b.get_stateInfo()), *this);
    build_pattern = build_pattern_old;
  }
  else assert(false);

  dmrginp.makeopsT -> stop();
pout << "CCD op\n";
pout << *this;
pout << "done building CreCreDes renormalized operator!\n";

}

////
////
////  // Sys has i,j only
////  else if (sysBlock->get_op_array(CRE_CRE).has(i,j)) {
//////pout << "maw sys(i,j)\n";
////    assert( dotBlock->get_op_array(CRE).has(k) );
////    const boost::shared_ptr<SparseMatrix>& opCC = sysBlock->get_op_rep(CRE_CRE, deltaQuantum12, i,j);
////    const boost::shared_ptr<SparseMatrix>& opC3 = dotBlock->get_op_rep(CRE, getSpinQuantum(k), k);
////    // Build ((CC)D)
////    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, *opCC, Transposeview(*opC3), &b, &(b.get_stateInfo()), *this, 1.0);
////  }
/////////////////////////////////////////////////
////  // Sys has j,k only
////  else if (sysBlock->get_op_array(CRE_DES).has(j,k)) {
//////pout << "maw sys(j,k)\n";
////////    assert( dotBlock->get_op_array(CRE).has(i) );
////    const boost::shared_ptr<SparseMatrix>& opC = dotBlock->get_op_rep(CRE, getSpinQuantum(i), i);
////    const boost::shared_ptr<SparseMatrix>& opCD0 = sysBlock->get_op_array(CRE_DES).get_element(j,k).at(0)->getworkingrepresentation(sysBlock);
//////pout << "maw opC\n";
//////pout << *opC;
//////pout << "maw opCD[0]\n";
//////pout << *opCD0;
////    const boost::shared_ptr<SparseMatrix>& opCD1 = sysBlock->get_op_array(CRE_DES).get_element(j,k).at(1)->getworkingrepresentation(sysBlock);
//////pout << "maw opCD[1]\n";
//////pout << *opCD1;
////    // Build (C(CD))
////    double parity = getCommuteParity( opC->get_deltaQuantum(), opCD0->get_deltaQuantum(), get_deltaQuantum() );
////    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *opC, *opCD0, &b, &(b.get_stateInfo()), *this, 1.0*parity);
/////////////////////////////
////
//////    const boost::shared_ptr<SparseMatrix>& opCC = b.get_op_rep(CRE_CRE, deltaQuantum12, i, j);
//////pout << "maw opCC\n";
//////pout << *opCC;
////
////
////
/////////////////////////////
////    const boost::shared_ptr<SparseMatrix>& opC1 = dotBlock->get_op_array(CRE).get_element(i).at(0)->getworkingrepresentation(dotBlock);
////    const boost::shared_ptr<SparseMatrix>& opC2 = sysBlock->get_op_array(CRE).get_element(j).at(0)->getworkingrepresentation(sysBlock);
////    const boost::shared_ptr<SparseMatrix>& opC3 = sysBlock->get_op_array(CRE).get_element(k).at(0)->getworkingrepresentation(sysBlock);
////pout << "maw opC1\n";
////pout << *opC1;
////pout << "maw opC2\n";
////pout << *opC2;
////pout << "maw opC3\n";
////pout << *opC3;
////////
/////////////////////////////
////////    assert( dotBlock->get_op_array(CRE).has(i) );
////////    assert( sysBlock->get_op_array(CRE).has(j) );
////////    assert( sysBlock->get_op_array(CRE).has(k) );
////////    // Build (CiCj) first
////////    CreCre opCC;
////////    opCC.set_orbs() = { i, j };
////////    opCC.set_initialised() = true;
////////    opCC.set_fermion() = false;
////////    opCC.set_deltaQuantum() = deltaQuantum12;
////////    opCC.allocate( b.get_stateInfo() );
////////    const boost::shared_ptr<SparseMatrix>& opC1 = dotBlock->get_op_rep(CRE, getSpinQuantum(i), i);
////////    const boost::shared_ptr<SparseMatrix>& opC2 = sysBlock->get_op_rep(CRE, getSpinQuantum(j), j);
////////
////////    parity = getCommuteParity( opC1->get_deltaQuantum(), opC2->get_deltaQuantum(), opCC.get_deltaQuantum() );
////////    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *opC1, *opC2, &b, &(b.get_stateInfo()), opCC, 1.0*parity);
////////    // Build C3 on full spinblock
////////    const boost::shared_ptr<SparseMatrix>& opC3 = sysBlock->get_op_rep(CRE, getSpinQuantum(k), k);
////////pout << "maw opC1\n";
////////pout << *opC1;
////////pout << "maw opC2\n";
////////pout << *opC2;
////////pout << "maw opC3\n";
////////pout << *opC3;
////////pout << "maw opCC(i,j)\n";
////////pout << opCC;
////////    CreCre opC3big;
////////    opC3big.set_orbs() = { k };
////////    opC3big.set_initialised() = true;
////////    opC3big.set_fermion() = true;
////////    opC3big.set_deltaQuantum() = getSpinQuantum(k);
////////    opC3big.allocate( b.get_stateInfo() );
////////    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *opC3, &b, &(b.get_stateInfo()), opC3big);
////////pout << "maw opC3big(k)\n";
////////pout << opC3big;
////////    // Build ((CC)D) on full spinblock
////////    SpinAdapted::operatorfunctions::Product(&b, opCC, Transposeview(opC3big), *this, 1.0 );
////
////
////
////
////
////  }
////  // Sys has i only
////  else if (sysBlock->get_op_array(CRE).has(i)) {
//////pout << "maw sys(i)\n";
////    assert( j == k );
////    assert( dotBlock->get_op_array(CRE_DES).has(j,k) );
////    assert( dotBlock->get_op_array(CRE).has(j) );
////    // Build (CiCj) first
////    CreCre opCC;
////    opCC.set_orbs() = { i, j };
////    opCC.set_initialised() = true;
////    opCC.set_fermion() = false;
////    opCC.set_deltaQuantum() = deltaQuantum12;
////    opCC.allocate( b.get_stateInfo() );
////    const boost::shared_ptr<SparseMatrix>& opC1 = sysBlock->get_op_rep(CRE, getSpinQuantum(i), i);
////    const boost::shared_ptr<SparseMatrix>& opC2 = dotBlock->get_op_rep(CRE, getSpinQuantum(j), j);
////    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, *opC1, *opC2, &b, &(b.get_stateInfo()), opCC, 1.0);
////    // Build C3 on full spinblock
////    CreCre opC3big;
////    opC3big.set_orbs() = { k };
////    opC3big.set_initialised() = true;
////    opC3big.set_fermion() = true;
////    opC3big.set_deltaQuantum() = getSpinQuantum(k);
////    opC3big.allocate( b.get_stateInfo() );
////    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *opC2, &b, &(b.get_stateInfo()), opC3big);
////    // Build ((CC)D) on full spinblock
////    SpinAdapted::operatorfunctions::Product(&b, opCC, Transposeview(opC3big), *this, 1.0 );
////  }
////  // Sys has k only
////  else if (sysBlock->get_op_array(CRE).has(k)) {
//////pout << "maw sys(k)\n";
////    assert( i == j );
////    assert( dotBlock->get_op_array(CRE_DES).has(i,j) );
////    const boost::shared_ptr<SparseMatrix>& opCC = dotBlock->get_op_rep(CRE_CRE, deltaQuantum12, i,j);
////    Transposeview opD = Transposeview( sysBlock->get_op_rep(CRE, getSpinQuantum(k), k) );
////    // Build ((CC)D)
////    double parity = getCommuteParity( opCC->get_deltaQuantum(), opD.get_deltaQuantum(), get_deltaQuantum() );
////    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *opCC, opD, &b, &(b.get_stateInfo()), *this, 1.0*parity);
////  }
////  // Dot has i,j,k
////  else if (dotBlock->get_op_array(CRE_CRE_DES).has(i,j,k)) {
////pout << "maw dot(i,j,k)\n";
////    build_pattern = op->get_build_pattern();
////    const boost::shared_ptr<SparseMatrix>& op = dotBlock->get_op_rep(CRE_CRE_DES, quantum_ladder, i,j,k);
////    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *op, &b, &(b.get_stateInfo()), *this);
////  }
////  else assert(false);
////
//////pout << "CCD op\n";
//////pout << *this;
//////pout << "done building CreCreDes renormalized operator!\n";

////-------------------------------------------------------------------------------------------------
/// BELOW IS WRONG!

////
////  //---------------------
////  //  Forwards sweep
////  //---------------------
////  // (xxx,0)
////  if (sysBlock->get_op_array(CRE_CRE_DES).has(i,j,k)) {
////pout << "maw forward(xxx,0)\n";
////    const boost::shared_ptr<SparseMatrix>& op = sysBlock->get_op_rep(CRE_CRE_DES, quantum_ladder, i,j,k);
////pout << "opCCD\n";
////pout << *op;
////    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *op, &b, &(b.get_stateInfo()), *this);
////  }
////  // (XX,X)
////  else if (sysBlock->get_op_array(CRE_CRE).has(i,j)) {
////pout << "maw forward(XX,X)\n";
////    assert( sysBlock->get_op_array(CRE).has(i) );
////    assert( sysBlock->get_op_array(CRE).has(j) );
////    assert( dotBlock->get_op_array(CRE).has(k) );
////    if ( i == j ) {
////      const boost::shared_ptr<SparseMatrix>& opCCref = sysBlock->get_op_rep(CRE_CRE, deltaQuantum12, i,j);
//////pout << "opCCref\n";
//////pout << *opCCref;
//////      const boost::shared_ptr<SparseMatrix>& opC3 = dotBlock->get_op_rep(CRE, getSpinQuantum(k), k);
//////      // Build ((CC)D) as tensor product
//////      SpinAdapted::operatorfunctions::TensorProduct(sysBlock, *opCC, Transposeview(*opC3), &b, &(b.get_stateInfo()), *this, 1.0);
////    }
//////    else {
////    {
////      const boost::shared_ptr<SparseMatrix>& opC1 = sysBlock->get_op_rep(CRE, getSpinQuantum(i), i);
////      const boost::shared_ptr<SparseMatrix>& opC2 = sysBlock->get_op_rep(CRE, getSpinQuantum(j), j);
////      const boost::shared_ptr<SparseMatrix>& opC3 = dotBlock->get_op_rep(CRE, getSpinQuantum(k), k);
////      // Build CC first as product on sysBlock
////      CreCre opCC;
////      opCC.set_orbs() = { i, j };
////      opCC.set_initialised() = true;
////      opCC.set_fermion() = false;
////      opCC.set_deltaQuantum() = deltaQuantum12;
////      opCC.allocate( sysBlock->get_stateInfo() );
//////pout << "opC1\n";
//////pout << *opC1;
//////pout << "opC2\n";
//////pout << *opC2;
//////pout << "opC3\n";
//////pout << *opC3;
////      SpinAdapted::operatorfunctions::Product(sysBlock, *opC1, *opC2, opCC, 1.0 );
//////pout << "opCC\n";
//////pout << opCC;
////      // Build ((CC)D) as tensor product
////      SpinAdapted::operatorfunctions::TensorProduct(sysBlock, opCC, Transposeview(*opC3), &b, &(b.get_stateInfo()), *this, 1.0);
//////pout << "opCCD\n";
//////pout << *this;
////    }
////  }
////  // (X,XX)
////  else if (sysBlock->get_op_array(CRE).has(i)) {
////pout << "maw forward(X,XX)\n";
////    assert( j == k );
////    assert( sysBlock->get_op_array(CRE).has(i) );
////    assert( dotBlock->get_op_array(CRE).has(j) );
////    const boost::shared_ptr<SparseMatrix>& opC1 = sysBlock->get_op_rep(CRE, getSpinQuantum(i), i);
////    const boost::shared_ptr<SparseMatrix>& opC2 = dotBlock->get_op_rep(CRE, getSpinQuantum(j), j);
////    // Build CC first as tensor product
////    CreCre opCC;
////    opCC.set_orbs() = { i, j };
////    opCC.set_initialised() = true;
////    opCC.set_fermion() = false;
////    opCC.set_deltaQuantum() = deltaQuantum12;
////    opCC.allocate( b.get_stateInfo() );
////    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, *opC1, *opC2, &b, &(b.get_stateInfo()), opCC, 1.0);
//////pout << "opCC\n";
//////pout << opCC;
////    // Get C3 operator
////    Cre opC3;
////    opC3.set_orbs() = { k };
////    opC3.set_initialised() = true;
////    opC3.set_fermion() = true;
////    opC3.set_deltaQuantum() = getSpinQuantum(k);
////    opC3.allocate( b.get_stateInfo() );
////    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *opC2, &b, &(b.get_stateInfo()), opC3);
//////pout << "opD3\n";
//////pout << opD3;
////    // Build ((CC)D) in space of full spinblock
////    SpinAdapted::operatorfunctions::Product(&b, opCC, Transposeview(opC3), *this, 1.0 );
//////pout << "opCCD\n";
//////pout << *this;
////  }
////  //---------------------
////  //  Backwards sweep
////  //---------------------
////  // (0,xxx)
////  else if (dotBlock->get_op_array(CRE_CRE_DES).has(i,j,k)) {
////pout << "maw backward(0,xxx)\n";
////    const boost::shared_ptr<SparseMatrix>& op = dotBlock->get_op_rep(CRE_CRE_DES, quantum_ladder, i,j,k);
////    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *op, &b, &(b.get_stateInfo()), *this);
////  }
////  // (X,XX)
////  else if (sysBlock->get_op_array(CRE_CRE).has(j,k)) {
////pout << "maw backward(X,XX)\n";
////    assert( dotBlock->get_op_array(CRE).has(i) );
////    assert( sysBlock->get_op_array(CRE).has(j) );
////    assert( sysBlock->get_op_array(CRE).has(k) );
////    const boost::shared_ptr<SparseMatrix>& opC1 = dotBlock->get_op_rep(CRE, getSpinQuantum(i), i);
////    const boost::shared_ptr<SparseMatrix>& opC2 = sysBlock->get_op_rep(CRE, getSpinQuantum(j), j);
////    const boost::shared_ptr<SparseMatrix>& opC3 = sysBlock->get_op_rep(CRE, getSpinQuantum(k), k);
//////pout << "opC1\n";
//////pout << *opC1;
////    // Build CC first as product on sysBlock
////    CreCre opCC;
////    opCC.set_orbs() = { i, j };
////    opCC.set_initialised() = true;
////    opCC.set_fermion() = false;
////    opCC.set_deltaQuantum() = deltaQuantum12;
////    opCC.allocate( b.get_stateInfo() );
////    double parity = getCommuteParity( opC1->get_deltaQuantum(), opC2->get_deltaQuantum(), opCC.get_deltaQuantum() );
////    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *opC1, *opC2, &b, &(b.get_stateInfo()), opCC, 1.0*parity);
//////pout << "opCC\n";
//////pout << opCC;
////    // Get C3 operator
////    Cre opC3big;
////    opC3big.set_orbs() = { k };
////    opC3big.set_initialised() = true;
////    opC3big.set_fermion() = true;
////    opC3big.set_deltaQuantum() = getSpinQuantum(k);
////    opC3big.allocate( b.get_stateInfo() );
//////FIXME can we put transpose here to make opD first?
////    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *opC3, &b, &(b.get_stateInfo()), opC3big);
//////pout << "opC3big\n";
//////pout << opC3big;
////    // Build ((CC)D) in space of full spinblock
////    SpinAdapted::operatorfunctions::Product(&b, opCC, Transposeview(opC3big), *this, 1.0 );
////  }
////  // (XX,X)
////  else if (sysBlock->get_op_array(CRE).has(k)) {
////pout << "maw backward(XX,X)\n";
////    assert( i == j );
////    assert( dotBlock->get_op_array(CRE).has(i) );
////    assert( sysBlock->get_op_array(CRE).has(k) );
////    const boost::shared_ptr<SparseMatrix>& opCC = dotBlock->get_op_rep(CRE_CRE, deltaQuantum12, i,j);
////    const boost::shared_ptr<SparseMatrix>& opC3 = sysBlock->get_op_rep(CRE, getSpinQuantum(k), k);
////    // Build ((CC)D) as tensor product
////    double parity = getCommuteParity( opCC->get_deltaQuantum(), opC3->get_deltaQuantum(), get_deltaQuantum() );
////    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *opCC, Transposeview(opC3), &b, &(b.get_stateInfo()), *this, 1.0*parity);
////  }
////  else assert(false);
////
//////pout << "CCD op\n";
//////pout << *this;
////  dmrginp.makeopsT -> stop();
////
////}

////-------------------------------------------------------------------------------------------------
/// BELOW IS WRONG!
///
///
///  //--------------------------
///  // Get as product of ops
///  //--------------------------
///
///  // Get predefined 2-index op and expand it
/////  if (leftBlock->get_op_array(CRE_CRE).has(i,j)) {
/////    const boost::shared_ptr<SparseMatrix> opCCref = leftBlock->get_op_rep(CRE_CRE, deltaQuantum12, i,j);
/////pout << "opCCref\n";
/////pout << *opCCref;
/////    CreCre opCCrefb;
/////    opCCrefb.set_orbs() = get_orbs();
/////    opCCrefb.set_initialised() = true;
/////    opCCrefb.set_fermion() = false;
/////    opCCrefb.set_deltaQuantum() = deltaQuantum12;
/////    opCCrefb.allocate( b.get_stateInfo() );
/////    // Expand small op to big_op
/////    SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *opCCref, &b, &(b.get_stateInfo()), opCCrefb);
/////pout << "opCCref_expanded\n";
/////pout << opCCrefb;
/////  }
///
///  // Get i-op
///  Cre opC1;
///  opC1.set_orbs() = { i };
///  opC1.set_initialised() = true;
///  opC1.set_fermion() = true;
///  opC1.set_deltaQuantum() = getSpinQuantum(i);
///  opC1.allocate( b.get_stateInfo() );
///  // Expand small op to big_op
///  if (leftBlock->get_op_array(CRE).has(i)) {
///    const boost::shared_ptr<SparseMatrix> opC = leftBlock->get_op_rep(CRE, getSpinQuantum(i), i);
///    SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *opC, &b, &(b.get_stateInfo()), opC1);
///  }
///  else if (rightBlock->get_op_array(CRE).has(i)) {
///    const boost::shared_ptr<SparseMatrix> opC = rightBlock->get_op_rep(CRE, getSpinQuantum(i), i);
///    SpinAdapted::operatorfunctions::TensorTrace(rightBlock, *opC, &b, &(b.get_stateInfo()), opC1);
///  }
///  else assert(false);
///
///  // Get j-op
///  Cre opC2;
///  opC2.set_orbs() = { j };
///  opC2.set_initialised() = true;
///  opC2.set_fermion() = true;
///  opC2.set_deltaQuantum() = getSpinQuantum(j);
///  opC2.allocate( b.get_stateInfo() );
///  // Expand small op to big_op
///  if (leftBlock->get_op_array(CRE).has(j)) {
///    const boost::shared_ptr<SparseMatrix> opC = leftBlock->get_op_rep(CRE, getSpinQuantum(j), j);
///    SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *opC, &b, &(b.get_stateInfo()), opC2);
///  }
///  else if (rightBlock->get_op_array(CRE).has(j)) {
///    const boost::shared_ptr<SparseMatrix> opC = rightBlock->get_op_rep(CRE, getSpinQuantum(j), j);
///    SpinAdapted::operatorfunctions::TensorTrace(rightBlock, *opC, &b, &(b.get_stateInfo()), opC2);
///  }
///  else assert(false);
///
///  // Get k-op
///  Cre opD3;
///  opD3.set_orbs() = { k };
///  opD3.set_initialised() = true;
///  opD3.set_fermion() = true;
///  opD3.set_deltaQuantum() = getSpinQuantum(k);
///  opD3.allocate( b.get_stateInfo() );
///  // Expand small op to big_op
///  if (leftBlock->get_op_array(CRE).has(k)) {
///    const boost::shared_ptr<SparseMatrix> opC = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
///    SpinAdapted::operatorfunctions::TensorTrace(leftBlock, Transposeview(*opC), &b, &(b.get_stateInfo()), opD3);
///  }
///  else if (rightBlock->get_op_array(CRE).has(k)) {
///    const boost::shared_ptr<SparseMatrix> opC = rightBlock->get_op_rep(CRE, getSpinQuantum(k), k);
///    SpinAdapted::operatorfunctions::TensorTrace(rightBlock, Transposeview(*opC), &b, &(b.get_stateInfo()), opD3);
///  }
///  else {
///    pout << "leftBlock->get_op_array(CRE).has(k) = "  << leftBlock->get_op_array(CRE).has(k) << std::endl;
///    pout << "rightBlock->get_op_array(CRE).has(k) = "  << rightBlock->get_op_array(CRE).has(k) << std::endl;
///
///    pout << "leftleftBlock->get_op_array(CRE).has(k) = "  << leftleftBlock->get_op_array(CRE).has(k) << std::endl;
///    pout << "leftrightBlock->get_op_array(CRE).has(k) = "  << leftrightBlock->get_op_array(CRE).has(k) << std::endl;
///
///    pout << "rightleftBlock->get_op_array(CRE).has(k) = "  << rightleftBlock->get_op_array(CRE).has(k) << std::endl;
///    pout << "rightrightBlock->get_op_array(CRE).has(k) = "  << rightrightBlock->get_op_array(CRE).has(k) << std::endl;
///    assert(false);
///  }
///
///  // Build 2-index CC operator
///  CreCre opCC;
///  opCC.set_orbs() = { i, j };
///  opCC.set_initialised() = true;
///  opCC.set_fermion() = false;
///  opCC.set_deltaQuantum() = deltaQuantum12;
///  opCC.allocate( b.get_stateInfo() );
///  SpinAdapted::operatorfunctions::Product(&b, opC1, opC2, opCC, 1.0 );
///pout << "CC op\n";
///pout << opCC;
///
///  // Build 3-index operator as ((CreCre)Des)
///  SpinAdapted::operatorfunctions::Product(&b, opCC, opD3, *this, 1.0 );
///pout << "CCD op\n";
///pout << *this;
///
/////////    Transposeview op2 = Transposeview(leftBlock->get_op_rep(CRE, getSpinQuantum(k), k));
/////////    double parity = getCommuteParity(op1->get_deltaQuantum(), op2.get_deltaQuantum(), get_deltaQuantum());
/////////    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, op2, &b, &(b.get_stateInfo()), *this, 1.0*parity);
///
///
///  dmrginp.makeopsT -> stop();
///
///}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//MAW debug alternative to REDMATRIX routine
void SpinAdapted::CreCreDes::build_in_csf_space(const SpinBlock& b)
{
//pout << "building CreCreDes in CSF space as a product..\n";
assert(false);
////  built = true;
////  allocate(b.get_stateInfo());
////  Sign = 1;
////
//////FIXME Expensive? Cheaper to use pre-defined 2-index ops to build ((CRE_CRE)*DES) ?
////  const int i = get_orbs()[0];
////  const int j = get_orbs()[1];
////  const int k = get_orbs()[2];
//////pout << "indices = " << i << " " << j << " " << k << std::endl;
//////FIXME THIS FAILS SOMETIMES... e.g. c2_d2h test example; don't know why.
////  assert( b.get_op_array(CRE).has(i) );
////  assert( b.get_op_array(CRE).has(j) );
////  assert( b.get_op_array(CRE).has(k) );
////  // Spin quanta
////  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder();
////  assert( quantum_ladder.size() == 2 );
////  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
////  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
////  assert( deltaQuantum == deltaQuantum123 );
////
//////FIXME
////  if (i>=j) {
////    // Use predefined 2-index operator
////    assert( b.get_op_array(CRE_CRE).has(i,j) );
////    const boost::shared_ptr<SparseMatrix> opCC = b.get_op_rep(CRE_CRE, deltaQuantum12, i, j);
////    // Build 3-index as ((CreCre)Des)
////    Transposeview opD3 = Transposeview( b.get_op_rep(CRE, getSpinQuantum(k), k) );
////    operatorfunctions::Product(&b, *opCC, opD3, *this, 1.0 );
////  }
////  else {
////    // Build 2-index as Product
////    const boost::shared_ptr<SparseMatrix> opC1 = b.get_op_rep(CRE, getSpinQuantum(i), i);
////    const boost::shared_ptr<SparseMatrix> opC2 = b.get_op_rep(CRE, getSpinQuantum(j), j);
////    CreCre opCC;
////    opCC.set_orbs() = get_orbs();
////    opCC.set_orbs().pop_back();
////    opCC.set_initialised() = true;
////    opCC.set_fermion() = false;
////    opCC.set_deltaQuantum() = deltaQuantum12;
////    opCC.allocate(b.get_stateInfo());
////    operatorfunctions::Product(&b, *opC1, *opC2, opCC, 1.0 );
////    // Build 3-index as ((CreCre)Des)
////    Transposeview opD3 = Transposeview( b.get_op_rep(CRE, getSpinQuantum(k), k) );
////    operatorfunctions::Product(&b, opCC, opD3, *this, 1.0 );
////  }
////
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

double SpinAdapted::CreCreDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
//pout << "building CreCreDes explicitly from CSF..\n";
  assert( build_pattern == "((CC)D)" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];

  // Must take into account how the 3-index is built from a combination of 2-index and 1-index
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("((CC)D)");
  assert( quantum_ladder.size() == 2 );
  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  deltaQuantum = deltaQuantum123;

  // Spin quantum data for first pair of operators combined
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s();
  // Spin quantum data for total operator
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123 = deltaQuantum123.get_s();

//pout << "redMatrixElement indices:\n";
//pout << I << "  " << J << "  " << K << std::endl;
//pout << "spin composition:\n";
//pout << spin12/2.0 << "  " << spin123/2.0 << std::endl;
  TensorOp C1(I, 1); 
  TensorOp C2(J, 1); 
  TensorOp D3(K,-1); 

  // Combine first two operators
//FIXME  TensorOp CC = C1.product(C2, spin12, irrep12);  I==J has no affect
  TensorOp CC = C1.product(C2, spin12, irrep12, I==J);
  // Combine with third operator
  TensorOp CCD = CC.product(D3, spin123, irrep123);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum, ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CCD, ladder[i]) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CreCreDes::getworkingrepresentation(const SpinBlock* block)
{
  assert(this->get_initialised());
  if (this->get_built()) {
    return boost::shared_ptr<CreCreDes>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
  }
  else {
    boost::shared_ptr<SparseMatrix> rep(new CreCreDes);
    *rep = *this;
    rep->build(*block);
    return rep;
  }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Des,Des)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinAdapted::CreDesDes::build(const SpinBlock& b)
{
pout << "building CreDesDes renormalized operator...\n";
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_stateInfo());
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
pout << "indices  " << i << " " << j << " " << k << std::endl;

  SpinBlock* sysBlock = b.get_leftBlock();
  SpinBlock* dotBlock = b.get_rightBlock();
pout << "sysBlock:\n";
pout << *sysBlock;
pout << "dotBlock:\n";
pout << *dotBlock;

//  SpinQuantum deltaQuantum12 = get_quantum_ladder().at("0").at(0);
//  assert( get_quantum_ladder().size() == 2 );

  // Sys has i,j,k
  if (sysBlock->get_op_array(CRE_DES_DES).has(i,j,k)) {
pout << "maw sys(i,j,k)\n";
    std::string build_pattern_old = sysBlock->get_op_array(CRE_DES_DES).get_element(i,j,k).at(0)->get_build_pattern();
    assert( build_pattern_old == sysBlock->get_op_array(CRE_DES_DES).get_element(i,j,k).at(1)->get_build_pattern() );
    assert( build_pattern_old == sysBlock->get_op_array(CRE_DES_DES).get_element(i,j,k).at(2)->get_build_pattern() );
    const boost::shared_ptr<SparseMatrix>& op = sysBlock->get_op_rep(CRE_DES_DES, quantum_ladder.at(build_pattern_old), i,j,k);
    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *op, &b, &(b.get_stateInfo()), *this);
    build_pattern = build_pattern_old;
  }
  // Dot has i,j,k
  else if (dotBlock->get_op_array(CRE_DES_DES).has(i,j,k)) {
pout << "maw dot(i,j,k)\n";
    std::string build_pattern_old = dotBlock->get_op_array(CRE_DES_DES).get_element(i,j,k).at(0)->get_build_pattern();
    assert( build_pattern_old == dotBlock->get_op_array(CRE_DES_DES).get_element(i,j,k).at(1)->get_build_pattern() );
    assert( build_pattern_old == dotBlock->get_op_array(CRE_DES_DES).get_element(i,j,k).at(2)->get_build_pattern() );
    const boost::shared_ptr<SparseMatrix>& op = dotBlock->get_op_rep(CRE_DES_DES, quantum_ladder.at(build_pattern_old), i,j,k);
    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *op, &b, &(b.get_stateInfo()), *this);
    build_pattern = build_pattern_old;
  }
  else assert(false);

  dmrginp.makeopsT -> stop();
pout << "CDD op\n";
pout << *this;
pout << "done building CreDesDes renormalized operator!\n";

}

////
////
////  // Sys has i,j,k
////  if (sysBlock->get_op_array(CRE_DES_DES).has(i,j,k)) {
////pout << "maw sys(i,j,k)\n";
//////    const boost::shared_ptr<SparseMatrix>& op = sysBlock->get_op_rep(CRE_DES_DES, quantum_ladder, i,j,k);
//////    build_pattern = op->get_build_pattern();
//////pout << "opCDD\n";
//////pout << *op;
//////    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *op, &b, &(b.get_stateInfo()), *this);
//////pout << "opCDD\n";
//////pout << *this;
////  }
////  // Sys has i,j only
////  else if (sysBlock->get_op_array(CRE_DES).has(i,j)) {
////pout << "maw sys(i,j)\n";
////    assert( dotBlock->get_op_array(CRE).has(k) );
////    const boost::shared_ptr<SparseMatrix>& opCD = sysBlock->get_op_rep(CRE_DES, deltaQuantum12, i,j);
////    const boost::shared_ptr<SparseMatrix>& opC3 = dotBlock->get_op_rep(CRE, getSpinQuantum(k), k);
////    // Build ((CD)D)
////    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, *opCD, Transposeview(*opC3), &b, &(b.get_stateInfo()), *this, 1.0);
////  }
////  // Sys has j,k only
////  else if (sysBlock->get_op_array(CRE_CRE).has(j,k)) {
////pout << "maw sys(j,k)\n";
////    assert( dotBlock->get_op_array(CRE).has(i) );
////    assert( sysBlock->get_op_array(CRE).has(j) );
////    assert( sysBlock->get_op_array(CRE).has(k) );
////    // Build (CiDj) first
////    CreCre opCD;
////    opCD.set_orbs() = { i, j };
////    opCD.set_initialised() = true;
////    opCD.set_fermion() = false;
////    opCD.set_deltaQuantum() = deltaQuantum12;
////    opCD.allocate( b.get_stateInfo() );
////    const boost::shared_ptr<SparseMatrix>& opC1 = dotBlock->get_op_rep(CRE, getSpinQuantum(i), i);
////    Transposeview opD2 = Transposeview( sysBlock->get_op_rep(CRE, getSpinQuantum(j), j) );
////    double parity = getCommuteParity( opC1->get_deltaQuantum(), opD2.get_deltaQuantum(), opCD.get_deltaQuantum() );
////    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *opC1, opD2, &b, &(b.get_stateInfo()), opCD, 1.0*parity);
////    // Build C3 on full spinblock
////    const boost::shared_ptr<SparseMatrix>& opC3 = sysBlock->get_op_rep(CRE, getSpinQuantum(k), k);
////    CreCre opC3big;
////    opC3big.set_orbs() = { k };
////    opC3big.set_initialised() = true;
////    opC3big.set_fermion() = true;
////    opC3big.set_deltaQuantum() = getSpinQuantum(k);
////    opC3big.allocate( b.get_stateInfo() );
////    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *opC3, &b, &(b.get_stateInfo()), opC3big);
////    // Build ((CD)D) on full spinblock
////    SpinAdapted::operatorfunctions::Product(&b, opCD, Transposeview(opC3big), *this, 1.0 );
////  }
////  // Sys has i only
////  else if (sysBlock->get_op_array(CRE).has(i)) {
////pout << "maw sys(i)\n";
////    assert( j == k );
////    assert( dotBlock->get_op_array(CRE_DES).has(j,k) );
////    assert( dotBlock->get_op_array(CRE).has(j) );
////    // Build (CiDj) first
////    CreCre opCD;
////    opCD.set_orbs() = { i, j };
////    opCD.set_initialised() = true;
////    opCD.set_fermion() = false;
////    opCD.set_deltaQuantum() = deltaQuantum12;
////    opCD.allocate( b.get_stateInfo() );
////    const boost::shared_ptr<SparseMatrix>& opC1 = sysBlock->get_op_rep(CRE, getSpinQuantum(i), i);
////    Transposeview opD2 = Transposeview( dotBlock->get_op_rep(CRE, getSpinQuantum(j), j) );
////    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, *opC1, opD2, &b, &(b.get_stateInfo()), opCD, 1.0);
////    // Build C3 on full spinblock
////    CreCre opC3big;
////    opC3big.set_orbs() = { k };
////    opC3big.set_initialised() = true;
////    opC3big.set_fermion() = true;
////    opC3big.set_deltaQuantum() = getSpinQuantum(k);
////    opC3big.allocate( b.get_stateInfo() );
////    const boost::shared_ptr<SparseMatrix>& opC3 = dotBlock->get_op_rep(CRE, getSpinQuantum(k), k);
////    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *opC3, &b, &(b.get_stateInfo()), opC3big);
////    // Build ((CD)D) on full spinblock
////    SpinAdapted::operatorfunctions::Product(&b, opCD, Transposeview(opC3big), *this, 1.0 );
////  }
////  // Sys has k only
////  else if (sysBlock->get_op_array(CRE).has(k)) {
////pout << "maw sys(k)\n";
////    assert( i == j );
////    assert( dotBlock->get_op_array(CRE_DES).has(i,j) );
////    const boost::shared_ptr<SparseMatrix>& opCD = dotBlock->get_op_rep(CRE_DES, deltaQuantum12, i,j);
////    Transposeview opD = Transposeview( sysBlock->get_op_rep(CRE, getSpinQuantum(k), k) );
////    // Build ((CD)D)
////    double parity = getCommuteParity( opCD->get_deltaQuantum(), opD.get_deltaQuantum(), get_deltaQuantum() );
////    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *opCD, opD, &b, &(b.get_stateInfo()), *this, 1.0*parity);
////  }
////  // Dot has i,j,k
////  else if (dotBlock->get_op_array(CRE_DES_DES).has(i,j,k)) {
////pout << "maw dot(i,j,k)\n";
//////    const boost::shared_ptr<SparseMatrix>& op = dotBlock->get_op_rep(CRE_DES_DES, quantum_ladder, i,j,k);
//////    build_pattern = op->get_build_pattern();
//////    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *op, &b, &(b.get_stateInfo()), *this);
////  }
////  else assert(false);
////
//////pout << "CDD op\n";
//////pout << *this;
////  dmrginp.makeopsT -> stop();
////pout << "done building CDD renormalized operator!\n";
////
////}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

double SpinAdapted::CreDesDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
//pout << "building CreDesDes explicitly from CSF..\n";
  assert( build_pattern == "((CD)D)" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];

  // Must take into account how the 3-index is built from a combination of 2-index and 1-index
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("((CD)D)");
  assert( quantum_ladder.size() == 2 );
  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  deltaQuantum = deltaQuantum123;

  // Spin quantum data for first pair of operators combined
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s();
  // Spin quantum data for total operator
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123 = deltaQuantum123.get_s();

//pout << "redMatrixElement indices:\n";
//pout << I << "  " << J << "  " << K << std::endl;
//pout << "spin composition:\n";
//pout << spin12/2.0 << "  " << spin123/2.0 << std::endl;
  TensorOp C1(I, 1); 
  TensorOp D2(J,-1); 
  TensorOp D3(K,-1); 

  // Combine first two operators
  TensorOp CD = C1.product(D2, spin12, irrep12);
  // Combine with third operator
  TensorOp CDD = CD.product(D3, spin123, irrep123);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum, ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CDD, ladder[i]) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CreDesDes::getworkingrepresentation(const SpinBlock* block)
{
  assert(this->get_initialised());
  if (this->get_built()) {
//pout << "maw getworkingrep: CDD get pre-built   " << get_build_pattern() << std::endl;
    return boost::shared_ptr<CreDesDes>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
  }
  else {
//pout << "maw getworkingrep: CDD build\n";
    boost::shared_ptr<SparseMatrix> rep(new CreDesDes);
//pout << "maw getworkingrep: build_pattern before =  " << rep->get_build_pattern() << std::endl;
    *rep = *this;
    rep->build(*block);
//pout << "maw getworkingrep: build_pattern after =  " << rep->get_build_pattern() << std::endl;

    return rep;
  }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Des,Cre)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinAdapted::CreDesCre::build(const SpinBlock& b)
{
assert(false);
//////pout << "building CreDesCre renormalized operator...\n";
////  dmrginp.makeopsT -> start();
////  built = true;
////  allocate(b.get_stateInfo());
////  Sign = 1;
////
////  const int i = get_orbs()[0];
////  const int j = get_orbs()[1];
////  const int k = get_orbs()[2];
//////pout << "indices  " << i << " " << j << " " << k << std::endl;
////
////  SpinBlock* leftBlock = b.get_leftBlock();
////  SpinBlock* rightBlock = b.get_rightBlock();
////
////  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder();
////  assert( quantum_ladder.size() == 2 );
////  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
////  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
////
////  // Below, note transposes for (i,j) since 2-index ops are built with i>=j 
////
////  // Forward
////  if (leftBlock->get_op_array(CRE_DES_CRE).has(i,j,k)) {
////    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE_DES_CRE, quantum_ladder, i,j,k);
////    SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
////  }
////  else if (rightBlock->get_op_array(CRE_DES_CRE).has(i,j,k)) {
////    const boost::shared_ptr<SparseMatrix>& op = rightBlock->get_op_rep(CRE_DES_CRE, quantum_ladder, i,j,k);
////    SpinAdapted::operatorfunctions::TensorTrace(rightBlock, *op, &b, &(b.get_stateInfo()), *this);
////  }
//////  else if (leftBlock->get_op_array(CRE_CRE).has(j,i)) {
//////    assert (rightBlock->get_op_array(CRE).has(k));
//////    const boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(CRE_CRE, deltaQuantum12, j,i);
//////    Transposeview op2 = Transposeview(rightBlock->get_op_rep(CRE, getSpinQuantum(k), k));
//////    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, op2, &b, &(b.get_stateInfo()), *this, 1.0);
//////  }
//////  // Want to build as ((CC)D), not (C(CD)), but FIXME this is relatively expensive!
//////  else if (b.get_op_array(CRE_CRE).has(j,i)) {
//////    // Get CC on big block
//////    const boost::shared_ptr<SparseMatrix> op1 = b.get_op_rep(CRE_CRE, deltaQuantum12, j,i);
//////    // Get D on big block
//////    assert (b.get_op_array(CRE).has(k));
//////    Transposeview op2 = Transposeview( b.get_op_rep(CRE, getSpinQuantum(k), k) );
//////    // Build tensor product operator directly in 4M*4M space
//////    boost::shared_ptr<SparseMatrix> newOp (new Cre);
//////    newOp->set_orbs().push_back(i);
//////    newOp->set_orbs().push_back(j);
//////    newOp->set_orbs().push_back(k);
//////    newOp->set_initialised() = true;
//////    newOp->set_fermion() = true; //is_fermion;
//////    operatorfunctions::Product(&b, *op1, op2, *newOp, 1.0 );
//////  }
////////  else if (rightBlock->get_op_array(CRE_CRE).has(j,i)) {
////////pout << "hello6\n";
////////    assert (false);
////////    assert (leftBlock->get_op_array(CRE).has(k));
////////    const boost::shared_ptr<SparseMatrix> op1 = rightBlock->get_op_rep(CRE_CRE, deltaQuantum12, i,j);
////////    Transposeview op2 = Transposeview(leftBlock->get_op_rep(CRE, getSpinQuantum(k), k));
////////    double parity = getCommuteParity(op1->get_deltaQuantum(), op2.get_deltaQuantum(), get_deltaQuantum());
////////    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, op2, &b, &(b.get_stateInfo()), *this, 1.0*parity);
////////  }
////////  else if (rightBlock->get_op_array(CRE).has(i)) {
////////pout << "hello7\n";
////////    assert (false);
////////  }
//////  else {
//////    assert (false);
//////  }
////
////  dmrginp.makeopsT -> stop();
////
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//MAW debug alternative to REDMATRIX routine
void SpinAdapted::CreDesCre::build_in_csf_space(const SpinBlock& b)
{
assert(false);
//////pout << "building CreDesCre in CSF space as a product..\n";
////  built = true;
////  allocate(b.get_stateInfo());
////  Sign = 1;
////
////  const int i = get_orbs()[0];
////  const int j = get_orbs()[1];
////  const int k = get_orbs()[2];
//////pout << "indices = " << i << " " << j << " " << k << std::endl;
////  assert( b.get_op_array(CRE).has(i) );
////  assert( b.get_op_array(CRE).has(j) );
////  assert( b.get_op_array(CRE).has(k) );
////  // Spin quanta
////  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder();
////  assert( quantum_ladder.size() == 2 );
////  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
////  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
////  assert( deltaQuantum == deltaQuantum123 );
////
////  // Build 2-index as Product first
//////  const boost::shared_ptr<SparseMatrix> opC = b.get_op_rep(CRE, getSpinQuantum(i), i);
//////  Transposeview opD = Transposeview( b.get_op_rep(CRE, getSpinQuantum(j), j) );
//////  CreDes opCD;
//////  opCD.set_orbs() = get_orbs();
//////  opCD.set_orbs().pop_back();
//////  opCD.set_initialised() = true;
//////  // 2-index operators are bosonic
//////  opCD.set_fermion() = false;
//////  opCD.set_deltaQuantum() = deltaQuantum12;
////////NOTE allocate must not come too early!!!
//////  opCD.allocate(b.get_stateInfo());
//////  operatorfunctions::Product(&b, *opC, opD, opCD, 1.0 );
////
//////  // Can do as above, or invoke pre-defined 2-index operator instead
////  const boost::shared_ptr<SparseMatrix> opCD = b.get_op_rep(CRE_DES, deltaQuantum12, i, j);
////
////  // Build 3-index as ((CreDes)Cre)
////  const boost::shared_ptr<SparseMatrix> opC = b.get_op_rep(CRE, getSpinQuantum(k), k);
////  operatorfunctions::Product(&b, *opCD, *opC, *this, 1.0 );
}


//-------------------------------------------------------------------------------------------------------------------------------------------------------------

double SpinAdapted::CreDesCre::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
assert(false);
//////pout << "building CreDesCre in CSF space explicitly..\n";
////  double element = 0.0;
////  int I = get_orbs()[0]; 
////  int J = get_orbs()[1];
////  int K = get_orbs()[2];
////
////  // Must take into account how the 3-index is built from a combination of 2-index and 1-index
////  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder();
////  assert( quantum_ladder.size() == 2 );
////  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
////  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
////  assert( deltaQuantum == deltaQuantum123 );
////
////  // Spin quantum data for first pair of operators combined
////  IrrepSpace sym12 = deltaQuantum12.get_symm();
////  int irrep12 = deltaQuantum12.get_symm().getirrep();
////  int spin12 = deltaQuantum12.get_s();
////  // Spin quantum data for total operator
////  IrrepSpace sym123 = deltaQuantum123.get_symm();
////  int irrep123 = deltaQuantum123.get_symm().getirrep();
////  int spin123 = deltaQuantum123.get_s();
////
//////pout << "redMatrixElement indices:\n";
//////pout << I << "  " << J << "  " << K << std::endl;
//////pout << "spin composition:\n";
//////pout << spin12/2.0 << "  " << spin123/2.0 << std::endl;
////  TensorOp C1(I, 1); 
////  TensorOp D2(J,-1); 
////  TensorOp C3(K, 1); 
////
////  // Combine first two operators
////  TensorOp CD = C1.product(D2, spin12, irrep12);
////  // Combine with third operator
////  TensorOp CDC = CD.product(C3, spin123, irrep123);
////
////  for (int i=0; i<ladder.size(); i++)
////  {
////    int index = 0; double cleb=0.0;
////    if (nonZeroTensorComponent(c1, deltaQuantum, ladder[i], index, cleb)) {
////      std::vector<double> MatElements = calcMatrixElements(c1, CDC, ladder[i]) ;
////      element = MatElements[index]/cleb;
////      break;
////    }
////    else
////      continue;
////  }
////  return element;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CreDesCre::getworkingrepresentation(const SpinBlock* block)
{
assert(false);
////  assert(this->get_initialised());
////  if (this->get_built()) {
////    return boost::shared_ptr<CreDesCre>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
////  }
////  else {
////    boost::shared_ptr<SparseMatrix> rep(new CreDesCre);
////    *rep = *this;
////    rep->build(*block);
////    return rep;
////  }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------


