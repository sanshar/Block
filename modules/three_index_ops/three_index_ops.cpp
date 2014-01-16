/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//
// Extension of Operators.C for 3-index operators
//FIXME there's a lot of duplication, especially in build_from_disk... Templates??
//
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "op_components.h"
#include "BaseOperator.h"
#include "spinblock.h"
#include "operatorfunctions.h"
#include "tensor_operator.h"
#include "three_index_ops.h"

//===========================================================================================================================================================
// 3PDM operators
//===========================================================================================================================================================

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Cre,Des)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinAdapted::CreCreDes::build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs)
{
  // Note that we should only allocate after deltaQuantum is finally determined
//cout << "building CreCreDes renormalized operator on disk...\n";
//FIXME timer
//  dmrginp.makeopsT -> start();
  built = true;
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
//cout << "i,j,k = " << i << " " << j << " " << k << endl;

  SpinBlock* sysBlock = b.get_leftBlock();
  SpinBlock* dotBlock = b.get_rightBlock();

  // Sys has i,j,k
  if (sysBlock->get_op_array(CRE_CRE_DES).has_local_index(i,j,k)) {
//cout << "maw disk sys(i,j,k)\n";
    assert( sysBlock->get_op_array(CRE_CRE_DES).get_element(i,j,k).at(0)->get_built_on_disk() );
    // Retrieve from disk 3-index operator on sys block
    boost::shared_ptr<SparseMatrix> op (new CreCreDes);
    boost::archive::binary_iarchive load_op(sysfs);
    load_op >> *op;
    // Assume this is the operator we want, and expand it to the big block
    assert( i == op->get_orbs()[0] );
    assert( j == op->get_orbs()[1] );
    assert( k == op->get_orbs()[2] );
    // Build according to previous build_pattern
    std::string build_pattern_old = op->get_build_pattern();
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  // Dot has i,j,k
  else if (dotBlock->get_op_array(CRE_CRE_DES).has_local_index(i,j,k)) {
//cout << "maw dot(i,j,k)\n";
    assert( dotBlock->get_op_array(CRE_CRE_DES).get_element(i,j,k).at(0)->get_built_on_disk() );
    // Retrieve from disk 3-index operator on dot block
    boost::shared_ptr<SparseMatrix> op (new CreCreDes);
    boost::archive::binary_iarchive load_op(dotfs);
    load_op >> *op;
    // Assume this is the operator we want, and expand it to the big block
    assert( i == op->get_orbs()[0] );
    assert( j == op->get_orbs()[1] );
    assert( k == op->get_orbs()[2] );
    // Build according to previous build_pattern
    std::string build_pattern_old = op->get_build_pattern();
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  else {
    // Build from in-core 1 and 2-index operators
    build(b);
  }

//  dmrginp.makeopsT -> stop();

}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinAdapted::CreCreDes::build(const SpinBlock& b)
{
//cout << "building CreCreDes renormalized operator...\n";
//build_in_csf_space(b); //fails

  dmrginp.makeopsT -> start();
  built = true;
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
//pout << "indices  " << i << " " << j << " " << k << std::endl;

  SpinBlock* sysBlock = b.get_leftBlock();
  SpinBlock* dotBlock = b.get_rightBlock();
//pout << "sysBlock:\n";
//pout << *sysBlock;
//pout << "dotBlock:\n";
//pout << *dotBlock;

//  SpinQuantum deltaQuantum12 = get_quantum_ladder().at("0").at(0);
//  assert( get_quantum_ladder().size() == 2 );

  // Sys has i,j,k
  if (sysBlock->get_op_array(CRE_CRE_DES).has_local_index(i,j,k)) {
//pout << "maw sys(i,j,k)\n";
    std::string build_pattern_old = sysBlock->get_op_array(CRE_CRE_DES).get_element(i,j,k).at(0)->get_build_pattern();
    assert( build_pattern_old == sysBlock->get_op_array(CRE_CRE_DES).get_element(i,j,k).at(1)->get_build_pattern() );
    assert( build_pattern_old == sysBlock->get_op_array(CRE_CRE_DES).get_element(i,j,k).at(2)->get_build_pattern() );
    const boost::shared_ptr<SparseMatrix>& op = sysBlock->get_op_rep(CRE_CRE_DES, quantum_ladder.at(build_pattern_old), i,j,k);
    // Build according to previous build_pattern
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  // Dot has i,j,k
  else if (dotBlock->get_op_array(CRE_CRE_DES).has_local_index(i,j,k)) {
//pout << "maw dot(i,j,k)\n";
    assert( i == j );
    assert( j == k );
    std::string build_pattern_old = dotBlock->get_op_array(CRE_CRE_DES).get_element(i,j,k).at(0)->get_build_pattern();
    assert( build_pattern_old == dotBlock->get_op_array(CRE_CRE_DES).get_element(i,j,k).at(1)->get_build_pattern() );
    assert( build_pattern_old == dotBlock->get_op_array(CRE_CRE_DES).get_element(i,j,k).at(2)->get_build_pattern() );
    const boost::shared_ptr<SparseMatrix>& op = dotBlock->get_op_rep(CRE_CRE_DES, quantum_ladder.at(build_pattern_old), i,j,k);
    // Build according to previous build_pattern
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  // Dot has j,k;
  else if (dotBlock->get_op_array(CRE_DES).has_local_index(j,k)) {
//pout << "maw dot(j,k)\n";
    assert( j == k );
    assert( sysBlock->get_op_array(CRE).has_local_index(i) );
    build_pattern = "(C(CD))";
    const boost::shared_ptr<SparseMatrix>& opCD = dotBlock->get_op_rep(CRE_DES, quantum_ladder.at(build_pattern).at(0), j,k);
    const boost::shared_ptr<SparseMatrix>& opC = sysBlock->get_op_rep(CRE, getSpinQuantum(i), i);
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, *opC, *opCD, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  // Sys has i,j
  else if (sysBlock->get_op_array(CRE_CRE).has_local_index(i,j)) {
//pout << "maw sys(i,j)\n";
    assert( dotBlock->get_op_array(CRE).has_local_index(k) );
    build_pattern = "((CC)D)";
    const boost::shared_ptr<SparseMatrix>& opCC = sysBlock->get_op_rep(CRE_CRE, quantum_ladder.at(build_pattern).at(0), i,j);
    Transposeview opD = Transposeview( dotBlock->get_op_rep(CRE, getSpinQuantum(k), k) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, *opCC, opD, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  // Sys has j,k;
  else if (sysBlock->get_op_array(CRE_DES).has_local_index(j,k)) {
//pout << "maw sys(j,k)\n";
    assert( dotBlock->get_op_array(CRE).has_local_index(i) );
    build_pattern = "(C(CD))";
    const boost::shared_ptr<SparseMatrix>& opC = dotBlock->get_op_rep(CRE, getSpinQuantum(i), i);
    const boost::shared_ptr<SparseMatrix>& opCD = sysBlock->get_op_rep(CRE_DES, quantum_ladder.at(build_pattern).at(0), j,k);
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    // Tensor product of dot*sys so need to take into account parity factors
    double parity = getCommuteParity( opC->get_deltaQuantum(), opCD->get_deltaQuantum(), get_deltaQuantum() );
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *opC, *opCD, &b, &(b.get_stateInfo()), *this, 1.0*parity);
  }
  // Dot has i,j
  else if (dotBlock->get_op_array(CRE_CRE).has_local_index(i,j)) {
//pout << "maw dot(i,j)\n";
    assert( i == j );
    assert( sysBlock->get_op_array(CRE).has_local_index(k) );
    build_pattern = "((CC)D)";
    const boost::shared_ptr<SparseMatrix>& opCC = dotBlock->get_op_rep(CRE_CRE, quantum_ladder.at(build_pattern).at(0), i,j);
    Transposeview opD = Transposeview( sysBlock->get_op_rep(CRE, getSpinQuantum(k), k) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    // Tensor product of dot*sys so need to take into account parity factors
    double parity = getCommuteParity( opCC->get_deltaQuantum(), opD.get_deltaQuantum(), get_deltaQuantum() );
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *opCC, opD, &b, &(b.get_stateInfo()), *this, 1.0*parity);
  }
  else assert(false);

  dmrginp.makeopsT -> stop();
//pout << "CCD op\n";
//pout << *this;
//pout << "done building CreCreDes renormalized operator!\n";

}

////
////
////  // Sys has i,j only
////  else if (sysBlock->get_op_array(CRE_CRE).has_local_index(i,j)) {
//////pout << "maw sys(i,j)\n";
////    assert( dotBlock->get_op_array(CRE).has_local_index(k) );
////    const boost::shared_ptr<SparseMatrix>& opCC = sysBlock->get_op_rep(CRE_CRE, deltaQuantum12, i,j);
////    const boost::shared_ptr<SparseMatrix>& opC3 = dotBlock->get_op_rep(CRE, getSpinQuantum(k), k);
////    // Build ((CC)D)
////    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, *opCC, Transposeview(*opC3), &b, &(b.get_stateInfo()), *this, 1.0);
////  }
/////////////////////////////////////////////////
////  // Sys has j,k only
////  else if (sysBlock->get_op_array(CRE_DES).has_local_index(j,k)) {
//////pout << "maw sys(j,k)\n";
////////    assert( dotBlock->get_op_array(CRE).has_local_index(i) );
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
////////    assert( dotBlock->get_op_array(CRE).has_local_index(i) );
////////    assert( sysBlock->get_op_array(CRE).has_local_index(j) );
////////    assert( sysBlock->get_op_array(CRE).has_local_index(k) );
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
////  else if (sysBlock->get_op_array(CRE).has_local_index(i)) {
//////pout << "maw sys(i)\n";
////    assert( j == k );
////    assert( dotBlock->get_op_array(CRE_DES).has_local_index(j,k) );
////    assert( dotBlock->get_op_array(CRE).has_local_index(j) );
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
////  else if (sysBlock->get_op_array(CRE).has_local_index(k)) {
//////pout << "maw sys(k)\n";
////    assert( i == j );
////    assert( dotBlock->get_op_array(CRE_DES).has_local_index(i,j) );
////    const boost::shared_ptr<SparseMatrix>& opCC = dotBlock->get_op_rep(CRE_CRE, deltaQuantum12, i,j);
////    Transposeview opD = Transposeview( sysBlock->get_op_rep(CRE, getSpinQuantum(k), k) );
////    // Build ((CC)D)
////    double parity = getCommuteParity( opCC->get_deltaQuantum(), opD.get_deltaQuantum(), get_deltaQuantum() );
////    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *opCC, opD, &b, &(b.get_stateInfo()), *this, 1.0*parity);
////  }
////  // Dot has i,j,k
////  else if (dotBlock->get_op_array(CRE_CRE_DES).has_local_index(i,j,k)) {
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
////  if (sysBlock->get_op_array(CRE_CRE_DES).has_local_index(i,j,k)) {
////pout << "maw forward(xxx,0)\n";
////    const boost::shared_ptr<SparseMatrix>& op = sysBlock->get_op_rep(CRE_CRE_DES, quantum_ladder, i,j,k);
////pout << "opCCD\n";
////pout << *op;
////    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *op, &b, &(b.get_stateInfo()), *this);
////  }
////  // (XX,X)
////  else if (sysBlock->get_op_array(CRE_CRE).has_local_index(i,j)) {
////pout << "maw forward(XX,X)\n";
////    assert( sysBlock->get_op_array(CRE).has_local_index(i) );
////    assert( sysBlock->get_op_array(CRE).has_local_index(j) );
////    assert( dotBlock->get_op_array(CRE).has_local_index(k) );
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
////  else if (sysBlock->get_op_array(CRE).has_local_index(i)) {
////pout << "maw forward(X,XX)\n";
////    assert( j == k );
////    assert( sysBlock->get_op_array(CRE).has_local_index(i) );
////    assert( dotBlock->get_op_array(CRE).has_local_index(j) );
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
////  else if (dotBlock->get_op_array(CRE_CRE_DES).has_local_index(i,j,k)) {
////pout << "maw backward(0,xxx)\n";
////    const boost::shared_ptr<SparseMatrix>& op = dotBlock->get_op_rep(CRE_CRE_DES, quantum_ladder, i,j,k);
////    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *op, &b, &(b.get_stateInfo()), *this);
////  }
////  // (X,XX)
////  else if (sysBlock->get_op_array(CRE_CRE).has_local_index(j,k)) {
////pout << "maw backward(X,XX)\n";
////    assert( dotBlock->get_op_array(CRE).has_local_index(i) );
////    assert( sysBlock->get_op_array(CRE).has_local_index(j) );
////    assert( sysBlock->get_op_array(CRE).has_local_index(k) );
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
////  else if (sysBlock->get_op_array(CRE).has_local_index(k)) {
////pout << "maw backward(XX,X)\n";
////    assert( i == j );
////    assert( dotBlock->get_op_array(CRE).has_local_index(i) );
////    assert( sysBlock->get_op_array(CRE).has_local_index(k) );
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
/////  if (leftBlock->get_op_array(CRE_CRE).has_local_index(i,j)) {
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
///  if (leftBlock->get_op_array(CRE).has_local_index(i)) {
///    const boost::shared_ptr<SparseMatrix> opC = leftBlock->get_op_rep(CRE, getSpinQuantum(i), i);
///    SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *opC, &b, &(b.get_stateInfo()), opC1);
///  }
///  else if (rightBlock->get_op_array(CRE).has_local_index(i)) {
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
///  if (leftBlock->get_op_array(CRE).has_local_index(j)) {
///    const boost::shared_ptr<SparseMatrix> opC = leftBlock->get_op_rep(CRE, getSpinQuantum(j), j);
///    SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *opC, &b, &(b.get_stateInfo()), opC2);
///  }
///  else if (rightBlock->get_op_array(CRE).has_local_index(j)) {
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
///  if (leftBlock->get_op_array(CRE).has_local_index(k)) {
///    const boost::shared_ptr<SparseMatrix> opC = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
///    SpinAdapted::operatorfunctions::TensorTrace(leftBlock, Transposeview(*opC), &b, &(b.get_stateInfo()), opD3);
///  }
///  else if (rightBlock->get_op_array(CRE).has_local_index(k)) {
///    const boost::shared_ptr<SparseMatrix> opC = rightBlock->get_op_rep(CRE, getSpinQuantum(k), k);
///    SpinAdapted::operatorfunctions::TensorTrace(rightBlock, Transposeview(*opC), &b, &(b.get_stateInfo()), opD3);
///  }
///  else {
///    pout << "leftBlock->get_op_array(CRE).has_local_index(k) = "  << leftBlock->get_op_array(CRE).has_local_index(k) << std::endl;
///    pout << "rightBlock->get_op_array(CRE).has_local_index(k) = "  << rightBlock->get_op_array(CRE).has_local_index(k) << std::endl;
///
///    pout << "leftleftBlock->get_op_array(CRE).has_local_index(k) = "  << leftleftBlock->get_op_array(CRE).has_local_index(k) << std::endl;
///    pout << "leftrightBlock->get_op_array(CRE).has_local_index(k) = "  << leftrightBlock->get_op_array(CRE).has_local_index(k) << std::endl;
///
///    pout << "rightleftBlock->get_op_array(CRE).has_local_index(k) = "  << rightleftBlock->get_op_array(CRE).has_local_index(k) << std::endl;
///    pout << "rightrightBlock->get_op_array(CRE).has_local_index(k) = "  << rightrightBlock->get_op_array(CRE).has_local_index(k) << std::endl;
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
//cout << "building CreCreDes in CSF space as a product..\n";
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
////  assert( b.get_op_array(CRE).has_local_index(i) );
////  assert( b.get_op_array(CRE).has_local_index(j) );
////  assert( b.get_op_array(CRE).has_local_index(k) );
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
////    assert( b.get_op_array(CRE_CRE).has_local_index(i,j) );
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
//cout << "building CreCreDes explicitly from CSF..\n";
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

void SpinAdapted::CreDesDes::build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs)
{
//cout << "building CreDesDes renormalized operator on disk...\n";
//FIXME timer
//  dmrginp.makeopsT -> start();
  built = true;
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
//cout << "i,j,k = " << i << " " << j << " " << k << endl;

  SpinBlock* sysBlock = b.get_leftBlock();
  SpinBlock* dotBlock = b.get_rightBlock();

  // Sys has i,j,k
  if (sysBlock->get_op_array(CRE_DES_DES).has_local_index(i,j,k)) {
//cout << "maw disk sys(i,j,k)\n";
    assert( sysBlock->get_op_array(CRE_DES_DES).get_element(i,j,k).at(0)->get_built_on_disk() );
    // Retrieve from disk 3-index operator on sys block
    boost::shared_ptr<SparseMatrix> op (new CreDesDes);
    boost::archive::binary_iarchive load_op(sysfs);
    load_op >> *op;
    // Assume this is the operator we want, and expand it to the big block
    assert( i == op->get_orbs()[0] );
    assert( j == op->get_orbs()[1] );
    assert( k == op->get_orbs()[2] );
    // Build according to previous build_pattern
    std::string build_pattern_old = op->get_build_pattern();
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  // Dot has i,j,k
  else if (dotBlock->get_op_array(CRE_DES_DES).has_local_index(i,j,k)) {
//cout << "maw dot(i,j,k)\n";
    assert( dotBlock->get_op_array(CRE_DES_DES).get_element(i,j,k).at(0)->get_built_on_disk() );
    // Retrieve from disk 3-index operator on dot block
    boost::shared_ptr<SparseMatrix> op (new CreDesDes);
    boost::archive::binary_iarchive load_op(dotfs);
    load_op >> *op;
    // Assume this is the operator we want, and expand it to the big block
    assert( i == op->get_orbs()[0] );
    assert( j == op->get_orbs()[1] );
    assert( k == op->get_orbs()[2] );
    // Build according to previous build_pattern
    std::string build_pattern_old = op->get_build_pattern();
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  else {
    // Build from in-core 1 and 2-index operators
    build(b);
  }

//  dmrginp.makeopsT -> stop();

}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinAdapted::CreDesDes::build(const SpinBlock& b)
{
//cout << "building CreDesDes renormalized operator...\n";
  dmrginp.makeopsT -> start();
  built = true;
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
//pout << "indices  " << i << " " << j << " " << k << std::endl;

  SpinBlock* sysBlock = b.get_leftBlock();
  SpinBlock* dotBlock = b.get_rightBlock();
//pout << "sysBlock:\n";
//pout << *sysBlock;
//pout << "dotBlock:\n";
//pout << *dotBlock;

  // Sys has i,j,k
  if (sysBlock->get_op_array(CRE_DES_DES).has_local_index(i,j,k)) {
//pout << "maw sys(i,j,k)\n";
    std::string build_pattern_old = sysBlock->get_op_array(CRE_DES_DES).get_element(i,j,k).at(0)->get_build_pattern();
    assert( build_pattern_old == sysBlock->get_op_array(CRE_DES_DES).get_element(i,j,k).at(1)->get_build_pattern() );
    assert( build_pattern_old == sysBlock->get_op_array(CRE_DES_DES).get_element(i,j,k).at(2)->get_build_pattern() );
    const boost::shared_ptr<SparseMatrix>& op = sysBlock->get_op_rep(CRE_DES_DES, quantum_ladder.at(build_pattern_old), i,j,k);
    // Build according to previous build_pattern
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  // Dot has i,j,k
  else if (dotBlock->get_op_array(CRE_DES_DES).has_local_index(i,j,k)) {
//pout << "maw dot(i,j,k)\n";
    assert( i == j );
    assert( j == k );
    std::string build_pattern_old = dotBlock->get_op_array(CRE_DES_DES).get_element(i,j,k).at(0)->get_build_pattern();
    assert( build_pattern_old == dotBlock->get_op_array(CRE_DES_DES).get_element(i,j,k).at(1)->get_build_pattern() );
    assert( build_pattern_old == dotBlock->get_op_array(CRE_DES_DES).get_element(i,j,k).at(2)->get_build_pattern() );
    const boost::shared_ptr<SparseMatrix>& op = dotBlock->get_op_rep(CRE_DES_DES, quantum_ladder.at(build_pattern_old), i,j,k);
    // Build according to previous build_pattern
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  // Dot has j,k;
  else if (dotBlock->get_op_array(CRE_CRE).has_local_index(j,k)) {
//pout << "maw dot(j,k)\n";
    assert( j == k );
    assert( sysBlock->get_op_array(CRE).has_local_index(i) );
    build_pattern = "(C(DD))";
    const boost::shared_ptr<SparseMatrix>& opC = sysBlock->get_op_rep(CRE, getSpinQuantum(i), i);
    Transposeview opDD = Transposeview( dotBlock->get_op_rep(CRE_CRE, quantum_ladder.at(build_pattern).at(0), j,k) );
    // Indices OK after transpose since j=k
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, *opC, opDD, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  // Sys has i,j
  else if (sysBlock->get_op_array(CRE_DES).has_local_index(i,j)) {
//pout << "maw dot(k)\n";
    assert( dotBlock->get_op_array(CRE).has_local_index(k) );
    build_pattern = "((CD)D)";
    const boost::shared_ptr<SparseMatrix>& opCD = sysBlock->get_op_rep(CRE_DES, quantum_ladder.at(build_pattern).at(0), i,j);
    Transposeview opD = Transposeview( dotBlock->get_op_rep(CRE, getSpinQuantum(k), k) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, *opCD, opD, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  // Sys has j,k;
  else if (sysBlock->get_op_array(CRE_CRE).has_local_index(j,k)) {
//pout << "maw sys(j,k)\n";
    assert( dotBlock->get_op_array(CRE).has_local_index(i) );
    build_pattern = "(C(DD))";
    const boost::shared_ptr<SparseMatrix>& opC = dotBlock->get_op_rep(CRE, getSpinQuantum(i), i);
    // We need to commute after transposing (j,k)
//FIXME minus sign if j!=k??
    Transposeview opDD = Transposeview( sysBlock->get_op_rep(CRE_CRE, quantum_ladder.at(build_pattern).at(0), j,k) );
    // Tensor product of dot*sys so need to take into account parity factors
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
//FIXME parity broken with opDD as we saw with CDC operator???  (c.f. minus signs in operator_wrappers...) ??
    double parity = getCommuteParity( opC->get_deltaQuantum(), opDD.get_deltaQuantum(), get_deltaQuantum() );
//    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *opC, opDD, &b, &(b.get_stateInfo()), *this, 1.0*parity);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *opC, opDD, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  // Dot has i,j
  else if (dotBlock->get_op_array(CRE_DES).has_local_index(i,j)) {
//pout << "maw dot(i,j)\n";
    assert( i == j );
    assert( sysBlock->get_op_array(CRE).has_local_index(k) );
    build_pattern = "((CD)D)";
    const boost::shared_ptr<SparseMatrix>& opCD = dotBlock->get_op_rep(CRE_DES, quantum_ladder.at(build_pattern).at(0), i,j);
    Transposeview opD = Transposeview( sysBlock->get_op_rep(CRE, getSpinQuantum(k), k) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    // Tensor product of dot*sys so need to take into account parity factors
    double parity = getCommuteParity( opCD->get_deltaQuantum(), opD.get_deltaQuantum(), get_deltaQuantum() );
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *opCD, opD, &b, &(b.get_stateInfo()), *this, 1.0*parity);
  }
  else assert(false);

  dmrginp.makeopsT -> stop();
//pout << "CDD op\n";
//pout << *this;
//pout << "done building CreDesDes renormalized operator!\n";

}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

double SpinAdapted::CreDesDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
//cout << "building CreDesDes explicitly from CSF..\n";
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
    return boost::shared_ptr<CreDesDes>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
  }
  else {
    boost::shared_ptr<SparseMatrix> rep(new CreDesDes);
    *rep = *this;
    rep->build(*block);

    return rep;
  }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Des,Cre)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinAdapted::CreDesCre::build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs)
{
//cout << "building CreDesCre renormalized operator on disk...\n";
//FIXME timer
//  dmrginp.makeopsT -> start();
  built = true;
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
//cout << "i,j,k = " << i << " " << j << " " << k << endl;

  SpinBlock* sysBlock = b.get_leftBlock();
  SpinBlock* dotBlock = b.get_rightBlock();

  // Sys has i,j,k
  if (sysBlock->get_op_array(CRE_DES_CRE).has_local_index(i,j,k)) {
//cout << "maw disk sys(i,j,k)\n";
    assert( sysBlock->get_op_array(CRE_DES_CRE).get_element(i,j,k).at(0)->get_built_on_disk() );
    // Retrieve from disk 3-index operator on sys block
    boost::shared_ptr<SparseMatrix> op (new CreDesCre);
    boost::archive::binary_iarchive load_op(sysfs);
    load_op >> *op;
    // Assume this is the operator we want, and expand it to the big block
    assert( i == op->get_orbs()[0] );
    assert( j == op->get_orbs()[1] );
    assert( k == op->get_orbs()[2] );
    // Build according to previous build_pattern
    std::string build_pattern_old = op->get_build_pattern();
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  // Dot has i,j,k
  else if (dotBlock->get_op_array(CRE_DES_CRE).has_local_index(i,j,k)) {
//cout << "maw dot(i,j,k)\n";
    assert( dotBlock->get_op_array(CRE_DES_CRE).get_element(i,j,k).at(0)->get_built_on_disk() );
    // Retrieve from disk 3-index operator on dot block
    boost::shared_ptr<SparseMatrix> op (new CreDesCre);
    boost::archive::binary_iarchive load_op(dotfs);
    load_op >> *op;
    // Assume this is the operator we want, and expand it to the big block
    assert( i == op->get_orbs()[0] );
    assert( j == op->get_orbs()[1] );
    assert( k == op->get_orbs()[2] );
    // Build according to previous build_pattern
    std::string build_pattern_old = op->get_build_pattern();
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  else {
    // Build from in-core 1 and 2-index operators
    build(b);
  }

//  dmrginp.makeopsT -> stop();

}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinAdapted::CreDesCre::build(const SpinBlock& b)
{
//cout << "building CreDesCre renormalized operator...\n";
  dmrginp.makeopsT -> start();
  built = true;
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
//pout << "indices  " << i << " " << j << " " << k << std::endl;

  SpinBlock* sysBlock = b.get_leftBlock();
  SpinBlock* dotBlock = b.get_rightBlock();
//pout << "sysBlock:\n";
//pout << *sysBlock;
//pout << "dotBlock:\n";
//pout << *dotBlock;

//  SpinQuantum deltaQuantum12 = get_quantum_ladder().at("0").at(0);

  // Sys has i,j,k
  if (sysBlock->get_op_array(CRE_DES_CRE).has_local_index(i,j,k)) {
//pout << "maw sys(i,j,k)\n";
    std::string build_pattern_old = sysBlock->get_op_array(CRE_DES_CRE).get_element(i,j,k).at(0)->get_build_pattern();
    assert( build_pattern_old == sysBlock->get_op_array(CRE_DES_CRE).get_element(i,j,k).at(1)->get_build_pattern() );
    assert( build_pattern_old == sysBlock->get_op_array(CRE_DES_CRE).get_element(i,j,k).at(2)->get_build_pattern() );
    const boost::shared_ptr<SparseMatrix>& op = sysBlock->get_op_rep(CRE_DES_CRE, quantum_ladder.at(build_pattern_old), i,j,k);
    // Build according to previous build_pattern
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  // Dot has i,j,k
  else if (dotBlock->get_op_array(CRE_DES_CRE).has_local_index(i,j,k)) {
//pout << "maw dot(i,j,k)\n";
    assert( i == j );
    assert( j == k );
    std::string build_pattern_old = dotBlock->get_op_array(CRE_DES_CRE).get_element(i,j,k).at(0)->get_build_pattern();
    assert( build_pattern_old == dotBlock->get_op_array(CRE_DES_CRE).get_element(i,j,k).at(1)->get_build_pattern() );
    assert( build_pattern_old == dotBlock->get_op_array(CRE_DES_CRE).get_element(i,j,k).at(2)->get_build_pattern() );
    const boost::shared_ptr<SparseMatrix>& op = dotBlock->get_op_rep(CRE_DES_CRE, quantum_ladder.at(build_pattern_old), i,j,k);
    // Build according to previous build_pattern
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  // Dot has j,k;
  else if (dotBlock->get_op_array(CRE_DES).has_local_index(j,k)) {
//pout << "maw dot(j,k)\n";
    assert( j == k );
    assert( sysBlock->get_op_array(CRE).has_local_index(i) );
    build_pattern = "(C(DC))";
    const boost::shared_ptr<SparseMatrix>& opC = sysBlock->get_op_rep(CRE, getSpinQuantum(i), i);
    const boost::shared_ptr<SparseMatrix>& opDC = dotBlock->get_op_rep(DES_CRE, quantum_ladder.at(build_pattern).at(0), j,k);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, *opC, *opDC, &b, &(b.get_stateInfo()), *this, 1.0);

//CANT USE TRANSPOSE!!!
//WRONG    Transposeview opDC = Transposeview( dotBlock->get_op_rep(CRE_DES, quantum_ladder.at(build_pattern).at(0), j,k) );
//    //FIXME  Tranpose OK since j=k (test with DES_CRE and see if same???)
//    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
//    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, *opC, opDC, &b, &(b.get_stateInfo()), *this, 1.0);

  }
  // Sys has i,j
  else if (sysBlock->get_op_array(CRE_DES).has_local_index(i,j)) {
//pout << "maw dot(k)\n";
    assert( dotBlock->get_op_array(CRE).has_local_index(k) );
    build_pattern = "((CD)C)";
    const boost::shared_ptr<SparseMatrix>& opCD = sysBlock->get_op_rep(CRE_DES, quantum_ladder.at(build_pattern).at(0), i,j);
    const boost::shared_ptr<SparseMatrix>& opC = dotBlock->get_op_rep(CRE, getSpinQuantum(k), k);
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, *opCD, *opC, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  // Sys has j,k;
  else if (sysBlock->get_op_array(CRE_DES).has_local_index(j,k)) {
//pout << "maw sys(j,k)\n";
    assert( dotBlock->get_op_array(CRE).has_local_index(i) );
    assert( sysBlock->get_op_array(DES_CRE).has_local_index(j,k) );
    build_pattern = "(C(DC))";
    const boost::shared_ptr<SparseMatrix>& opC = dotBlock->get_op_rep(CRE, getSpinQuantum(i), i);
    const boost::shared_ptr<SparseMatrix>& opDC = sysBlock->get_op_rep(DES_CRE, quantum_ladder.at(build_pattern).at(0), j,k);
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    double parity = getCommuteParity( opC->get_deltaQuantum(), opDC->get_deltaQuantum(), get_deltaQuantum() );
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *opC, *opDC, &b, &(b.get_stateInfo()), *this, 1.0*parity);

//CANT USE TRANSPOSE!!!
////    // We need to commute after transposing (j,k)
//////FIXME minus sign if j!=k??
////    Transposeview opDC = Transposeview( sysBlock->get_op_rep(CRE_DES, quantum_ladder.at(build_pattern).at(0), j,k) );
////    // Tensor product of dot*sys so need to take into account parity factors
////    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
////    double parity = getCommuteParity( opC->get_deltaQuantum(), opDC.get_deltaQuantum(), get_deltaQuantum() );
//////    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *opC, opDC, &b, &(b.get_stateInfo()), *this, 1.0*parity);
//////FIXME!!! WHY DON'T NEED PARITY????
////    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *opC, opDC, &b, &(b.get_stateInfo()), *this, 1.0);

  }
  // Dot has i,j
  else if (dotBlock->get_op_array(CRE_DES).has_local_index(i,j)) {
//pout << "maw dot(i,j)\n";
    assert( i == j );
    assert( sysBlock->get_op_array(CRE).has_local_index(k) );
    build_pattern = "((CD)C)";
    const boost::shared_ptr<SparseMatrix>& opCD = dotBlock->get_op_rep(CRE_DES, quantum_ladder.at(build_pattern).at(0), i,j);
    const boost::shared_ptr<SparseMatrix>& opC = sysBlock->get_op_rep(CRE, getSpinQuantum(k), k);
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    // Tensor product of dot*sys so need to take into account parity factors
    double parity = getCommuteParity( opCD->get_deltaQuantum(), opC->get_deltaQuantum(), get_deltaQuantum() );
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *opCD, *opC, &b, &(b.get_stateInfo()), *this, 1.0*parity);
  }
  else assert(false);

  dmrginp.makeopsT -> stop();
//pout << "CDC op\n";
//pout << *this;
//pout << "done building CreDesCre renormalized operator!\n";

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
////  assert( b.get_op_array(CRE).has_local_index(i) );
////  assert( b.get_op_array(CRE).has_local_index(j) );
////  assert( b.get_op_array(CRE).has_local_index(k) );
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
//cout << "building CreDesCre in CSF space explicitly..\n";
  assert( build_pattern == "((CD)C)" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];

  // Must take into account how the 3-index is built from a combination of 2-index and 1-index
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("((CD)C)");
  assert( quantum_ladder.size() == 2 );
  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  assert( deltaQuantum == deltaQuantum123 );

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
  TensorOp C3(K, 1); 

  // Combine first two operators
  TensorOp CD = C1.product(D2, spin12, irrep12);
  // Combine with third operator
  TensorOp CDC = CD.product(C3, spin123, irrep123);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum, ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CDC, ladder[i]) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CreDesCre::getworkingrepresentation(const SpinBlock* block)
{
  assert(this->get_initialised());
  if (this->get_built()) {
    return boost::shared_ptr<CreDesCre>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
  }
  else {
    boost::shared_ptr<SparseMatrix> rep(new CreDesCre);
    *rep = *this;
    rep->build(*block);
    return rep;
  }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Cre,Cre)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinAdapted::CreCreCre::build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs)
{
//cout << "building CreCreCre renormalized operator on disk...\n";
//FIXME timer
//  dmrginp.makeopsT -> start();
  built = true;
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
//cout << "i,j,k = " << i << " " << j << " " << k << endl;

  SpinBlock* sysBlock = b.get_leftBlock();
  SpinBlock* dotBlock = b.get_rightBlock();

  // Sys has i,j,k
  if (sysBlock->get_op_array(CRE_CRE_CRE).has_local_index(i,j,k)) {
//cout << "maw disk sys(i,j,k)\n";
    assert( sysBlock->get_op_array(CRE_CRE_CRE).get_element(i,j,k).at(0)->get_built_on_disk() );
    // Retrieve from disk 3-index operator on sys block
    boost::shared_ptr<SparseMatrix> op (new CreCreCre);
    boost::archive::binary_iarchive load_op(sysfs);
    load_op >> *op;
    // Assume this is the operator we want, and expand it to the big block
    assert( i == op->get_orbs()[0] );
    assert( j == op->get_orbs()[1] );
    assert( k == op->get_orbs()[2] );
    // Build according to previous build_pattern
    std::string build_pattern_old = op->get_build_pattern();
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  // Dot has i,j,k
  else if (dotBlock->get_op_array(CRE_CRE_CRE).has_local_index(i,j,k)) {
//cout << "maw dot(i,j,k)\n";
    assert( dotBlock->get_op_array(CRE_CRE_CRE).get_element(i,j,k).at(0)->get_built_on_disk() );
    // Retrieve from disk 3-index operator on dot block
    boost::shared_ptr<SparseMatrix> op (new CreCreCre);
    boost::archive::binary_iarchive load_op(dotfs);
    load_op >> *op;
    // Assume this is the operator we want, and expand it to the big block
    assert( i == op->get_orbs()[0] );
    assert( j == op->get_orbs()[1] );
    assert( k == op->get_orbs()[2] );
    // Build according to previous build_pattern
    std::string build_pattern_old = op->get_build_pattern();
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  else {
    // Build from in-core 1 and 2-index operators
    build(b);
  }

//  dmrginp.makeopsT -> stop();

}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinAdapted::CreCreCre::build(const SpinBlock& b)
{
//cout << "building CreCreCre renormalized operator...\n";
  dmrginp.makeopsT -> start();
  built = true;
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
//cout << "i,j,k = " << i << " " << j << " " << k << endl;

  SpinBlock* sysBlock = b.get_leftBlock();
  SpinBlock* dotBlock = b.get_rightBlock();

  // Sys has i,j,k
  if (sysBlock->get_op_array(CRE_CRE_CRE).has_local_index(i,j,k)) {
//cout << "maw sys(i,j,k)\n";
    std::string build_pattern_old = sysBlock->get_op_array(CRE_CRE_CRE).get_element(i,j,k).at(0)->get_build_pattern();
    assert( build_pattern_old == sysBlock->get_op_array(CRE_CRE_CRE).get_element(i,j,k).at(1)->get_build_pattern() );
    assert( build_pattern_old == sysBlock->get_op_array(CRE_CRE_CRE).get_element(i,j,k).at(2)->get_build_pattern() );
    const boost::shared_ptr<SparseMatrix>& op = sysBlock->get_op_rep(CRE_CRE_CRE, quantum_ladder.at(build_pattern_old), i,j,k);
    // Build according to previous build_pattern
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  // Dot has i,j,k
  else if (dotBlock->get_op_array(CRE_CRE_CRE).has_local_index(i,j,k)) {
//cout << "maw dot(i,j,k)\n";
    assert( i == j );
    assert( j == k );
    std::string build_pattern_old = dotBlock->get_op_array(CRE_CRE_CRE).get_element(i,j,k).at(0)->get_build_pattern();
    assert( build_pattern_old == dotBlock->get_op_array(CRE_CRE_CRE).get_element(i,j,k).at(1)->get_build_pattern() );
    assert( build_pattern_old == dotBlock->get_op_array(CRE_CRE_CRE).get_element(i,j,k).at(2)->get_build_pattern() );
    const boost::shared_ptr<SparseMatrix>& op = dotBlock->get_op_rep(CRE_CRE_CRE, quantum_ladder.at(build_pattern_old), i,j,k);
    // Build according to previous build_pattern
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  // Dot has j,k;
  else if (dotBlock->get_op_array(CRE_CRE).has_local_index(j,k)) {
//cout << "maw dot(j,k)\n";
    assert( j == k );
    assert( sysBlock->get_op_array(CRE).has_local_index(i) );
    build_pattern = "(C(CC))";
    const boost::shared_ptr<SparseMatrix>& opC  = sysBlock->get_op_rep(CRE, getSpinQuantum(i), i);
    const boost::shared_ptr<SparseMatrix>& opCC = dotBlock->get_op_rep(CRE_CRE, quantum_ladder.at(build_pattern).at(0), j,k);
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, *opC, *opCC, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  // Sys has i,j
  else if (sysBlock->get_op_array(CRE_CRE).has_local_index(i,j)) {
//cout << "maw dot(k)\n";
    assert( dotBlock->get_op_array(CRE).has_local_index(k) );
    build_pattern = "((CC)C)";
    const boost::shared_ptr<SparseMatrix>& opCC = sysBlock->get_op_rep(CRE_CRE, quantum_ladder.at(build_pattern).at(0), i,j);
    const boost::shared_ptr<SparseMatrix>& opC  = dotBlock->get_op_rep(CRE, getSpinQuantum(k), k);
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, *opCC, *opC, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  // Sys has j,k;
  else if (sysBlock->get_op_array(CRE_CRE).has_local_index(j,k)) {
//cout << "maw sys(j,k)\n";
    assert( dotBlock->get_op_array(CRE).has_local_index(i) );
    build_pattern = "(C(CC))";
    const boost::shared_ptr<SparseMatrix>& opC  = dotBlock->get_op_rep(CRE, getSpinQuantum(i), i);
    const boost::shared_ptr<SparseMatrix>& opCC = sysBlock->get_op_rep(CRE_CRE, quantum_ladder.at(build_pattern).at(0), j,k);
    // Tensor product of dot*sys so need to take into account parity factors
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
//FIXME parity broken ??
    double parity = getCommuteParity( opC->get_deltaQuantum(), opCC->get_deltaQuantum(), get_deltaQuantum() );
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *opC, *opCC, &b, &(b.get_stateInfo()), *this, 1.0*parity);
  }
  // Dot has i,j
  else if (dotBlock->get_op_array(CRE_CRE).has_local_index(i,j)) {
//cout << "maw dot(i,j)\n";
    assert( i == j );
    assert( sysBlock->get_op_array(CRE).has_local_index(k) );
    build_pattern = "((CC)C)";
    const boost::shared_ptr<SparseMatrix>& opCC = dotBlock->get_op_rep(CRE_CRE, quantum_ladder.at(build_pattern).at(0), i,j);
    const boost::shared_ptr<SparseMatrix>& opC  = sysBlock->get_op_rep(CRE, getSpinQuantum(k), k);
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    // Tensor product of dot*sys so need to take into account parity factors
//FIXME parity broken ??
    double parity = getCommuteParity( opCC->get_deltaQuantum(), opC->get_deltaQuantum(), get_deltaQuantum() );
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *opCC, *opC, &b, &(b.get_stateInfo()), *this, 1.0*parity);
  }
  else assert(false);

  dmrginp.makeopsT -> stop();
//pout << "CCC op\n";
//pout << *this;
//pout << "done building CreCreCre renormalized operator!\n";

}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

double SpinAdapted::CreCreCre::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
//cout << "building CreCreCre explicitly from CSF..\n";
//cout << "mpirank = " << mpigetrank() << endl;
  assert( build_pattern == "((CC)C)" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
//cout << "i,j,k = " << I << " " << J << " " << K << endl;

  // Must take into account how the 3-index is built from a combination of 2-index and 1-index
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("((CC)C)");
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

  TensorOp C1(I,1); 
  TensorOp C2(J,1); 
  TensorOp C3(K,1); 

  // Combine first two operators
//FIXME I=J argument
  TensorOp CC = C1.product(C2, spin12, irrep12);
  // Combine with third operator
  TensorOp CCC = CC.product(C3, spin123, irrep123);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum, ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CCC, ladder[i]) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CreCreCre::getworkingrepresentation(const SpinBlock* block)
{
//pout << "CreCreCre::getworkingrepresentation\n";
  assert(this->get_initialised());
  if (this->get_built()) {
//pout << "get CCC from memory\n";
    return boost::shared_ptr<CreCreCre>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
  }
  else {
//pout << "build CCC\n";
    boost::shared_ptr<SparseMatrix> rep(new CreCreCre);
    *rep = *this;
    rep->build(*block);

    return rep;
  }
}

//===========================================================================================================================================================
// 4PDM operators
//===========================================================================================================================================================

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Des,Cre,Des)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinAdapted::DesCreDes::build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs)
{
//cout << "building DesCreDes renormalized operator on disk...\n";
//FIXME timer
//  dmrginp.makeopsT -> start();
  built = true;
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
//cout << "i,j,k = " << i << " " << j << " " << k << endl;

  SpinBlock* sysBlock = b.get_leftBlock();
  SpinBlock* dotBlock = b.get_rightBlock();

  // Sys has i,j,k
  if (sysBlock->get_op_array(DES_CRE_DES).has_local_index(i,j,k)) {
//cout << "maw disk sys(i,j,k)\n";
    assert( sysBlock->get_op_array(DES_CRE_DES).get_element(i,j,k).at(0)->get_built_on_disk() );
    // Retrieve from disk 3-index operator on sys block
    boost::shared_ptr<SparseMatrix> op (new DesCreDes);
    boost::archive::binary_iarchive load_op(sysfs);
    load_op >> *op;
    // Assume this is the operator we want, and expand it to the big block
    assert( i == op->get_orbs()[0] );
    assert( j == op->get_orbs()[1] );
    assert( k == op->get_orbs()[2] );
    // Build according to previous build_pattern
    std::string build_pattern_old = op->get_build_pattern();
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  // Dot has i,j,k
  else if (dotBlock->get_op_array(DES_CRE_DES).has_local_index(i,j,k)) {
//cout << "maw dot(i,j,k)\n";
    assert( dotBlock->get_op_array(DES_CRE_DES).get_element(i,j,k).at(0)->get_built_on_disk() );
    // Retrieve from disk 3-index operator on dot block
    boost::shared_ptr<SparseMatrix> op (new DesCreDes);
    boost::archive::binary_iarchive load_op(dotfs);
    load_op >> *op;
    // Assume this is the operator we want, and expand it to the big block
    assert( i == op->get_orbs()[0] );
    assert( j == op->get_orbs()[1] );
    assert( k == op->get_orbs()[2] );
    // Build according to previous build_pattern
    std::string build_pattern_old = op->get_build_pattern();
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  else {
    // Build from in-core 1 and 2-index operators
    build(b);
  }

//  dmrginp.makeopsT -> stop();

}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinAdapted::DesCreDes::build(const SpinBlock& b)
{
//cout << "building DesCreDes renormalized operator...\n";
  dmrginp.makeopsT -> start();
  built = true;
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];

  SpinBlock* sysBlock = b.get_leftBlock();
  SpinBlock* dotBlock = b.get_rightBlock();

  // Sys has i,j,k
  if (sysBlock->get_op_array(DES_CRE_DES).has_local_index(i,j,k)) {
    std::string build_pattern_old = sysBlock->get_op_array(DES_CRE_DES).get_element(i,j,k).at(0)->get_build_pattern();
    assert( build_pattern_old == sysBlock->get_op_array(DES_CRE_DES).get_element(i,j,k).at(1)->get_build_pattern() );
    assert( build_pattern_old == sysBlock->get_op_array(DES_CRE_DES).get_element(i,j,k).at(2)->get_build_pattern() );
    const boost::shared_ptr<SparseMatrix>& op = sysBlock->get_op_rep(DES_CRE_DES, quantum_ladder.at(build_pattern_old), i,j,k);
    // Build according to previous build_pattern
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  // Dot has i,j,k
  else if (dotBlock->get_op_array(DES_CRE_DES).has_local_index(i,j,k)) {
    assert( i == j );
    assert( j == k );
    std::string build_pattern_old = dotBlock->get_op_array(DES_CRE_DES).get_element(i,j,k).at(0)->get_build_pattern();
    assert( build_pattern_old == dotBlock->get_op_array(DES_CRE_DES).get_element(i,j,k).at(1)->get_build_pattern() );
    assert( build_pattern_old == dotBlock->get_op_array(DES_CRE_DES).get_element(i,j,k).at(2)->get_build_pattern() );
    const boost::shared_ptr<SparseMatrix>& op = dotBlock->get_op_rep(DES_CRE_DES, quantum_ladder.at(build_pattern_old), i,j,k);
    // Build according to previous build_pattern
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  // Dot has j,k;
  else if (dotBlock->get_op_array(CRE_DES).has_local_index(j,k)) {
    assert( j == k );
    assert( sysBlock->get_op_array(CRE).has_local_index(i) );
    build_pattern = "(D(CD))";
    Transposeview opD = Transposeview( sysBlock->get_op_rep(CRE, getSpinQuantum(i), i) );
    const boost::shared_ptr<SparseMatrix>& opCD = dotBlock->get_op_rep(CRE_DES, quantum_ladder.at(build_pattern).at(0), j,k);
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, opD, *opCD, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  // Sys has i,j
  else if (sysBlock->get_op_array(DES_CRE).has_local_index(i,j)) {
    assert( dotBlock->get_op_array(CRE).has_local_index(k) );
    build_pattern = "((DC)D)";
    const boost::shared_ptr<SparseMatrix>& opDC = sysBlock->get_op_rep(DES_CRE, quantum_ladder.at(build_pattern).at(0), i,j);
    Transposeview opD = Transposeview( dotBlock->get_op_rep(CRE, getSpinQuantum(k), k) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, *opDC, opD, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  // Sys has j,k;
  else if (sysBlock->get_op_array(CRE_DES).has_local_index(j,k)) {
    assert( dotBlock->get_op_array(CRE).has_local_index(i) );
    build_pattern = "(D(CD))";
    Transposeview opD = Transposeview( dotBlock->get_op_rep(CRE, getSpinQuantum(i), i) );
    const boost::shared_ptr<SparseMatrix>& opCD = sysBlock->get_op_rep(CRE_DES, quantum_ladder.at(build_pattern).at(0), j,k);
    // Tensor product of dot*sys so need to take into account parity factors
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
//FIXME parity broken ??
    double parity = getCommuteParity( opD.get_deltaQuantum(), opCD->get_deltaQuantum(), get_deltaQuantum() );
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, opD, *opCD, &b, &(b.get_stateInfo()), *this, 1.0*parity);
  }
  // Dot has i,j
  else if (dotBlock->get_op_array(DES_CRE).has_local_index(i,j)) {
    assert( i == j );
    assert( sysBlock->get_op_array(CRE).has_local_index(k) );
    build_pattern = "((DC)D)";
    const boost::shared_ptr<SparseMatrix>& opDC = dotBlock->get_op_rep(DES_CRE, quantum_ladder.at(build_pattern).at(0), i,j);
    Transposeview opD = Transposeview( sysBlock->get_op_rep(CRE, getSpinQuantum(k), k) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    // Tensor product of dot*sys so need to take into account parity factors
//FIXME parity broken ??
    double parity = getCommuteParity( opDC->get_deltaQuantum(), opD.get_deltaQuantum(), get_deltaQuantum() );
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *opDC, opD, &b, &(b.get_stateInfo()), *this, 1.0*parity);
  }
  else assert(false);

  dmrginp.makeopsT -> stop();

}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

double SpinAdapted::DesCreDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  assert( build_pattern == "((DC)D)" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];

  // Must take into account how the 3-index is built from a combination of 2-index and 1-index
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("((DC)D)");
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

  TensorOp D1(I,-1); 
  TensorOp C2(J, 1); 
  TensorOp D3(K,-1); 

  // Combine first two operators
  TensorOp DC = D1.product(C2, spin12, irrep12);
  // Combine with third operator
  TensorOp DCD = DC.product(D3, spin123, irrep123);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum, ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, DCD, ladder[i]) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::DesCreDes::getworkingrepresentation(const SpinBlock* block)
{
  assert(this->get_initialised());
  if (this->get_built()) {
    return boost::shared_ptr<DesCreDes>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
  }
  else {
    boost::shared_ptr<SparseMatrix> rep(new DesCreDes);
    *rep = *this;
    rep->build(*block);

    return rep;
  }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Des,Des,Cre)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinAdapted::DesDesCre::build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs)
{
//cout << "building DesDesCre renormalized operator on disk...\n";
//FIXME timer
//  dmrginp.makeopsT -> start();
  built = true;
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
//cout << "i,j,k = " << i << " " << j << " " << k << endl;

  SpinBlock* sysBlock = b.get_leftBlock();
  SpinBlock* dotBlock = b.get_rightBlock();

  // Sys has i,j,k
  if (sysBlock->get_op_array(DES_DES_CRE).has_local_index(i,j,k)) {
//cout << "maw disk sys(i,j,k)\n";
    assert( sysBlock->get_op_array(DES_DES_CRE).get_element(i,j,k).at(0)->get_built_on_disk() );
    // Retrieve from disk 3-index operator on sys block
    boost::shared_ptr<SparseMatrix> op (new DesDesCre);
    boost::archive::binary_iarchive load_op(sysfs);
    load_op >> *op;
    // Assume this is the operator we want, and expand it to the big block
    assert( i == op->get_orbs()[0] );
    assert( j == op->get_orbs()[1] );
    assert( k == op->get_orbs()[2] );
    // Build according to previous build_pattern
    std::string build_pattern_old = op->get_build_pattern();
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  // Dot has i,j,k
  else if (dotBlock->get_op_array(DES_DES_CRE).has_local_index(i,j,k)) {
//cout << "maw dot(i,j,k)\n";
    assert( dotBlock->get_op_array(DES_DES_CRE).get_element(i,j,k).at(0)->get_built_on_disk() );
    // Retrieve from disk 3-index operator on dot block
    boost::shared_ptr<SparseMatrix> op (new DesDesCre);
    boost::archive::binary_iarchive load_op(dotfs);
    load_op >> *op;
    // Assume this is the operator we want, and expand it to the big block
    assert( i == op->get_orbs()[0] );
    assert( j == op->get_orbs()[1] );
    assert( k == op->get_orbs()[2] );
    // Build according to previous build_pattern
    std::string build_pattern_old = op->get_build_pattern();
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  else {
    // Build from in-core 1 and 2-index operators
    build(b);
  }

//  dmrginp.makeopsT -> stop();

}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinAdapted::DesDesCre::build(const SpinBlock& b)
{
//cout << "building DesDesCre renormalized operator...\n";
  dmrginp.makeopsT -> start();
  built = true;
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];

  SpinBlock* sysBlock = b.get_leftBlock();
  SpinBlock* dotBlock = b.get_rightBlock();

  // Sys has i,j,k
  if (sysBlock->get_op_array(DES_DES_CRE).has_local_index(i,j,k)) {
//pout << "maw sys(i,j,k)\n";
    std::string build_pattern_old = sysBlock->get_op_array(DES_DES_CRE).get_element(i,j,k).at(0)->get_build_pattern();
    assert( build_pattern_old == sysBlock->get_op_array(DES_DES_CRE).get_element(i,j,k).at(1)->get_build_pattern() );
    assert( build_pattern_old == sysBlock->get_op_array(DES_DES_CRE).get_element(i,j,k).at(2)->get_build_pattern() );
    const boost::shared_ptr<SparseMatrix>& op = sysBlock->get_op_rep(DES_DES_CRE, quantum_ladder.at(build_pattern_old), i,j,k);
    // Build according to previous build_pattern
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  // Dot has i,j,k
  else if (dotBlock->get_op_array(DES_DES_CRE).has_local_index(i,j,k)) {
//pout << "maw dot(i,j,k)\n";
    assert( i == j );
    assert( j == k );
    std::string build_pattern_old = dotBlock->get_op_array(DES_DES_CRE).get_element(i,j,k).at(0)->get_build_pattern();
    assert( build_pattern_old == dotBlock->get_op_array(DES_DES_CRE).get_element(i,j,k).at(1)->get_build_pattern() );
    assert( build_pattern_old == dotBlock->get_op_array(DES_DES_CRE).get_element(i,j,k).at(2)->get_build_pattern() );
    const boost::shared_ptr<SparseMatrix>& op = dotBlock->get_op_rep(DES_DES_CRE, quantum_ladder.at(build_pattern_old), i,j,k);
    // Build according to previous build_pattern
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  // Dot has j,k;
  else if (dotBlock->get_op_array(DES_CRE).has_local_index(j,k)) {
//pout << "maw dot(j,k)\n";
    assert( j == k );
    assert( sysBlock->get_op_array(CRE).has_local_index(i) );
    build_pattern = "(D(DC))";
    Transposeview opD = Transposeview( sysBlock->get_op_rep(CRE, getSpinQuantum(i), i) );
    const boost::shared_ptr<SparseMatrix>& opDC = dotBlock->get_op_rep(DES_CRE, quantum_ladder.at(build_pattern).at(0), j,k);
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, opD, *opDC, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  // Sys has i,j
  else if (sysBlock->get_op_array(CRE_CRE).has_local_index(i,j)) {
//pout << "maw sys(i,j)\n";
    assert( dotBlock->get_op_array(CRE).has_local_index(k) );
    build_pattern = "((DD)C)";
//FIXME transpose indices / commute ?? minus sign if i!=j ??  c.f. CDD case above?
    Transposeview opDD = Transposeview( sysBlock->get_op_rep(CRE_CRE, quantum_ladder.at(build_pattern).at(0), i,j) );
    const boost::shared_ptr<SparseMatrix>& opC  = dotBlock->get_op_rep(CRE, getSpinQuantum(k), k);
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, opDD, *opC, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  // Sys has j,k;
  else if (sysBlock->get_op_array(DES_CRE).has_local_index(j,k)) {
//pout << "maw sys(j,k)\n";
    assert( dotBlock->get_op_array(CRE).has_local_index(i) );
    build_pattern = "(D(DC))";
    Transposeview opD = Transposeview( dotBlock->get_op_rep(CRE, getSpinQuantum(i), i) );
    const boost::shared_ptr<SparseMatrix>& opDC = sysBlock->get_op_rep(DES_CRE, quantum_ladder.at(build_pattern).at(0), j,k);
    // Tensor product of dot*sys so need to take into account parity factors
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
//FIXME parity broken ??
    double parity = getCommuteParity( opD.get_deltaQuantum(), opDC->get_deltaQuantum(), get_deltaQuantum() );
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, opD, *opDC, &b, &(b.get_stateInfo()), *this, 1.0*parity);
  }
  // Dot has i,j
  else if (dotBlock->get_op_array(CRE_CRE).has_local_index(i,j)) {
//pout << "maw dot(i,j)\n";
    assert( i == j );
    assert( sysBlock->get_op_array(CRE).has_local_index(k) );
    build_pattern = "((DD)C)";
    Transposeview opDD = Transposeview( dotBlock->get_op_rep(CRE_CRE, quantum_ladder.at(build_pattern).at(0), i,j) );
//FIXME transpose indices / commute ?? minus sign if i!=j ??  c.f. CDD case above?
    const boost::shared_ptr<SparseMatrix>& opC  = sysBlock->get_op_rep(CRE, getSpinQuantum(k), k);
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    // Tensor product of dot*sys so need to take into account parity factors
//FIXME parity broken ??
    double parity = getCommuteParity( opDD.get_deltaQuantum(), opC->get_deltaQuantum(), get_deltaQuantum() );
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, opDD, *opC, &b, &(b.get_stateInfo()), *this, 1.0*parity);
  }
  else assert(false);

  dmrginp.makeopsT -> stop();

}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

double SpinAdapted::DesDesCre::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  assert( build_pattern == "((DD)C)" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];

  // Must take into account how the 3-index is built from a combination of 2-index and 1-index
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("((DD)C)");
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

  TensorOp D1(I,-1); 
  TensorOp D2(J,-1); 
  TensorOp C3(K, 1); 

  // Combine first two operators
//FIXME I=J argument??
  TensorOp DD = D1.product(D2, spin12, irrep12);
  // Combine with third operator
  TensorOp DDC = DD.product(C3, spin123, irrep123);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum, ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, DDC, ladder[i]) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::DesDesCre::getworkingrepresentation(const SpinBlock* block)
{
  assert(this->get_initialised());
  if (this->get_built()) {
    return boost::shared_ptr<DesDesCre>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
  }
  else {
    boost::shared_ptr<SparseMatrix> rep(new DesDesCre);
    *rep = *this;
    rep->build(*block);

    return rep;
  }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

