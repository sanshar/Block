/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//
// Extension of Operators.C for 4-index operators
//FIXME there's a lot of duplication, especially in build_from_disk... Templates??
//
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "op_components.h"
#include "BaseOperator.h"
#include "spinblock.h"
#include "operatorfunctions.h"
#include "tensor_operator.h"
#include "four_index_ops.h"

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Cre,Des,Des)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
// Note it is crucial to allocate the new operator only after deltaQuantum is well-defined.

void SpinAdapted::CreCreDesDes::build(const SpinBlock& b)
{
cout << "building CreCreDesDes renormalized operator...\n";
  dmrginp.makeopsT -> start();
  built = true;
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
  const int l = get_orbs()[3];
cout << "indices = " << i << "," << j << "," << k << "," << l << endl;

  SpinBlock* sysBlock = b.get_leftBlock();
  SpinBlock* dotBlock = b.get_rightBlock();

  // Sys has i,j,k,l
  if (sysBlock->get_op_array(CRE_CRE_DES_DES).has_local_index(i,j,k,l)) {
pout << "maw sys(i,j,k,l)\n";
    std::string build_pattern_old = sysBlock->get_op_array(CRE_CRE_DES_DES).get_element(i,j,k,l).at(0)->get_build_pattern();
cout << "build pattern = " << build_pattern_old << endl;
    const boost::shared_ptr<SparseMatrix>& op = sysBlock->get_op_rep(CRE_CRE_DES_DES, quantum_ladder.at(build_pattern_old), i,j,k,l);
    // Build according to previous build_pattern
    build_pattern = build_pattern_old;
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(2);
//cout << "maw op before TensorTrace\n";
//cout << *op;
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *op, &b, &(b.get_stateInfo()), *this);
//cout << "maw op after TensorTrace\n";
//cout << *this;
  }
  // Dot has i,j,k,l
  else if (dotBlock->get_op_array(CRE_CRE_DES_DES).has_local_index(i,j,k,l)) {
pout << "maw dot(i,j,k,l)\n";
    assert( i == j );
    assert( j == k );
    assert( k == l );
    std::string build_pattern_old = dotBlock->get_op_array(CRE_CRE_DES_DES).get_element(i,j,k,l).at(0)->get_build_pattern();
cout << "build pattern = " << build_pattern_old << endl;
    const boost::shared_ptr<SparseMatrix>& op = dotBlock->get_op_rep(CRE_CRE_DES_DES, quantum_ladder.at(build_pattern_old), i,j,k,l);
    // Build according to previous build_pattern
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(0) == op->get_quantum_ladder().at( build_pattern ).at(0) );
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    assert( get_quantum_ladder().at( build_pattern ).at(2) == op->get_quantum_ladder().at( build_pattern ).at(2) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(2);
//cout << "maw op before TensorTrace\n";
//cout << *op;
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *op, &b, &(b.get_stateInfo()), *this);
//cout << "maw op after TensorTrace\n";
//cout << *this;
//pout << "done dot(i,j,k,l)\n";
  }
  // Dot has j,k,l;
  else if (dotBlock->get_op_array(CRE_DES_DES).has_local_index(j,k,l)) {
pout << "maw dot(j,k,l)\n";
    assert( j == k );
    assert( k == l );
    assert( sysBlock->get_op_array(CRE).has_local_index(i) );
    // 3-index op
    std::string build_234 = dotBlock->get_op_array(CRE_DES_DES).get_element(j,k,l).at(0)->get_build_pattern();
    build_pattern = "((C)" + build_234 + ")";
cout << "build pattern = " << build_pattern << endl;
    std::vector<SpinQuantum> spin_234 = { quantum_ladder.at(build_pattern).at(0), quantum_ladder.at(build_pattern).at(1) };
    const boost::shared_ptr<SparseMatrix>& op234 = dotBlock->get_op_rep(CRE_DES_DES, spin_234, j,k,l);
    // 1-index op
    const boost::shared_ptr<SparseMatrix>& op1 = sysBlock->get_op_rep(CRE, getSpinQuantum(i), i);
    // 4-index op
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(2); //FIXME .at(2) brittle
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, *op1, *op234, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  // Sys has i,j,k
  else if (sysBlock->get_op_array(CRE_CRE_DES).has_local_index(i,j,k)) {
pout << "maw sys(i,j,k)\n";
    assert( dotBlock->get_op_array(CRE).has_local_index(l) );
    // 3-index op
    std::string build_123 = sysBlock->get_op_array(CRE_CRE_DES).get_element(i,j,k).at(0)->get_build_pattern();
    build_pattern = "(" + build_123 + "(D))";
cout << "build pattern = " << build_pattern << endl;
    std::vector<SpinQuantum> spin_123 = { quantum_ladder.at(build_pattern).at(0), quantum_ladder.at(build_pattern).at(1) };
    const boost::shared_ptr<SparseMatrix>& op123 = sysBlock->get_op_rep(CRE_CRE_DES, spin_123, i,j,k);
    // 1-index op
    Transposeview op4 = Transposeview( dotBlock->get_op_rep(CRE, getSpinQuantum(l), l) );
    // 4-index op
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(2);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, *op123, op4, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  // Sys has j,k,l;
  else if (sysBlock->get_op_array(CRE_DES_DES).has_local_index(j,k,l)) {
pout << "maw sys(j,k,l)\n";
    assert( dotBlock->get_op_array(CRE).has_local_index(i) );
    // 3-index op
    std::string build_234 = sysBlock->get_op_array(CRE_DES_DES).get_element(j,k,l).at(0)->get_build_pattern();
    build_pattern = "((C)" + build_234 + ")";
cout << "build pattern = " << build_pattern << endl;
    std::vector<SpinQuantum> spin_234 = { quantum_ladder.at(build_pattern).at(0), quantum_ladder.at(build_pattern).at(1) };
    const boost::shared_ptr<SparseMatrix>& op234 = sysBlock->get_op_rep(CRE_DES_DES, spin_234, j,k,l);
    // 1-index op
    const boost::shared_ptr<SparseMatrix>& op1 = dotBlock->get_op_rep(CRE, getSpinQuantum(i), i);
    // 4-index op
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(2);
    // Tensor product of dot*sys so need to take into account parity factors
    double parity = getCommuteParity( op1->get_deltaQuantum(), op234->get_deltaQuantum(), get_deltaQuantum() );
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *op1, *op234, &b, &(b.get_stateInfo()), *this, 1.0*parity);
  }
  // Dot has i,j,k
  else if (dotBlock->get_op_array(CRE_CRE_DES).has_local_index(i,j,k)) {
pout << "maw dot(i,j,k)\n";
    assert( i == j );
    assert( j == k );
    assert( sysBlock->get_op_array(CRE).has_local_index(l) );
    // 3-index op
    std::string build_123 = dotBlock->get_op_array(CRE_CRE_DES).get_element(i,j,k).at(0)->get_build_pattern();
    build_pattern = "(" + build_123 + "(D))";
cout << "build pattern = " << build_pattern << endl;
    std::vector<SpinQuantum> spin_123 = { quantum_ladder.at(build_pattern).at(0), quantum_ladder.at(build_pattern).at(1) };
    const boost::shared_ptr<SparseMatrix>& op123 = dotBlock->get_op_rep(CRE_CRE_DES, spin_123, i,j,k);
    // 1-index op
    Transposeview op4 = Transposeview( sysBlock->get_op_rep(CRE, getSpinQuantum(l), l) );
    // 4-index op
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(2);
    // Tensor product of dot*sys so need to take into account parity factors
    double parity = getCommuteParity( op123->get_deltaQuantum(), op4.get_deltaQuantum(), get_deltaQuantum() );
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *op123, op4, &b, &(b.get_stateInfo()), *this, 1.0*parity);
  }
  // Sys has i,j Dot has k,l;  
  else if (dotBlock->get_op_array(CRE_CRE).has_local_index(k,l)) {
pout << "maw sys(i,j); dot(k,l)\n";
    assert( sysBlock->get_op_array(CRE_CRE).has_local_index(i,j) );
    build_pattern = "((CC)(DD))";
cout << "build pattern = " << build_pattern << endl;
    Transposeview op34 = Transposeview( dotBlock->get_op_rep(CRE_CRE, quantum_ladder.at(build_pattern).at(1), k,l) );
    const boost::shared_ptr<SparseMatrix>& op12 = sysBlock->get_op_rep(CRE_CRE, quantum_ladder.at(build_pattern).at(0), i,j);
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(2);
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, *op12, op34, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  // Dot has i,j Sys has k,l;  
  else if (sysBlock->get_op_array(CRE_CRE).has_local_index(k,l)) {
pout << "maw dot(i,j); sys(k,l)\n";
    assert( dotBlock->get_op_array(CRE_CRE).has_local_index(i,j) );
    build_pattern = "((CC)(DD))";
cout << "build pattern = " << build_pattern << endl;
    Transposeview op34 = Transposeview( sysBlock->get_op_rep(CRE_CRE, quantum_ladder.at(build_pattern).at(1), k,l) );
    const boost::shared_ptr<SparseMatrix>& op12 = dotBlock->get_op_rep(CRE_CRE, quantum_ladder.at(build_pattern).at(0), i,j);
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(2);
cout << "Spin12 = " << get_quantum_ladder().at( build_pattern ).at(0).get_s() << endl;
cout << "Spin34 = " << get_quantum_ladder().at( build_pattern ).at(1).get_s() << endl;
cout << "totalSpin = " << get_deltaQuantum().get_s() << endl;
    // Tensor product of dot*sys so need to take into account parity factors
    double parity = getCommuteParity( op12->get_deltaQuantum(), op34.get_deltaQuantum(), get_deltaQuantum() );
    allocate(b.get_stateInfo());
    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *op12, op34, &b, &(b.get_stateInfo()), *this, 1.0*parity);
cout << "maw op after TensorProduct\n";
cout << *this;
  }
  else assert(false);

  dmrginp.makeopsT -> stop();

}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

double SpinAdapted::CreCreDesDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
////  assert( build_pattern == "((CC)(DD))" );
////  double element = 0.0;
////  int I = get_orbs()[0]; 
////  int J = get_orbs()[1];
////  int K = get_orbs()[2];
////  int L = get_orbs()[3];
////
////  // Must take into account how the 4-index is built from a combination of the 2-index ops
////  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("((CC)(DD))");
////  assert( quantum_ladder.size() == 3 );
////
////  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
////  SpinQuantum deltaQuantum34 = quantum_ladder.at(1);
////  SpinQuantum deltaQuantum1234 = quantum_ladder.at(2);
////  deltaQuantum = deltaQuantum1234;
////
////  // Spin quantum data for CC
////  IrrepSpace sym12 = deltaQuantum12.get_symm();
////  int irrep12 = deltaQuantum12.get_symm().getirrep();
////  int spin12 = deltaQuantum12.get_s();
////  // Spin quantum data for DD
////  IrrepSpace sym34 = deltaQuantum34.get_symm();
////  int irrep34 = deltaQuantum34.get_symm().getirrep();
////  int spin34 = deltaQuantum34.get_s();
////  // Spin quantum data for total operator
////  IrrepSpace sym1234 = deltaQuantum1234.get_symm();
////  int irrep1234 = deltaQuantum1234.get_symm().getirrep();
////  int spin1234 = deltaQuantum1234.get_s();
////
////  TensorOp C1(I, 1); 
////  TensorOp C2(J, 1); 
////  TensorOp D3(K,-1); 
////  TensorOp D4(L,-1); 
////
////  TensorOp CC = C1.product(C2, spin12, irrep12);   // I=J argument??
////  TensorOp DD = D3.product(D4, spin34, irrep34);
////  TensorOp CCDD = CC.product(DD, spin1234, irrep1234);


  assert( build_pattern == "(((CC)D)(D))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
  int L = get_orbs()[3];

  // Must take into account how the 4-index is built from a combination of the 2-index ops
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("(((CC)D)(D))");
  assert( quantum_ladder.size() == 3 );

  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  SpinQuantum deltaQuantum1234 = quantum_ladder.at(2);
  deltaQuantum = deltaQuantum1234;

  // Spin quantum data for CC
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s();
  // Spin quantum data for (CC)D
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123= deltaQuantum123.get_s();
  // Spin quantum data for total operator
  IrrepSpace sym1234 = deltaQuantum1234.get_symm();
  int irrep1234 = deltaQuantum1234.get_symm().getirrep();
  int spin1234 = deltaQuantum1234.get_s();

  TensorOp C1(I, 1); 
  TensorOp C2(J, 1); 
  TensorOp D3(K,-1); 
  TensorOp D4(L,-1); 

  TensorOp CC = C1.product(C2, spin12, irrep12);
  TensorOp CCD = CC.product(D3, spin123, irrep123);
  TensorOp CCDD = CCD.product(D4, spin1234, irrep1234);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum, ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CCDD, ladder[i]) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CreCreDesDes::getworkingrepresentation(const SpinBlock* block)
{
  assert(this->get_initialised());
  if (this->get_built()) {
    return boost::shared_ptr<CreCreDesDes>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
  }
  else {
    boost::shared_ptr<SparseMatrix> rep(new CreCreDesDes);
    *rep = *this;
    rep->build(*block);

    return rep;
  }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

