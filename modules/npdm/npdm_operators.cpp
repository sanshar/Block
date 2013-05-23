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
//#include "csf.h"
#include "spinblock.h"
//#include "couplingCoeffs.h"
#include "operatorfunctions.h"
//#include "opxop.h"
//#include "operatorloops.h"
//#include "distribute.h"
#include "tensor_operator.h"
//#include "SpinQuantum.h"


//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Cre,Des)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinAdapted::CreCreDes::build(const SpinBlock& b)
{
pout << "building CreCreDes renormalized operator...\n";
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_stateInfo());
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];

  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder();
  assert( quantum_ladder.size() == 2 );
  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);

  // Below, note transposes for (i,j) since 2-index ops are built with i>=j 

pout << "indices  " << i << " " << j << " " << k << std::endl;
  // Forward
  if (leftBlock->get_op_array(CRE_CRE_DES).has(i,j,k)) {
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE_CRE_DES, quantum_ladder, i,j,k);
    SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  else if (rightBlock->get_op_array(CRE_CRE_DES).has(i,j,k)) {
    const boost::shared_ptr<SparseMatrix>& op = rightBlock->get_op_rep(CRE_CRE_DES, quantum_ladder, i,j,k);
    SpinAdapted::operatorfunctions::TensorTrace(rightBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  else if (leftBlock->get_op_array(CRE_CRE).has(j,i)) {
    assert (rightBlock->get_op_array(CRE).has(k));
    const boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(CRE_CRE, deltaQuantum12, j,i);
    Transposeview op2 = Transposeview(rightBlock->get_op_rep(CRE, getSpinQuantum(k), k));
    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, op2, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  // Want to build as ((CC)D), not (C(CD)), but FIXME this is relatively expensive!
//FIXME commute indices???
  else if (b.get_op_array(CRE_CRE).has(j,i)) {
    // Get CC on big block
    const boost::shared_ptr<SparseMatrix> op1 = b.get_op_rep(CRE_CRE, deltaQuantum12, j,i);
    // Get D on big block
    assert (b.get_op_array(CRE).has(k));
    Transposeview op2 = Transposeview( b.get_op_rep(CRE, getSpinQuantum(k), k) );
    // Build tensor product operator directly in 4M*4M space
//    boost::shared_ptr<SparseMatrix> newOp (new Cre);
//    newOp->set_orbs().push_back(i);
//    newOp->set_orbs().push_back(j);
//    newOp->set_orbs().push_back(k);
//    newOp->set_initialised() = true;
//    newOp->set_fermion() = true; //is_fermion;
//    operatorfunctions::Product(&b, *op1, op2, *newOp, 1.0 );
//    operatorfunctions::Product(&b, *op1, op2, *this, 1.0 );
  }
//  else if (rightBlock->get_op_array(CRE_CRE).has(j,i)) {
//pout << "hello6\n";
//    assert (false);
//    assert (leftBlock->get_op_array(CRE).has(k));
//    const boost::shared_ptr<SparseMatrix> op1 = rightBlock->get_op_rep(CRE_CRE, deltaQuantum12, i,j);
//    Transposeview op2 = Transposeview(leftBlock->get_op_rep(CRE, getSpinQuantum(k), k));
//    double parity = getCommuteParity(op1->get_deltaQuantum(), op2.get_deltaQuantum(), get_deltaQuantum());
//    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, op2, &b, &(b.get_stateInfo()), *this, 1.0*parity);
//  }
//  else if (rightBlock->get_op_array(CRE).has(i)) {
//pout << "hello7\n";
//    assert (false);
//  }
  else {
    assert (false);
  }
  dmrginp.makeopsT -> stop();
pout << "done!\n";

}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

double SpinAdapted::CreCreDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];

  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder();
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

pout << "redMatrixElement indices:\n";
pout << I << "  " << J << "  " << K << std::endl;
pout << "spin composition:\n";
pout << spin12/2.0 << "  " << spin123/2.0 << std::endl;
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

  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder();
  assert( quantum_ladder.size() == 2 );
  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);

  // Below, note transposes for (i,j) since 2-index ops are built with i>=j 

pout << "indices  " << i << " " << j << " " << k << std::endl;
pout <<  (leftBlock->get_op_array(CRE_CRE).has(j,i)) << std::endl;
  // Forward
  if (leftBlock->get_op_array(CRE_CRE_DES).has(i,j,k)) {
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE_CRE_DES, quantum_ladder, i,j,k);
    SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  else if (rightBlock->get_op_array(CRE_CRE_DES).has(i,j,k)) {
    const boost::shared_ptr<SparseMatrix>& op = rightBlock->get_op_rep(CRE_CRE_DES, quantum_ladder, i,j,k);
    SpinAdapted::operatorfunctions::TensorTrace(rightBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  else if (leftBlock->get_op_array(CRE_CRE).has(j,i)) {
    assert (rightBlock->get_op_array(CRE).has(k));
    const boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(CRE_CRE, deltaQuantum12, j,i);
    Transposeview op2 = Transposeview(rightBlock->get_op_rep(CRE, getSpinQuantum(k), k));
    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, op2, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  // Want to build as ((CC)D), not (C(CD)), but FIXME this is relatively expensive!
  else if (b.get_op_array(CRE_CRE).has(j,i)) {
    // Get CC on big block
    const boost::shared_ptr<SparseMatrix> op1 = b.get_op_rep(CRE_CRE, deltaQuantum12, j,i);
    // Get D on big block
    assert (b.get_op_array(CRE).has(k));
    Transposeview op2 = Transposeview( b.get_op_rep(CRE, getSpinQuantum(k), k) );
    // Build tensor product operator directly in 4M*4M space
    boost::shared_ptr<SparseMatrix> newOp (new Cre);
    newOp->set_orbs().push_back(i);
    newOp->set_orbs().push_back(j);
    newOp->set_orbs().push_back(k);
    newOp->set_initialised() = true;
    newOp->set_fermion() = true; //is_fermion;
    operatorfunctions::Product(&b, *op1, op2, *newOp, 1.0 );
  }
//  else if (rightBlock->get_op_array(CRE_CRE).has(j,i)) {
//pout << "hello6\n";
//    assert (false);
//    assert (leftBlock->get_op_array(CRE).has(k));
//    const boost::shared_ptr<SparseMatrix> op1 = rightBlock->get_op_rep(CRE_CRE, deltaQuantum12, i,j);
//    Transposeview op2 = Transposeview(leftBlock->get_op_rep(CRE, getSpinQuantum(k), k));
//    double parity = getCommuteParity(op1->get_deltaQuantum(), op2.get_deltaQuantum(), get_deltaQuantum());
//    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, op2, &b, &(b.get_stateInfo()), *this, 1.0*parity);
//  }
//  else if (rightBlock->get_op_array(CRE).has(i)) {
//pout << "hello7\n";
//    assert (false);
//  }
  else {
    assert (false);
  }
  dmrginp.makeopsT -> stop();
pout << "done!\n";

}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

double SpinAdapted::CreDesDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];

  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder();
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

pout << "redMatrixElement indices:\n";
pout << I << "  " << J << "  " << K << std::endl;
pout << "spin composition:\n";
pout << spin12/2.0 << "  " << spin123/2.0 << std::endl;
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

