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
#include "fourpdm_operators.h"
#include "spinblock.h"
#include "operatorfunctions.h"
#include "tensor_operator.h"


//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Cre,Des,Des)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//
//void SpinAdapted::CreCreDesDes::build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs)
//{
//assert(false);
////cout << "building CreCreDesDes renormalized operator on disk...\n";
////FIXME timer
////  dmrginp.makeopsT -> start();
//  built = true;
////FIXME allocation/deallocation
//  allocate(b.get_stateInfo());
//  Sign = 1;
//
//  const int i = get_orbs()[0];
//  const int j = get_orbs()[1];
//  const int k = get_orbs()[2];
////cout << "i,j,k = " << i << " " << j << " " << k << endl;
//
//  SpinBlock* sysBlock = b.get_leftBlock();
//  SpinBlock* dotBlock = b.get_rightBlock();
//
//  // Sys has i,j,k
//  if (sysBlock->get_op_array(DES_DES_CRE).has_local_index(i,j,k)) {
////cout << "maw disk sys(i,j,k)\n";
//    assert( sysBlock->get_op_array(DES_DES_CRE).get_element(i,j,k).at(0)->get_built_on_disk() );
//    // Retrieve from disk 3-index operator on sys block
//    boost::shared_ptr<SparseMatrix> op (new CreCreDesDes);
//    boost::archive::binary_iarchive load_op(sysfs);
//    load_op >> *op;
//    // Assume this is the operator we want, and expand it to the big block
//    assert( i == op->get_orbs()[0] );
//    assert( j == op->get_orbs()[1] );
//    assert( k == op->get_orbs()[2] );
//    // Build according to previous build_pattern
//    std::string build_pattern_old = op->get_build_pattern();
//    build_pattern = build_pattern_old;
//    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
//    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
//    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *op, &b, &(b.get_stateInfo()), *this);
//  }
//  // Dot has i,j,k
//  else if (dotBlock->get_op_array(DES_DES_CRE).has_local_index(i,j,k)) {
////cout << "maw dot(i,j,k)\n";
//    assert( dotBlock->get_op_array(DES_DES_CRE).get_element(i,j,k).at(0)->get_built_on_disk() );
//    // Retrieve from disk 3-index operator on dot block
//    boost::shared_ptr<SparseMatrix> op (new CreCreDesDes);
//    boost::archive::binary_iarchive load_op(dotfs);
//    load_op >> *op;
//    // Assume this is the operator we want, and expand it to the big block
//    assert( i == op->get_orbs()[0] );
//    assert( j == op->get_orbs()[1] );
//    assert( k == op->get_orbs()[2] );
//    // Build according to previous build_pattern
//    std::string build_pattern_old = op->get_build_pattern();
//    build_pattern = build_pattern_old;
//    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
//    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
//    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *op, &b, &(b.get_stateInfo()), *this);
//  }
//  else {
//    // Build from in-core 1 and 2-index operators
//    build(b);
//  }
//
////  dmrginp.makeopsT -> stop();
//
//}
//
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinAdapted::CreCreDesDes::build(const SpinBlock& b)
{
assert(false);
//cout << "building CreCreDesDes renormalized operator...\n";
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_stateInfo());
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];

  SpinBlock* sysBlock = b.get_leftBlock();
  SpinBlock* dotBlock = b.get_rightBlock();

  // Sys has i,j,k
  if (sysBlock->get_op_array(DES_DES_CRE).has_local_index(i,j,k)) {
pout << "maw sys(i,j,k)\n";
    std::string build_pattern_old = sysBlock->get_op_array(DES_DES_CRE).get_element(i,j,k).at(0)->get_build_pattern();
    assert( build_pattern_old == sysBlock->get_op_array(DES_DES_CRE).get_element(i,j,k).at(1)->get_build_pattern() );
    assert( build_pattern_old == sysBlock->get_op_array(DES_DES_CRE).get_element(i,j,k).at(2)->get_build_pattern() );
    const boost::shared_ptr<SparseMatrix>& op = sysBlock->get_op_rep(DES_DES_CRE, quantum_ladder.at(build_pattern_old), i,j,k);
    // Build according to previous build_pattern
    build_pattern = build_pattern_old;
    assert( get_quantum_ladder().at( build_pattern ).at(1) == op->get_quantum_ladder().at( build_pattern ).at(1) );
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  // Dot has i,j,k
  else if (dotBlock->get_op_array(DES_DES_CRE).has_local_index(i,j,k)) {
pout << "maw dot(i,j,k)\n";
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
    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  // Dot has j,k;
  else if (dotBlock->get_op_array(DES_CRE).has_local_index(j,k)) {
pout << "maw dot(j,k)\n";
    assert( j == k );
    assert( sysBlock->get_op_array(CRE).has_local_index(i) );
    build_pattern = "(D(DC))";
    Transposeview opD = Transposeview( sysBlock->get_op_rep(CRE, getSpinQuantum(i), i) );
    const boost::shared_ptr<SparseMatrix>& opDC = dotBlock->get_op_rep(DES_CRE, quantum_ladder.at(build_pattern).at(0), j,k);
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, opD, *opDC, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  // Sys has i,j
  else if (sysBlock->get_op_array(CRE_CRE).has_local_index(i,j)) {
pout << "maw sys(i,j)\n";
    assert( dotBlock->get_op_array(CRE).has_local_index(k) );
    build_pattern = "((DD)C)";
//FIXME transpose indices / commute ?? minus sign if i!=j ??  c.f. CDD case above?
    Transposeview opDD = Transposeview( sysBlock->get_op_rep(CRE_CRE, quantum_ladder.at(build_pattern).at(0), i,j) );
    const boost::shared_ptr<SparseMatrix>& opC  = dotBlock->get_op_rep(CRE, getSpinQuantum(k), k);
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, opDD, *opC, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  // Sys has j,k;
  else if (sysBlock->get_op_array(DES_CRE).has_local_index(j,k)) {
pout << "maw sys(j,k)\n";
    assert( dotBlock->get_op_array(CRE).has_local_index(i) );
    build_pattern = "(D(DC))";
    Transposeview opD = Transposeview( dotBlock->get_op_rep(CRE, getSpinQuantum(i), i) );
    const boost::shared_ptr<SparseMatrix>& opDC = sysBlock->get_op_rep(DES_CRE, quantum_ladder.at(build_pattern).at(0), j,k);
    // Tensor product of dot*sys so need to take into account parity factors
    set_deltaQuantum() = get_quantum_ladder().at( build_pattern ).at(1);
//FIXME parity broken ??
    double parity = getCommuteParity( opD.get_deltaQuantum(), opDC->get_deltaQuantum(), get_deltaQuantum() );
    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, opD, *opDC, &b, &(b.get_stateInfo()), *this, 1.0*parity);
  }
  // Dot has i,j
  else if (dotBlock->get_op_array(CRE_CRE).has_local_index(i,j)) {
pout << "maw dot(i,j)\n";
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
    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, opDD, *opC, &b, &(b.get_stateInfo()), *this, 1.0*parity);
  }
  else assert(false);

  dmrginp.makeopsT -> stop();

}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

double SpinAdapted::CreCreDesDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
assert(false);
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

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CreCreDesDes::getworkingrepresentation(const SpinBlock* block)
{
assert(false);
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



