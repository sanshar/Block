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
 void SpinAdapted::CreCreDesDes::build(const SpinBlock& b) { 
      dmrginp.makeopsT -> start();
      built = true;
      allocate(b.get_braStateInfo(), b.get_ketStateInfo());
    
      const int i = get_orbs()[0];
      const int j = get_orbs()[1];
      const int k = get_orbs()[2];
      const int l = get_orbs()[3];
      SpinBlock* leftBlock = b.get_leftBlock();
      SpinBlock* rightBlock = b.get_rightBlock();
    
      if (leftBlock->get_op_array(CRE_CRE_DES_DES).has(i,j,k,l))
      {      
        const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE_CRE_DES_DES, quantum_ladder, i,j,k,l);
        if (rightBlock->get_sites().size() == 0) 
          SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
        dmrginp.makeopsT -> stop();
        return;
      }
      assert(false && "Only build CRECREDESDES in the starting block when spin-embeding is used");
    }

double SpinAdapted::CreCreDesDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  assert( build_pattern == "(((CC)(D))(D))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
  int L = get_orbs()[3];

  // Must take into account how the 4-index is built from a combination of the 2-index ops
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("(((CC)(D))(D))");
  assert( quantum_ladder.size() == 3 );

  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  SpinQuantum deltaQuantum1234 = quantum_ladder.at(2);
//FIXME components != 0
  deltaQuantum[0] = deltaQuantum1234;

  // Spin quantum data for CC
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s().getirrep();
  // Spin quantum data for (CC)D
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123= deltaQuantum123.get_s().getirrep();
  // Spin quantum data for total operator
  IrrepSpace sym1234 = deltaQuantum1234.get_symm();
  int irrep1234 = deltaQuantum1234.get_symm().getirrep();
  int spin1234 = deltaQuantum1234.get_s().getirrep();

  TensorOp C1(I, 1); 
  TensorOp C2(J, 1); 
  TensorOp D3(K,-1); 
  TensorOp D4(L,-1); 

  TensorOp CC = C1.product(C2, spin12, irrep12);
  TensorOp CCD = CC.product(D3, spin123, irrep123);
  TensorOp CCDD = CCD.product(D4, spin1234, irrep1234);

//FIXME loop over deltaQuantum components
  int j = 0;
  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
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
//  (Cre,Des,Cre,Des)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
void SpinAdapted::CreDesCreDes::build(const SpinBlock& b) { 
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
  const int l = get_orbs()[3];
  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_DES_CRE_DES).has(i,j,k,l))
  {      
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE_DES_CRE_DES, quantum_ladder, i,j,k,l);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    dmrginp.makeopsT -> stop();
    return;
  }
  assert(false && "Only build CREDESCREDES in the starting block when spin-embeding is used");
}

double SpinAdapted::CreDesCreDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  assert( build_pattern == "(((CD)(C))(D))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
  int L = get_orbs()[3];

  // Must take into account how the 4-index is built from a combination of the 2-index ops
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("(((CD)(C))(D))");
  assert( quantum_ladder.size() == 3 );

  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  SpinQuantum deltaQuantum1234 = quantum_ladder.at(2);
  deltaQuantum[0] = deltaQuantum1234;

  // Spin quantum data for CD
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s().getirrep();
  // Spin quantum data for (CD)C
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123= deltaQuantum123.get_s().getirrep();
  // Spin quantum data for total operator
  IrrepSpace sym1234 = deltaQuantum1234.get_symm();
  int irrep1234 = deltaQuantum1234.get_symm().getirrep();
  int spin1234 = deltaQuantum1234.get_s().getirrep();

  TensorOp C1(I, 1); 
  TensorOp D2(J,-1); 
  TensorOp C3(K, 1); 
  TensorOp D4(L,-1); 

  TensorOp CD = C1.product(D2, spin12, irrep12);
  TensorOp CDC = CD.product(C3, spin123, irrep123);
  TensorOp CDCD = CDC.product(D4, spin1234, irrep1234);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[0], ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CDCD, ladder[i]) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CreDesCreDes::getworkingrepresentation(const SpinBlock* block)
{
  assert(this->get_initialised());
  if (this->get_built()) {
    return boost::shared_ptr<CreDesCreDes>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
  }
  else {
    boost::shared_ptr<SparseMatrix> rep(new CreDesCreDes);
    *rep = *this;
    rep->build(*block);

    return rep;
  }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Des,Des,Cre)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
void SpinAdapted::CreDesDesCre::build(const SpinBlock& b) { 
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
  const int l = get_orbs()[3];
  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_DES_DES_CRE).has(i,j,k,l))
  {      
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE_DES_DES_CRE, quantum_ladder, i,j,k,l);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    dmrginp.makeopsT -> stop();
    return;
  }
  assert(false && "Only build CREDESDESCRE in the starting block when spin-embeding is used");
}

double SpinAdapted::CreDesDesCre::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  assert( build_pattern == "(((CD)(D))(C))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
  int L = get_orbs()[3];

  // Must take into account how the 4-index is built from a combination of the 2-index ops
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("(((CD)(D))(C))");
  assert( quantum_ladder.size() == 3 );

  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  SpinQuantum deltaQuantum1234 = quantum_ladder.at(2);
  deltaQuantum[0] = deltaQuantum1234;

  // Spin quantum data for CD
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s().getirrep();
  // Spin quantum data for (CD)D
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123= deltaQuantum123.get_s().getirrep();
  // Spin quantum data for total operator
  IrrepSpace sym1234 = deltaQuantum1234.get_symm();
  int irrep1234 = deltaQuantum1234.get_symm().getirrep();
  int spin1234 = deltaQuantum1234.get_s().getirrep();

  TensorOp C1(I, 1); 
  TensorOp D2(J,-1); 
  TensorOp D3(K,-1); 
  TensorOp C4(L, 1); 

  TensorOp CD = C1.product(D2, spin12, irrep12);
  TensorOp CDD = CD.product(D3, spin123, irrep123);
  TensorOp CDDC = CDD.product(C4, spin1234, irrep1234);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[0], ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CDDC, ladder[i]) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CreDesDesCre::getworkingrepresentation(const SpinBlock* block)
{
  assert(this->get_initialised());
  if (this->get_built()) {
    return boost::shared_ptr<CreDesDesCre>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
  }
  else {
    boost::shared_ptr<SparseMatrix> rep(new CreDesDesCre);
    *rep = *this;
    rep->build(*block);

    return rep;
  }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Des,Des,Des)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
void SpinAdapted::CreDesDesDes::build(const SpinBlock& b) { 
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
  const int l = get_orbs()[3];
  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_DES_DES_DES).has(i,j,k,l))
  {      
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE_DES_DES_DES, quantum_ladder, i,j,k,l);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    dmrginp.makeopsT -> stop();
    return;
  }
  assert(false && "Only build CREDESDESDES in the starting block when spin-embeding is used");
}

double SpinAdapted::CreDesDesDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  assert( build_pattern == "(((CD)(D))(D))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
  int L = get_orbs()[3];

  // Must take into account how the 4-index is built from a combination of the 2-index ops
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("(((CD)(D))(D))");
  assert( quantum_ladder.size() == 3 );

  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  SpinQuantum deltaQuantum1234 = quantum_ladder.at(2);
  deltaQuantum[0] = deltaQuantum1234;

  // Spin quantum data for CD
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s().getirrep();
  // Spin quantum data for (CD)D
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123= deltaQuantum123.get_s().getirrep();
  // Spin quantum data for total operator
  IrrepSpace sym1234 = deltaQuantum1234.get_symm();
  int irrep1234 = deltaQuantum1234.get_symm().getirrep();
  int spin1234 = deltaQuantum1234.get_s().getirrep();

  TensorOp C1(I, 1); 
  TensorOp D2(J,-1); 
  TensorOp D3(K,-1); 
  TensorOp D4(L,-1); 

  TensorOp CD = C1.product(D2, spin12, irrep12);
  TensorOp CDD = CD.product(D3, spin123, irrep123);
  TensorOp CDDD = CDD.product(D4, spin1234, irrep1234);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[0], ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CDDD, ladder[i]) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CreDesDesDes::getworkingrepresentation(const SpinBlock* block)
{
  assert(this->get_initialised());
  if (this->get_built()) {
    return boost::shared_ptr<CreDesDesDes>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
  }
  else {
    boost::shared_ptr<SparseMatrix> rep(new CreDesDesDes);
    *rep = *this;
    rep->build(*block);

    return rep;
  }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Cre,Cre,Des)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
void SpinAdapted::CreCreCreDes::build(const SpinBlock& b) { 
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
  const int l = get_orbs()[3];
  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_CRE_CRE_DES).has(i,j,k,l))
  {      
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE_CRE_CRE_DES, quantum_ladder, i,j,k,l);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    dmrginp.makeopsT -> stop();
    return;
  }
  assert(false && "Only build CRECRECREDES in the starting block when spin-embeding is used");
}

double SpinAdapted::CreCreCreDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  assert( build_pattern == "(((CC)(C))(D))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
  int L = get_orbs()[3];

  // Must take into account how the 4-index is built from a combination of the 2-index ops
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("(((CC)(C))(D))");
  assert( quantum_ladder.size() == 3 );

  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  SpinQuantum deltaQuantum1234 = quantum_ladder.at(2);
  deltaQuantum[0] = deltaQuantum1234;

  // Spin quantum data for CC
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s().getirrep();
  // Spin quantum data for (CC)C
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123= deltaQuantum123.get_s().getirrep();
  // Spin quantum data for total operator
  IrrepSpace sym1234 = deltaQuantum1234.get_symm();
  int irrep1234 = deltaQuantum1234.get_symm().getirrep();
  int spin1234 = deltaQuantum1234.get_s().getirrep();

  TensorOp C1(I, 1); 
  TensorOp C2(J, 1); 
  TensorOp C3(K, 1); 
  TensorOp D4(L,-1); 

  TensorOp CC = C1.product(C2, spin12, irrep12);
  TensorOp CCC = CC.product(C3, spin123, irrep123);
  TensorOp CCCD = CCC.product(D4, spin1234, irrep1234);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[0], ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CCCD, ladder[i]) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CreCreCreDes::getworkingrepresentation(const SpinBlock* block)
{
  assert(this->get_initialised());
  if (this->get_built()) {
    return boost::shared_ptr<CreCreCreDes>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
  }
  else {
    boost::shared_ptr<SparseMatrix> rep(new CreCreCreDes);
    *rep = *this;
    rep->build(*block);

    return rep;
  }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Cre,Des,Cre)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
void SpinAdapted::CreCreDesCre::build(const SpinBlock& b) { 
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
  const int l = get_orbs()[3];
  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_CRE_DES_CRE).has(i,j,k,l))
  {      
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE_CRE_DES_CRE, quantum_ladder, i,j,k,l);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    dmrginp.makeopsT -> stop();
    return;
  }
  assert(false && "Only build CRECREDESCRE in the starting block when spin-embeding is used");
}

double SpinAdapted::CreCreDesCre::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  assert( build_pattern == "(((CC)(D))(C))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
  int L = get_orbs()[3];

  // Must take into account how the 4-index is built from a combination of the 2-index ops
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("(((CC)(D))(C))");
  assert( quantum_ladder.size() == 3 );

  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  SpinQuantum deltaQuantum1234 = quantum_ladder.at(2);
  deltaQuantum[0] = deltaQuantum1234;

  // Spin quantum data for CC
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s().getirrep();
  // Spin quantum data for (CC)D
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123= deltaQuantum123.get_s().getirrep();
  // Spin quantum data for total operator
  IrrepSpace sym1234 = deltaQuantum1234.get_symm();
  int irrep1234 = deltaQuantum1234.get_symm().getirrep();
  int spin1234 = deltaQuantum1234.get_s().getirrep();

  TensorOp C1(I, 1); 
  TensorOp C2(J, 1); 
  TensorOp D3(K,-1); 
  TensorOp C4(L, 1); 

  TensorOp CC = C1.product(C2, spin12, irrep12);
  TensorOp CCD = CC.product(D3, spin123, irrep123);
  TensorOp CCDC = CCD.product(C4, spin1234, irrep1234);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[0], ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CCDC, ladder[i]) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CreCreDesCre::getworkingrepresentation(const SpinBlock* block)
{
  assert(this->get_initialised());
  if (this->get_built()) {
    return boost::shared_ptr<CreCreDesCre>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
  }
  else {
    boost::shared_ptr<SparseMatrix> rep(new CreCreDesCre);
    *rep = *this;
    rep->build(*block);

    return rep;
  }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Des,Cre,Cre)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
void SpinAdapted::CreDesCreCre::build(const SpinBlock& b) { 
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
  const int l = get_orbs()[3];
  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_DES_CRE_CRE).has(i,j,k,l))
  {      
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE_DES_CRE_CRE, quantum_ladder, i,j,k,l);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    dmrginp.makeopsT -> stop();
    return;
  }
  assert(false && "Only build CREDESCRECRE in the starting block when spin-embeding is used");
}

double SpinAdapted::CreDesCreCre::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  assert( build_pattern == "(((CD)(C))(C))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
  int L = get_orbs()[3];

  // Must take into account how the 4-index is built from a combination of the 2-index ops
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("(((CD)(C))(C))");
  assert( quantum_ladder.size() == 3 );

  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  SpinQuantum deltaQuantum1234 = quantum_ladder.at(2);
  deltaQuantum[0] = deltaQuantum1234;

  // Spin quantum data for CD
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s().getirrep();
  // Spin quantum data for (CD)C
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123= deltaQuantum123.get_s().getirrep();
  // Spin quantum data for total operator
  IrrepSpace sym1234 = deltaQuantum1234.get_symm();
  int irrep1234 = deltaQuantum1234.get_symm().getirrep();
  int spin1234 = deltaQuantum1234.get_s().getirrep();

  TensorOp C1(I, 1); 
  TensorOp D2(J,-1); 
  TensorOp C3(K, 1); 
  TensorOp C4(L, 1); 

  TensorOp CD = C1.product(D2, spin12, irrep12);
  TensorOp CDC = CD.product(C3, spin123, irrep123);
  TensorOp CDCC = CDC.product(C4, spin1234, irrep1234);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[0], ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CDCC, ladder[i]) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CreDesCreCre::getworkingrepresentation(const SpinBlock* block)
{
  assert(this->get_initialised());
  if (this->get_built()) {
    return boost::shared_ptr<CreDesCreCre>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
  }
  else {
    boost::shared_ptr<SparseMatrix> rep(new CreDesCreCre);
    *rep = *this;
    rep->build(*block);

    return rep;
  }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Cre,Cre,Cre)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
void SpinAdapted::CreCreCreCre::build(const SpinBlock& b) { 
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
  const int l = get_orbs()[3];
  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_CRE_CRE_CRE).has(i,j,k,l))
  {      
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE_CRE_CRE_CRE, quantum_ladder, i,j,k,l);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    dmrginp.makeopsT -> stop();
    return;
  }
  assert(false && "Only build CRECRECRECRE in the starting block when spin-embeding is used");
}

double SpinAdapted::CreCreCreCre::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  assert( build_pattern == "(((CC)(C))(C))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
  int L = get_orbs()[3];

  // Must take into account how the 4-index is built from a combination of the 2-index ops
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("(((CC)(C))(C))");
  assert( quantum_ladder.size() == 3 );

  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  SpinQuantum deltaQuantum1234 = quantum_ladder.at(2);
  deltaQuantum[0] = deltaQuantum1234;

  // Spin quantum data for CC
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s().getirrep();
  // Spin quantum data for (CC)C
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123= deltaQuantum123.get_s().getirrep();
  // Spin quantum data for total operator
  IrrepSpace sym1234 = deltaQuantum1234.get_symm();
  int irrep1234 = deltaQuantum1234.get_symm().getirrep();
  int spin1234 = deltaQuantum1234.get_s().getirrep();

  TensorOp C1(I, 1); 
  TensorOp C2(J, 1); 
  TensorOp C3(K, 1); 
  TensorOp C4(L, 1); 

  TensorOp CC = C1.product(C2, spin12, irrep12);
  TensorOp CCC = CC.product(C3, spin123, irrep123);
  TensorOp CCCC = CCC.product(C4, spin1234, irrep1234);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[0], ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CCCC, ladder[i]) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CreCreCreCre::getworkingrepresentation(const SpinBlock* block)
{
  assert(this->get_initialised());
  if (this->get_built()) {
    return boost::shared_ptr<CreCreCreCre>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
  }
  else {
    boost::shared_ptr<SparseMatrix> rep(new CreCreCreCre);
    *rep = *this;
    rep->build(*block);

    return rep;
  }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
