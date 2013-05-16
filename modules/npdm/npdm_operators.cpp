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
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_stateInfo());
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];

  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_DES).has(i, j))
  {      
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE_DES, deltaQuantum, i,j);
    SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    dmrginp.makeopsT -> stop();
    return;
  }
  if (rightBlock->get_op_array(CRE_DES).has(i, j))
  {
    const boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_rep(CRE_DES, deltaQuantum, i,j);
    SpinAdapted::operatorfunctions::TensorTrace(rightBlock, *op, &b, &(b.get_stateInfo()), *this);
    dmrginp.makeopsT -> stop();
    return;
  }  
  if (leftBlock->get_op_array(CRE).has(i))
  {
    const boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(i), i);
    Transposeview top2 = Transposeview(rightBlock->get_op_rep(CRE, getSpinQuantum(j), j));
    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, top2, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  else if (rightBlock->get_op_array(CRE).has(i))
  {
    const boost::shared_ptr<SparseMatrix> op1 = rightBlock->get_op_rep(CRE, getSpinQuantum(i), i);
    Transposeview top2 = Transposeview(leftBlock->get_op_rep(CRE, getSpinQuantum(j), j));
    double parity = getCommuteParity(op1->get_deltaQuantum(), top2.get_deltaQuantum(), get_deltaQuantum());
    SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, top2, &b, &(b.get_stateInfo()), *this, 1.0*parity);
  }
  else
    abort();  
  dmrginp.makeopsT -> stop();
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

double SpinAdapted::CreCreDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0], 
    J = get_orbs()[1]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum.get_symm();
  int irrep = deltaQuantum.get_symm().getirrep();
  int spin = deltaQuantum.get_s();

  TensorOp C(I, 1), D(J, -1);
  TensorOp CD = C.product(D, spin, irrep);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum, ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CD, ladder[i]) ;
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
  //assert(this->get_initialised());
  if (this->get_built())
    {
      return boost::shared_ptr<CreCreDes>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
    }
  else
    {
      boost::shared_ptr<SparseMatrix> rep(new CreDes);
      *rep = *this;
      rep->build(*block);
      return rep;
    }

}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
