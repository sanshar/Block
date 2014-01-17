/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//
// Extension of Operators.C for additional 2-index operators needed for NPDMs
//FIXME there's a lot of duplication, especially in build_from_disk... Templates??
//
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "op_components.h"
#include "BaseOperator.h"
#include "spinblock.h"
#include "operatorfunctions.h"
#include "tensor_operator.h"
#include "two_index_ops.h"

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Des,Des)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinAdapted::DesDes::build(const SpinBlock& b)
{
//cout << "building DesDes renormalized operator...\n";
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_stateInfo());
  Sign = 1;

  SpinBlock* sysBlock = b.get_leftBlock();
  SpinBlock* dotBlock = b.get_rightBlock();

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
//cout << "indices  " << i << " " << j << std::endl;
//cout << "mpirank = " << mpigetrank() << endl;

  if (sysBlock->get_op_array(DES_DES).has_local_index(i, j)) {      
//cout << "maw: sys(i,j)\n";
    const boost::shared_ptr<SparseMatrix>& op = sysBlock->get_op_rep(DES_DES, deltaQuantum, i,j);
    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  else if (dotBlock->get_op_array(DES_DES).has_local_index(i, j)) {
//cout << "maw: dot(i,j)\n";
    const boost::shared_ptr<SparseMatrix> op = dotBlock->get_op_rep(DES_DES, deltaQuantum, i,j);
    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *op, &b, &(b.get_stateInfo()), *this);
  }  
  else if (sysBlock->get_op_array(CRE).has_local_index(i)) {
//cout << "maw: sys(i)\n";
    Transposeview opD1 = Transposeview( sysBlock->get_op_rep(CRE, getSpinQuantum(i), i) );
    Transposeview opD2 = Transposeview( dotBlock->get_op_rep(CRE, getSpinQuantum(j), j) );
    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, opD1, opD2, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  else if (dotBlock->get_op_array(CRE).has_local_index(i)) {
//cout << "maw: dot(i)\n";
    Transposeview opD1 = Transposeview( dotBlock->get_op_rep(CRE, getSpinQuantum(i), i) );
    Transposeview opD2 = Transposeview( sysBlock->get_op_rep(CRE, getSpinQuantum(j), j) );
    double parity = getCommuteParity( opD1.get_deltaQuantum(), opD2.get_deltaQuantum(), get_deltaQuantum() );
    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, opD1, opD2, &b, &(b.get_stateInfo()), *this, 1.0*parity);
  }
  else
    assert(false);

  dmrginp.makeopsT -> stop();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

double SpinAdapted::DesDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0], 
      J = get_orbs()[1]; //convert spatial id to spin id because slaters need that

  IrrepSpace sym = deltaQuantum.get_symm();
  int irrep = deltaQuantum.get_symm().getirrep();
  int spin = deltaQuantum.get_s();

  TensorOp D1(I, -1), D2(J, -1);
  TensorOp DD = D1.product(D2, spin, irrep);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum, ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, DD, ladder[i]) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::DesDes::getworkingrepresentation(const SpinBlock* block)
{
  assert(this->get_initialised());
  if (this->get_built())
    {
      return boost::shared_ptr<DesDes>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
    }
  else
    {
      boost::shared_ptr<SparseMatrix> rep(new DesDes);
      *rep = *this;
      rep->build(*block);
      return rep;
    }

}

//===========================================================================================================================================================
// 3PDM operators
//===========================================================================================================================================================

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Des,Cre)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinAdapted::DesCre::build(const SpinBlock& b)
{
//cout << "building DesCre renormalized operator...\n";
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_stateInfo());
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
//cout << "indices  " << i << " " << j << std::endl;
//cout << "mpirank = " << mpigetrank() << endl;

  SpinBlock* sysBlock = b.get_leftBlock();
  SpinBlock* dotBlock = b.get_rightBlock();

  if (sysBlock->get_op_array(DES_CRE).has_local_index(i, j))
  {      
    const boost::shared_ptr<SparseMatrix>& op = sysBlock->get_op_rep(DES_CRE, deltaQuantum, i,j);
    SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  else if (dotBlock->get_op_array(DES_CRE).has_local_index(i, j))
  {
    const boost::shared_ptr<SparseMatrix> op = dotBlock->get_op_rep(DES_CRE, deltaQuantum, i,j);
    SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *op, &b, &(b.get_stateInfo()), *this);
  }  
  else if (sysBlock->get_op_array(CRE).has_local_index(i))
  {
    Transposeview opD = Transposeview( sysBlock->get_op_rep(CRE, getSpinQuantum(i), i) );
    const boost::shared_ptr<SparseMatrix> opC = dotBlock->get_op_rep(CRE, getSpinQuantum(j), j);
    SpinAdapted::operatorfunctions::TensorProduct(sysBlock, opD, *opC, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  else if (dotBlock->get_op_array(CRE).has_local_index(i))
  {
    Transposeview opD = Transposeview( dotBlock->get_op_rep(CRE, getSpinQuantum(i), i) );
    const boost::shared_ptr<SparseMatrix> opC = sysBlock->get_op_rep(CRE, getSpinQuantum(j), j);
    double parity = getCommuteParity( opD.get_deltaQuantum(), opC->get_deltaQuantum(), get_deltaQuantum() );
    SpinAdapted::operatorfunctions::TensorProduct(dotBlock, opD, *opC, &b, &(b.get_stateInfo()), *this, 1.0*parity);
  }
  else
    assert(false);
//pout << "opDC\n";
//pout << *this;

  dmrginp.makeopsT -> stop();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//MAW debug
void SpinAdapted::DesCre::build_in_csf_space(const SpinBlock& b) 
{
//pout << "building DesCre in CSF space as a product..\n";
assert(false);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

double SpinAdapted::DesCre::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0], 
    J = get_orbs()[1]; //convert spatial id to spin id because slaters need that

  IrrepSpace sym = deltaQuantum.get_symm();
  int irrep = deltaQuantum.get_symm().getirrep();
  int spin = deltaQuantum.get_s();

  TensorOp D(I, -1), C(J, 1);
  TensorOp DC = D.product(C, spin, irrep);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum, ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, DC, ladder[i]) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::DesCre::getworkingrepresentation(const SpinBlock* block)
{
  assert(this->get_initialised());
  if (this->get_built())
    {
      return boost::shared_ptr<DesCre>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
    }
  else
    {
      boost::shared_ptr<SparseMatrix> rep(new DesCre);
      *rep = *this;
      rep->build(*block);
      return rep;
    }

}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
