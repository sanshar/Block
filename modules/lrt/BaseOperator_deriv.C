/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


// Written by N.N. with reference to Jon's LRT code

#include "BaseOperator.h"
#include "MatrixBLAS.h"
#include "spinblock.h"
#include "distribute.h"
#include "tensor_operator.h"
#include "blas_calls.h"

namespace SpinAdapted{

//Renormalization functions for core and virtual operators                                                                                
void SparseMatrix::renormalise_transform_deriv(
		 const SparseMatrix& op,
		 const std::vector<Matrix>& rotate_matrix,
		 const std::vector<Matrix>& rotate_matrix_deriv,
		 bool is_stateTrans,
		 const StateInfo *stateinfo)
{
  ObjectMatrix<Matrix> tmp       = op.operatorMatrix;
  ObjectMatrix<Matrix> tmp_deriv =    operatorMatrix;

  this->allocate(*stateinfo); // new allocations

  int newQ = 0;
  for (int Q = 0; Q < rotate_matrix.size (); ++Q) {
    if (rotate_matrix[Q].Ncols () != 0) {
      int newQPrime = 0;
      for (int QPrime = 0; QPrime < rotate_matrix.size (); ++QPrime) {
        if (rotate_matrix[QPrime].Ncols () != 0) {
          if (this->allowedQuantaMatrix (newQ, newQPrime)) {
            MatrixRotate (rotate_matrix[Q], tmp_deriv(Q, QPrime), rotate_matrix[QPrime], this->operatorMatrix (newQ, newQPrime) );
            if (!is_stateTrans) {
              MatrixRotate (rotate_matrix[Q], tmp(Q, QPrime), rotate_matrix_deriv[QPrime], this->operatorMatrix (newQ, newQPrime) );
            }
            else {
              MatrixRotate (rotate_matrix_deriv[Q], tmp(Q, QPrime), rotate_matrix[QPrime], this->operatorMatrix (newQ, newQPrime) );
            }
          }
          ++newQPrime;
        }
      }
      ++newQ;
    }
  }
}

void SparseMatrix::build_and_renormalise_transform_deriv(
		 const SparseMatrix& op,
		 SpinBlock *big,
		 const opTypes &ot_deriv,
		 const opTypes &ot,
		 const std::vector<Matrix>& rotate_matrix,
		 const std::vector<Matrix>& rotate_matrix_deriv,
		 bool is_stateTrans,
		 const StateInfo *newStateInfo)
{
  boost::shared_ptr<SparseMatrix> tmp;
  if (orbs.size() == 0)
    tmp = big->get_op_rep(ot, deltaQuantum);
  if (orbs.size() == 1)
    tmp = big->get_op_rep(ot, deltaQuantum, orbs[0]);
  if (orbs.size() == 2)
    tmp = big->get_op_rep(ot, deltaQuantum, orbs[0], orbs[1]);

  boost::shared_ptr<SparseMatrix> tmp_deriv;
  if (orbs.size() == 0)
    tmp_deriv = big->get_op_rep(ot_deriv, deltaQuantum);
  if (orbs.size() == 1)
    tmp_deriv = big->get_op_rep(ot_deriv, deltaQuantum, orbs[0]);
  if (orbs.size() == 2)
    tmp_deriv = big->get_op_rep(ot_deriv, deltaQuantum, orbs[0], orbs[1]);

  tmp_deriv->built = true;

  this->allocate(*newStateInfo);
  this->built = true;

  int newQ = 0;
  for (int Q = 0; Q < rotate_matrix.size (); ++Q) {
    if (rotate_matrix[Q].Ncols () != 0) {
      int newQPrime = 0;
      for (int QPrime = 0; QPrime < rotate_matrix.size (); ++QPrime) {
        if (rotate_matrix[QPrime].Ncols () != 0) {
          if (this->allowedQuantaMatrix (newQ, newQPrime)) {
            MatrixRotate (rotate_matrix[Q], tmp_deriv->operatorMatrix(Q, QPrime), rotate_matrix[QPrime], this->operatorMatrix (newQ, newQPrime) );
            if (!is_stateTrans) {
              MatrixRotate (rotate_matrix[Q], tmp->operatorMatrix(Q, QPrime), rotate_matrix_deriv[QPrime], this->operatorMatrix (newQ, newQPrime) );
            }
            else {
              MatrixRotate (rotate_matrix_deriv[Q], tmp->operatorMatrix(Q, QPrime), rotate_matrix[QPrime], this->operatorMatrix (newQ, newQPrime) );
            }
          }
          ++newQPrime;
        }
      }
      ++newQ;
    }
  }
}

}
