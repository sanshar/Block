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
void SparseMatrix::renormalise_transform_lrt
(const SparseMatrix& op_0,
 const std::vector<Matrix>& rotate_matrix_0,
 const std::vector<Matrix>& rotate_matrix_i, bool doTrans,
 const StateInfo *stateinfo)
{
  ObjectMatrix<Matrix> tmp_0 = op_0.operatorMatrix;
  ObjectMatrix<Matrix> tmp_i =      operatorMatrix;

  this->allocate(*stateinfo); // new allocations

  int newQ = 0;
  for (int Q = 0; Q < rotate_matrix_0.size (); ++Q) {
    if (rotate_matrix_0[Q].Ncols () != 0) {
      int newQPrime = 0;
      for (int QPrime = 0; QPrime < rotate_matrix_0.size (); ++QPrime) {
        if (rotate_matrix_0[QPrime].Ncols () != 0) {
          if (this->allowedQuantaMatrix (newQ, newQPrime)) {

              MatrixRotate (rotate_matrix_0[Q], tmp_i(Q, QPrime), rotate_matrix_0[QPrime], this->operatorMatrix (newQ, newQPrime) );

            if (!doTrans) {
              MatrixRotate (rotate_matrix_0[Q], tmp_0(Q, QPrime), rotate_matrix_i[QPrime], this->operatorMatrix (newQ, newQPrime) );
            }
            else {
              MatrixRotate (rotate_matrix_i[Q], tmp_0(Q, QPrime), rotate_matrix_0[QPrime], this->operatorMatrix (newQ, newQPrime) );
            }
          }
          ++newQPrime;
        }
      }
      ++newQ;
    }
  }
}

void SparseMatrix::build_and_renormalise_transform_lrt
(SpinBlock *big, const opTypes &ot_0, const opTypes &ot_i,
 const std::vector<Matrix>& rotate_matrix_0,
 const std::vector<Matrix>& rotate_matrix_i, bool doTrans,
 const StateInfo *newStateInfo)
{
  boost::shared_ptr<SparseMatrix> tmp_0;
  if (orbs.size() == 0)
    tmp_0 = big->get_op_rep(ot_0, deltaQuantum);
  if (orbs.size() == 1)
    tmp_0 = big->get_op_rep(ot_0, deltaQuantum, orbs[0]);
  if (orbs.size() == 2)
    tmp_0 = big->get_op_rep(ot_0, deltaQuantum, orbs[0], orbs[1]);

  boost::shared_ptr<SparseMatrix> tmp_i;
  if (orbs.size() == 0)
    tmp_i = big->get_op_rep(ot_i, deltaQuantum);
  if (orbs.size() == 1)
    tmp_i = big->get_op_rep(ot_i, deltaQuantum, orbs[0]);
  if (orbs.size() == 2)
    tmp_i = big->get_op_rep(ot_i, deltaQuantum, orbs[0], orbs[1]);

  tmp_i->built = true;

  this->allocate(*newStateInfo);
  this->built = true;


  int newQ = 0;
  for (int Q = 0; Q < rotate_matrix_0.size (); ++Q) {
    if (rotate_matrix_0[Q].Ncols () != 0) {
      int newQPrime = 0;
      for (int QPrime = 0; QPrime < rotate_matrix_0.size (); ++QPrime) {
        if (rotate_matrix_0[QPrime].Ncols () != 0) {
          if (this->allowedQuantaMatrix (newQ, newQPrime)) {

              MatrixRotate (rotate_matrix_0[Q], tmp_i->operatorMatrix(Q, QPrime), rotate_matrix_0[QPrime],
                            this->operatorMatrix (newQ, newQPrime) );

            if (!doTrans) {
              MatrixRotate (rotate_matrix_0[Q], tmp_0->operatorMatrix(Q, QPrime), rotate_matrix_i[QPrime],
                            this->operatorMatrix (newQ, newQPrime) );
            }
            else {
              MatrixRotate (rotate_matrix_i[Q], tmp_0->operatorMatrix(Q, QPrime), rotate_matrix_0[QPrime],
                            this->operatorMatrix (newQ, newQPrime) );
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
