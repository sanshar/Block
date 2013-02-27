/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "spinblock.h"
#include "op_components.h"
#include "operatorfunctions.h"
#include "opxop.h"
#include "wavefunction.h"
#include <boost/format.hpp>
#include "distribute.h"

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
#include "pario.h"

namespace SpinAdapted {
using namespace operatorfunctions;

// for DMRG-LRT: ( [ L(0i) x R(00) ] + [ L(i0)^(T) x R(00)^(T) ] ) * c
void SpinBlock::multiplyH_lrt_left(Wavefunction& c, Wavefunction* v, int iState, int num_threads) const
{

  Wavefunction *v_array=0, *v_distributed=0, *v_add=0;
  
  opTypes state_index_0i = make_state_index(0, iState);
  opTypes state_index_i0 = make_state_index(iState, 0);

  int maxt = 1;
  initiateMultiThread(v, v_array, v_distributed, MAX_THRD);
  dmrginp.oneelecT -> start();
  dmrginp.s0time -> start();
  boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_array(HAM+state_index_0i).get_local_element(0)[0];
  TensorMultiply(leftBlock, *op, this, c, *v, dmrginp.effective_molecule_quantum() ,1.0, MAX_THRD);

  op = rightBlock->get_op_array(HAM+state_index_0i).get_local_element(0)[0];
  TensorMultiply(rightBlock, *op, this, c, *v, dmrginp.effective_molecule_quantum(), 1.0, MAX_THRD);  

  dmrginp.s0time -> stop();
#ifndef SERIAL
  boost::mpi::communicator world;
  int size = world.size();
#endif

  dmrginp.s1time -> start();
  v_add =  leftBlock->get_op_array(CRE_CRE_DESCOMP+state_index_0i).is_local() ? v_array : v_distributed;
  Functor f = boost::bind(&opxop::generic::cxcddcomp, leftBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), state_index_0i); 
  for_all_multithread(rightBlock->get_op_array(CRE), rightBlock->get_op_array(CRE), f);

  v_add =  rightBlock->get_op_array(CRE_CRE_DESCOMP).is_local() ? v_array : v_distributed;
  f = boost::bind(&opxop::generic::cxcddcomp, rightBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), 0); 
  for_all_multithread(leftBlock->get_op_array(CRE+state_index_0i), leftBlock->get_op_array(CRE+state_index_i0), f);  

  dmrginp.s1time -> stop();

  dmrginp.oneelecT -> stop();

  dmrginp.twoelecT -> start();

  if (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY) {
    if(leftBlock->is_loopblock()) {
      dmrginp.s0time -> start();
      v_add =  rightBlock->get_op_array(CRE_DESCOMP).is_local() ? v_array : v_distributed;
      f = boost::bind(&opxop::generic::cdxcdcomp, rightBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), 0);
      for_all_multithread(leftBlock->get_op_array(CRE_DES+state_index_0i), leftBlock->get_op_array(CRE_DES+state_index_i0), f);

      v_add =  rightBlock->get_op_array(DES_DESCOMP).is_local() ? v_array : v_distributed;
      f = boost::bind(&opxop::generic::ddxcccomp, rightBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), 0);
      for_all_multithread(leftBlock->get_op_array(CRE_CRE+state_index_0i), leftBlock->get_op_array(CRE_CRE+state_index_i0), f);
      dmrginp.s0time -> stop();
    }
    else {
      dmrginp.s0time -> start();
      v_add =  leftBlock->get_op_array(CRE_DESCOMP+state_index_0i).is_local() ? v_array : v_distributed;
      f = boost::bind(&opxop::generic::cdxcdcomp, leftBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), state_index_0i);
      for_all_multithread(rightBlock->get_op_array(CRE_DES), rightBlock->get_op_array(CRE_DES), f);

      v_add =  leftBlock->get_op_array(DES_DESCOMP+state_index_0i).is_local() ? v_array : v_distributed;
      f = boost::bind(&opxop::generic::ddxcccomp, leftBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), state_index_0i);
      for_all_multithread(rightBlock->get_op_array(CRE_CRE), rightBlock->get_op_array(CRE_CRE), f);
      dmrginp.s0time -> stop();
    }
  }
  dmrginp.twoelecT -> stop();

  accumulateMultiThread(v, v_array, v_distributed, MAX_THRD);

}

void SpinBlock::multiplyH_lrt_total(Wavefunction& c, Wavefunction* v, int iState, int num_threads) const
{

  SpinBlock* loopBlock=(leftBlock->is_loopblock()) ? leftBlock : rightBlock;
  SpinBlock* otherBlock = loopBlock == leftBlock ? rightBlock : leftBlock;

  Wavefunction *v_array=0, *v_distributed=0, *v_add=0;
  
  opTypes state_index_0i = make_state_index(0, iState);
  opTypes state_index_i0 = make_state_index(iState, 0);

  int maxt = 1;
  initiateMultiThread(v, v_array, v_distributed, MAX_THRD);
  dmrginp.oneelecT -> start();
  dmrginp.s0time -> start();
  boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_array(HAM+state_index_0i).get_local_element(0)[0];
  TensorMultiply(leftBlock, *op, this, c, *v, dmrginp.effective_molecule_quantum() ,1.0, MAX_THRD);

  op = rightBlock->get_op_array(HAM+state_index_0i).get_local_element(0)[0];
  TensorMultiply(rightBlock, *op, this, c, *v, dmrginp.effective_molecule_quantum(), 1.0, MAX_THRD);  

  dmrginp.s0time -> stop();
#ifndef SERIAL
  boost::mpi::communicator world;
  int size = world.size();
#endif

  dmrginp.s1time -> start();
  v_add =  leftBlock->get_op_array(CRE_CRE_DESCOMP+state_index_0i).is_local() ? v_array : v_distributed;
  Functor f = boost::bind(&opxop::generic::cxcddcomp, leftBlock, _1, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), state_index_0i); 
  for_all_multithread(rightBlock->get_op_array(CRE), rightBlock->get_op_array(CRE), f);

  v_add =  rightBlock->get_op_array(CRE_CRE_DESCOMP).is_local() ? v_array : v_distributed;
  f = boost::bind(&opxop::generic::cxcddcomp, rightBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), 0); 
  for_all_multithread(leftBlock->get_op_array(CRE+state_index_0i), leftBlock->get_op_array(CRE+state_index_i0), f);  

  if(iState > 0) {

  v_add =  leftBlock->get_op_array(CRE_CRE_DESCOMP).is_local() ? v_array : v_distributed;
  f = boost::bind(&opxop::generic::cxcddcomp, leftBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), 0); 
  for_all_multithread(rightBlock->get_op_array(CRE+state_index_0i), rightBlock->get_op_array(CRE+state_index_i0), f);

  v_add =  rightBlock->get_op_array(CRE_CRE_DESCOMP+state_index_0i).is_local() ? v_array : v_distributed;
  f = boost::bind(&opxop::generic::cxcddcomp, rightBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), state_index_0i); 
  for_all_multithread(leftBlock->get_op_array(CRE), leftBlock->get_op_array(CRE), f);  

  }

  dmrginp.s1time -> stop();

  dmrginp.oneelecT -> stop();

  dmrginp.twoelecT -> start();

  if (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY) {
    if(leftBlock->is_loopblock()) {
      dmrginp.s0time -> start();
      v_add =  rightBlock->get_op_array(CRE_DESCOMP).is_local() ? v_array : v_distributed;
      f = boost::bind(&opxop::generic::cdxcdcomp, rightBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), 0);
      for_all_multithread(leftBlock->get_op_array(CRE_DES+state_index_0i), leftBlock->get_op_array(CRE_DES+state_index_i0), f);

      v_add =  rightBlock->get_op_array(DES_DESCOMP).is_local() ? v_array : v_distributed;
      f = boost::bind(&opxop::generic::ddxcccomp, rightBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), 0);
      for_all_multithread(leftBlock->get_op_array(CRE_CRE+state_index_0i), leftBlock->get_op_array(CRE_CRE+state_index_i0), f);

      if(iState > 0) {

      v_add =  leftBlock->get_op_array(CRE_DESCOMP+state_index_0i).is_local() ? v_array : v_distributed;
      f = boost::bind(&opxop::generic::cdxcdcomp, leftBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), state_index_0i);
      for_all_multithread(rightBlock->get_op_array(CRE_DES), rightBlock->get_op_array(CRE_DES), f);

      v_add =  leftBlock->get_op_array(DES_DESCOMP+state_index_0i).is_local() ? v_array : v_distributed;
      f = boost::bind(&opxop::generic::ddxcccomp, leftBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), state_index_0i);
      for_all_multithread(rightBlock->get_op_array(CRE_CRE), rightBlock->get_op_array(CRE_CRE), f);

      }

      dmrginp.s0time -> stop();
    }
    else {
      dmrginp.s0time -> start();
      v_add =  leftBlock->get_op_array(CRE_DESCOMP+state_index_0i).is_local() ? v_array : v_distributed;
      f = boost::bind(&opxop::generic::cdxcdcomp, leftBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), state_index_0i);
      for_all_multithread(rightBlock->get_op_array(CRE_DES), rightBlock->get_op_array(CRE_DES), f);

      v_add =  leftBlock->get_op_array(DES_DESCOMP+state_index_0i).is_local() ? v_array : v_distributed;
      f = boost::bind(&opxop::generic::ddxcccomp, leftBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), state_index_0i);
      for_all_multithread(rightBlock->get_op_array(CRE_CRE), rightBlock->get_op_array(CRE_CRE), f);

      if(iState > 0) {

      v_add =  rightBlock->get_op_array(CRE_DESCOMP).is_local() ? v_array : v_distributed;
      f = boost::bind(&opxop::generic::cdxcdcomp, rightBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), 0);
      for_all_multithread(leftBlock->get_op_array(CRE_DES+state_index_0i), leftBlock->get_op_array(CRE_DES+state_index_i0), f);

      v_add =  rightBlock->get_op_array(DES_DESCOMP).is_local() ? v_array : v_distributed;
      f = boost::bind(&opxop::generic::ddxcccomp, rightBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), 0);
      for_all_multithread(leftBlock->get_op_array(CRE_CRE+state_index_0i), leftBlock->get_op_array(CRE_CRE+state_index_i0), f);

      }
      dmrginp.s0time -> stop();
    }
  }
  dmrginp.twoelecT -> stop();

  accumulateMultiThread(v, v_array, v_distributed, MAX_THRD);

}

void SpinBlock::rotatebyRitzVectors(const Matrix& alpha, int nroot)
{
  if(has(HAM))
    RotateOps::rotate_operator_arrays(alpha, *this, HAM, nroots);
  if(has(CRE))
    RotateOps::rotate_operator_arrays(alpha, *this, CRE, nroots);
  if(has(CRE_DES))
    RotateOps::rotate_operator_arrays(alpha, *this, CRE_DES, nroots);
  if(has(CRE_CRE))
    RotateOps::rotate_operator_arrays(alpha, *this, CRE_CRE, nroots);
  if(has(CRE_DESCOMP))
    RotateOps::rotate_operator_arrays(alpha, *this, CRE_DESCOMP, nroots);
  if(has(DES_DESCOMP))
    RotateOps::rotate_operator_arrays(alpha, *this, DES_DESCOMP, nroots);
  if(has(CRE_CRE_DESCOMP))
    RotateOps::rotate_operator_arrays(alpha, *this, CRE_CRE_DESCOMP, nroots);
}

};
