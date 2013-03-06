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
#include "opxop_generic.h"
#include "wavefunction.h"
#include <boost/format.hpp>
#include "distribute.h"

#include "modules/lrt/rotateoperators.h"

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
#include "pario.h"

namespace SpinAdapted {
using namespace operatorfunctions;

// for DMRG-LRT: ( [ L(0i) x R(00) ] + [ L(i0)^(T) x R(00)^(T) ] ) * c
void SpinBlock::multiplyH_lrt_left(Wavefunction& c, Wavefunction* v, int iState, int num_threads) const
{
//pout << "DEBUG @ SpinBlock::multiplyH_lrt_left: state = " << iState << endl;
//pout << "DEBUG @ SpinBlock::multiplyH_lrt_left: leftBlock: " << endl << *leftBlock << endl;
//pout << "DEBUG @ SpinBlock::multiplyH_lrt_left: rightBlock: " << endl << *rightBlock << endl;

  Wavefunction *v_array=0, *v_distributed=0, *v_add=0;
  
  opTypes state_index_0i = make_state_index(0, iState);
  opTypes state_index_i0 = make_state_index(iState, 0);

  int maxt = 1;
  initiateMultiThread(v, v_array, v_distributed, MAX_THRD);
  dmrginp.oneelecT -> start();
  dmrginp.s0time -> start();

  if(leftBlock->has(HAM+state_index_0i)) {
    // iState >  0: Ham(L)[0, 1] x 1(R)[0, 0] * C[0]
    // iState == 0: Ham(L)[0, 0] x 1(R)[0, 0] * C[1]
    boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_array(HAM+state_index_0i).get_local_element(0)[0];
    TensorMultiply(leftBlock, *op, this, c, *v, dmrginp.effective_molecule_quantum() ,1.0, MAX_THRD);
  }

  if(iState == 0) {
    // iState >  0: 0
    // iState == 0: 1(R)[0, 0] x Ham(L)[0, 0] * C[1]
    boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_array(HAM).get_local_element(0)[0];
    TensorMultiply(rightBlock, *op, this, c, *v, dmrginp.effective_molecule_quantum(), 1.0, MAX_THRD);  
  }

  dmrginp.s0time -> stop();
#ifndef SERIAL
  boost::mpi::communicator world;
  int size = world.size();
#endif

  dmrginp.s1time -> start();

  Functor2 f;
  if(leftBlock->has(CRE_CRE_DESCOMP+state_index_0i)) {
    v_add =  leftBlock->get_op_array(CRE_CRE_DESCOMP+state_index_0i).is_local() ? v_array : v_distributed;
    f = boost::bind(&opxop::generic::cxcddcomp, leftBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), state_index_0i); 
    for_all_multithread(rightBlock->get_op_array(CRE), rightBlock->get_op_array(CRE), f);
  }

  if(leftBlock->has(CRE+state_index_0i)) {
    v_add =  rightBlock->get_op_array(CRE_CRE_DESCOMP).is_local() ? v_array : v_distributed;
    f = boost::bind(&opxop::generic::cxcddcomp, rightBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), 0); 
    for_all_multithread(leftBlock->get_op_array(CRE+state_index_0i), leftBlock->get_op_array(CRE+state_index_i0), f);  
  }

  dmrginp.s1time -> stop();

  dmrginp.oneelecT -> stop();

  dmrginp.twoelecT -> start();

  if (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY) {
    if(leftBlock->is_loopblock()) {
      dmrginp.s0time -> start();
      if(leftBlock->has(CRE_DES+state_index_0i)) {
        v_add =  rightBlock->get_op_array(CRE_DESCOMP).is_local() ? v_array : v_distributed;
        f = boost::bind(&opxop::generic::cdxcdcomp, rightBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), 0);
        for_all_multithread(leftBlock->get_op_array(CRE_DES+state_index_0i), leftBlock->get_op_array(CRE_DES+state_index_i0), f);
      }

      if(leftBlock->has(CRE_CRE+state_index_0i)) {
        v_add =  rightBlock->get_op_array(DES_DESCOMP).is_local() ? v_array : v_distributed;
        f = boost::bind(&opxop::generic::ddxcccomp, rightBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), 0);
        for_all_multithread(leftBlock->get_op_array(CRE_CRE+state_index_0i), leftBlock->get_op_array(CRE_CRE+state_index_i0), f);
      }
      dmrginp.s0time -> stop();
    }
    else {
      dmrginp.s0time -> start();
      if(leftBlock->has(CRE_DESCOMP+state_index_0i)) {
        v_add =  leftBlock->get_op_array(CRE_DESCOMP+state_index_0i).is_local() ? v_array : v_distributed;
        f = boost::bind(&opxop::generic::cdxcdcomp, leftBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), state_index_0i);
        for_all_multithread(rightBlock->get_op_array(CRE_DES), rightBlock->get_op_array(CRE_DES), f);
      }

      if(leftBlock->has(DES_DESCOMP+state_index_0i)) {
        v_add =  leftBlock->get_op_array(DES_DESCOMP+state_index_0i).is_local() ? v_array : v_distributed;
        f = boost::bind(&opxop::generic::ddxcccomp, leftBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), state_index_0i);
        for_all_multithread(rightBlock->get_op_array(CRE_CRE), rightBlock->get_op_array(CRE_CRE), f);
      }
      dmrginp.s0time -> stop();
    }
  }
  dmrginp.twoelecT -> stop();

  accumulateMultiThread(v, v_array, v_distributed, MAX_THRD);

}

void SpinBlock::multiplyH_lrt_total(Wavefunction& c, Wavefunction* v, int iState, int num_threads) const
{

  SpinBlock* loopBlock = (leftBlock->is_loopblock()) ? leftBlock : rightBlock;
  SpinBlock* otherBlock = loopBlock == leftBlock ? rightBlock : leftBlock;

  Wavefunction *v_array=0, *v_distributed=0, *v_add=0;
  
  opTypes state_index_0i = make_state_index(0, iState);
  opTypes state_index_i0 = make_state_index(iState, 0);

  int maxt = 1;
  initiateMultiThread(v, v_array, v_distributed, MAX_THRD);
  dmrginp.oneelecT -> start();
  dmrginp.s0time -> start();

  if(leftBlock->has(HAM+state_index_0i)) {
    boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_array(HAM+state_index_0i).get_local_element(0)[0];
    TensorMultiply(leftBlock, *op, this, c, *v, dmrginp.effective_molecule_quantum() ,1.0, MAX_THRD);
  }

  if(rightBlock->has(HAM+state_index_0i)) {
    boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_array(HAM+state_index_0i).get_local_element(0)[0];
    TensorMultiply(rightBlock, *op, this, c, *v, dmrginp.effective_molecule_quantum(), 1.0, MAX_THRD);  
  }

  dmrginp.s0time -> stop();
#ifndef SERIAL
  boost::mpi::communicator world;
  int size = world.size();
#endif

  dmrginp.s1time -> start();

  Functor2 f;

  if(leftBlock->has(CRE_CRE_DESCOMP+state_index_0i)) {
    v_add =  leftBlock->get_op_array(CRE_CRE_DESCOMP+state_index_0i).is_local() ? v_array : v_distributed;
    f = boost::bind(&opxop::generic::cxcddcomp, leftBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), state_index_0i); 
    for_all_multithread(rightBlock->get_op_array(CRE), rightBlock->get_op_array(CRE), f);
  }

  if(leftBlock->has(CRE+state_index_0i)) {
    v_add =  rightBlock->get_op_array(CRE_CRE_DESCOMP).is_local() ? v_array : v_distributed;
    f = boost::bind(&opxop::generic::cxcddcomp, rightBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), 0); 
    for_all_multithread(leftBlock->get_op_array(CRE+state_index_0i), leftBlock->get_op_array(CRE+state_index_i0), f);  
  }

  if(iState > 0) {

  if(rightBlock->has(CRE+state_index_0i)) {
    v_add =  leftBlock->get_op_array(CRE_CRE_DESCOMP).is_local() ? v_array : v_distributed;
    f = boost::bind(&opxop::generic::cxcddcomp, leftBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), 0); 
    for_all_multithread(rightBlock->get_op_array(CRE+state_index_0i), rightBlock->get_op_array(CRE+state_index_i0), f);
  }

  if(rightBlock->has(CRE_CRE_DESCOMP+state_index_0i)) {
    v_add =  rightBlock->get_op_array(CRE_CRE_DESCOMP+state_index_0i).is_local() ? v_array : v_distributed;
    f = boost::bind(&opxop::generic::cxcddcomp, rightBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), state_index_0i); 
    for_all_multithread(leftBlock->get_op_array(CRE), leftBlock->get_op_array(CRE), f);  
  }

  }

  dmrginp.s1time -> stop();

  dmrginp.oneelecT -> stop();

  dmrginp.twoelecT -> start();

  if (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY) {
    dmrginp.s0time -> start();

    if(loopBlock->has(CRE_DES+state_index_0i)) {
      v_add =  otherBlock->get_op_array(CRE_DESCOMP).is_local() ? v_array : v_distributed;
      f = boost::bind(&opxop::generic::cdxcdcomp, otherBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), 0);
      for_all_multithread(loopBlock->get_op_array(CRE_DES+state_index_0i), loopBlock->get_op_array(CRE_DES+state_index_i0), f);
    }

    if(loopBlock->has(CRE_CRE+state_index_0i)) {
      v_add =  otherBlock->get_op_array(DES_DESCOMP).is_local() ? v_array : v_distributed;
      f = boost::bind(&opxop::generic::ddxcccomp, otherBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), 0);
      for_all_multithread(loopBlock->get_op_array(CRE_CRE+state_index_0i), loopBlock->get_op_array(CRE_CRE+state_index_i0), f);
    }

    if(iState > 0) {

    if(otherBlock->has(CRE_DESCOMP+state_index_0i)) {
      v_add =  otherBlock->get_op_array(CRE_DESCOMP+state_index_0i).is_local() ? v_array : v_distributed;
      f = boost::bind(&opxop::generic::cdxcdcomp, otherBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), state_index_0i);
      for_all_multithread(loopBlock->get_op_array(CRE_DES), loopBlock->get_op_array(CRE_DES), f);
    }

    if(otherBlock->has(DES_DESCOMP+state_index_0i)) {
      v_add =  otherBlock->get_op_array(DES_DESCOMP+state_index_0i).is_local() ? v_array : v_distributed;
      f = boost::bind(&opxop::generic::ddxcccomp, otherBlock, _1, _2, this, ref(c), v_add, dmrginp.effective_molecule_quantum(), state_index_0i);
      for_all_multithread(loopBlock->get_op_array(CRE_CRE), loopBlock->get_op_array(CRE_CRE), f);
    }

    }

    dmrginp.s0time -> stop();
  }
  dmrginp.twoelecT -> stop();

  accumulateMultiThread(v, v_array, v_distributed, MAX_THRD);
}

void SpinBlock::rotatebyRitzVectors(const Matrix& alpha, int nroots)
{
// DEBUG
//for(int i = 0; i < nroots; ++i) {
//  pout << "DEBUG @ rotatebyRitzVectors: state = " << i << endl;
//  if(has(HAM            +make_state_index(0, i))) pout << "\t\t has HAM" << endl;
//  if(has(CRE            +make_state_index(0, i))) pout << "\t\t has CRE" << endl;
//  if(has(CRE_DES        +make_state_index(0, i))) pout << "\t\t has CRE_DES" << endl;
//  if(has(CRE_CRE        +make_state_index(0, i))) pout << "\t\t has CRE_CRE" << endl;
//  if(has(CRE_DESCOMP    +make_state_index(0, i))) pout << "\t\t has CRE_DESCOMP" << endl;
//  if(has(DES_DESCOMP    +make_state_index(0, i))) pout << "\t\t has DES_DESCOMP" << endl;
//  if(has(CRE_CRE_DESCOMP+make_state_index(0, i))) pout << "\t\t has CRE_CRE_DESCOMP" << endl;
//}
// DEBUG
  if(has(HAM)) {
    RotateOps::rotate_operator_arrays(alpha, *this, HAM, nroots, false);
  }
  if(has(CRE)) {
    RotateOps::rotate_operator_arrays(alpha, *this, CRE, nroots, false);
    RotateOps::rotate_operator_arrays(alpha, *this, CRE, nroots, true);
  }
  if(has(CRE_CRE_DESCOMP)) {
    RotateOps::rotate_operator_arrays(alpha, *this, CRE_CRE_DESCOMP, nroots, false);
    RotateOps::rotate_operator_arrays(alpha, *this, CRE_CRE_DESCOMP, nroots, true);
  }

  if(dmrginp.hamiltonian() == QUANTUM_CHEMISTRY) {

    if(has(CRE_DES)) {
      RotateOps::rotate_operator_arrays(alpha, *this, CRE_DES, nroots, false);
      RotateOps::rotate_operator_arrays(alpha, *this, CRE_DES, nroots, true);
    }
    if(has(CRE_CRE)) {
      RotateOps::rotate_operator_arrays(alpha, *this, CRE_CRE, nroots, false);
      RotateOps::rotate_operator_arrays(alpha, *this, CRE_CRE, nroots, true);
    }
    if(has(CRE_DESCOMP)) {
      RotateOps::rotate_operator_arrays(alpha, *this, CRE_DESCOMP, nroots, false);
      RotateOps::rotate_operator_arrays(alpha, *this, CRE_DESCOMP, nroots, true);
    }
    if(has(DES_DESCOMP)) {
      RotateOps::rotate_operator_arrays(alpha, *this, DES_DESCOMP, nroots, false);
      RotateOps::rotate_operator_arrays(alpha, *this, DES_DESCOMP, nroots, true);
    }

  }
}

};
