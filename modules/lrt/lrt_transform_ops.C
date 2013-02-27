/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

// Written by N.N. with reference to Jon's LRT code

#include "spinblock.h"
#include "wavefunction.h"
#include <boost/format.hpp>
#include <fstream>
#include <stdio.h>

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
#include "pario.h"
namespace SpinAdapted{


//void SpinBlock::remove_normal_ops()
//{
//  ops.erase(CRE_CRE);
//  ops.erase(CRE_DES);
//}

// rotateMatrices: rotation matrix for each state
void SpinBlock::transform_operators_lrt(std::vector< std::vector<Matrix> >& rotateMatrices)
{
  StateInfo oldStateInfo = stateInfo;
  std::vector<SpinQuantum> newQuanta;
  std::vector<int> newQuantaStates;
  std::vector<int> newQuantaMap;
  for (int Q = 0; Q < rotateMatrices[0].size (); ++Q)
    {
      if (rotateMatrices[0] [Q].Ncols () != 0)
        {
          newQuanta.push_back (stateInfo.quanta [Q]);
          newQuantaStates.push_back (rotateMatrices[0] [Q].Ncols ());
          newQuantaMap.push_back (Q);
        }
    }
  StateInfo newStateInfo = StateInfo (newQuanta, newQuantaStates, newQuantaMap);

  for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it) {
    int iState = (it->first & BRAROOT_MASK) >> MASROOT_BITS;
    int jState = (it->first & KETROOT_MASK);
    opTypes ot_no_deriv = it->first & OP_TYPE_MASK;
    if      (iState == 0 && jState >  0) {
      if (!it->second->is_core()) {
        for_all_operators_multithread(*it->second, this->get_op_array(ot_no_deriv),
          bind(&SparseMatrix::build_and_renormalise_transform_deriv, _1, _2, this, it->first, ot_no_deriv,
               ref(rotateMatrices[iState]), ref(rotateMatrices[jState]), false, &newStateInfo));
      }
    }
    else if (iState >  0 && jState == 0) {
      if (!it->second->is_core()) {
        for_all_operators_multithread(*it->second, this->get_op_array(ot_no_deriv),
          bind(&SparseMatrix::build_and_renormalise_transform_deriv, _1, _2, this, it->first, ot_no_deriv,
               ref(rotateMatrices[jState]), ref(rotateMatrices[iState]), true,  &newStateInfo));
      }
    }
  }
  stateInfo = newStateInfo;
  stateInfo.AllocatePreviousStateInfo ();
  *stateInfo.previousStateInfo = oldStateInfo;

  for (int i = 0; i < newQuantaMap.size (); ++i)
    assert (stateInfo.quanta [i] == oldStateInfo.quanta [newQuantaMap [i]]);

//if (dmrginp.outputlevel() > 0) {
//  pout << "\t\t\t total elapsed time " << globaltimer.totalwalltime() << " " << globaltimer.totalcputime() << " ... " 
//       << globaltimer.elapsedwalltime() << " " << globaltimer.elapsedcputime() << endl;
//  pout << "\t\t\t Transforming to new basis " << endl;
//}
//Timer transformtimer;

  for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it) {
    int iState = (it->first & BRAROOT_MASK) >> MASROOT_BITS;
    int jState = (it->first & KETROOT_MASK);
    opTypes ot_no_deriv = it->first & OP_TYPE_MASK;
    if      (iState == 0 && jState >  0) {
      if ( it->second->is_core()) {
        for_all_operators_multithread(*it->second, this->get_op_array(ot_no_deriv),
          bind(&SparseMatrix::renormalise_transform_deriv, _1, _2,
               ref(rotateMatrices[iState]), ref(rotateMatrices[jState]), false, &newStateInfo));
      }
    }
    else if (iState >  0 && jState == 0) {
      if ( it->second->is_core()) {
        for_all_operators_multithread(*it->second, this->get_op_array(ot_no_deriv),
          bind(&SparseMatrix::renormalise_transform_deriv, _1, _2,
               ref(rotateMatrices[jState]), ref(rotateMatrices[iState]), true,  &newStateInfo));
      }
    }
  }

//for (std::map<opTypes, boost::shared_ptr<Op_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
//  if (! it->second->is_core())
//    ops[it->first]->set_core(true);

//this->direct = false;
//if (dmrginp.outputlevel() > 0)
//  pout << "\t\t\t transform time " << transformtimer.elapsedwalltime() << " " << transformtimer.elapsedcputime() << endl;

//if (leftBlock)
//  leftBlock->clear();
//if (rightBlock)
//  rightBlock->clear();
}

}
