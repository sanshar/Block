/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "npdm_expectations_engine.h"
#include "MatrixBLAS.h"
#include <boost/format.hpp>
#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include "operatorfunctions.h"
#include "execinfo.h"

namespace SpinAdapted{
namespace Npdm{

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

double spinExpectation(Wavefunction& wave1, Wavefunction& wave2, SparseMatrix& leftOp, SparseMatrix& dotOp, SparseMatrix& rightOp, const SpinBlock& big)
{

  //calculating <wave1| Oa*Ob | wave2>
  // do transpose specifies if we want  <wave1| Oa^T*Ob |wave2> separately. This can be avoided in some sitations if wave1 and wave2 are the same functions
  int leftindices=0, dotindices=0, rightindices=0;

  leftindices = &leftOp ? leftOp.get_orbs().size() : 0;
  dotindices = &dotOp ? dotOp.get_orbs().size() : 0;
  rightindices = &rightOp ? rightOp.get_orbs().size() : 0;

  int Aindices, Bindices;
  Aindices = leftindices+dotindices;
  Bindices = rightindices;

  Wavefunction opw2;
  vector<SpinQuantum> dQ = wave1.get_deltaQuantum();
  opw2.initialise(dQ, &big, true);

  SpinBlock* leftBlock = big.get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();

  Cre AOp; //This is just an example class
  int totalspin;
  // In spin-adapted, getirrep is the value of S
  // In non spin-adapted, getirrep is the value of S_z
  if(dmrginp.spinAdapted())
    totalspin = (&rightOp) ? rightOp.get_spin().getirrep() : 0;
  else
    totalspin = (&rightOp) ? -(rightOp.get_spin().getirrep()) : 0;

  if (Aindices != 0) {
    FormLeftOp(leftBlock, leftOp, dotOp, AOp, totalspin);
  }
  
  //different cases
  if (Aindices == 0 && Bindices != 0)
  {
    operatorfunctions::TensorMultiply(rightBlock, rightOp, &big, wave2, opw2, dQ[0], 1.0);
  }
  else if (Aindices != 0 && Bindices == 0)
  { 
    operatorfunctions::TensorMultiply(leftBlock, AOp, &big, wave2, opw2, dQ[0], 1.0);
  }
  else if (Aindices != 0 && Bindices != 0)
  { 
    operatorfunctions::TensorMultiply(leftBlock, AOp, rightOp, &big, wave2, opw2, dQ[0], 1.0);
  }
  else abort();

  return DotProduct(wave1, opw2, dmrginp.Sz(), big);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void FormLeftOp(const SpinBlock* leftBlock, const SparseMatrix& leftOp, const SparseMatrix& dotOp, SparseMatrix& Aop, int totalspin)
{
  //Cre is just a class..it is not actually cre
  int leftindices=0, dotindices=0, rightindices=0;

  leftindices = &leftOp ? leftOp.get_orbs().size() : 0;
  dotindices = &dotOp ? dotOp.get_orbs().size() : 0;
  
  int Aindices, Bindices;
  Aindices = leftindices+dotindices;

  Aop.CleanUp();
  Aop.set_initialised() = true;
  Aop.set_fermion() = Aindices%2==1 ? true : false;
  if (dotindices == 0)
    { //FIXME
      // Why set fermion false? 
      //Aop.set_fermion() = false;
      Aop.set_orbs() = leftOp.get_orbs();
      Aop.set_deltaQuantum(1, leftOp.get_deltaQuantum(0)); // FIXME does leftOp always has only one dQ?
      Aop.allocate(leftBlock->get_stateInfo());
      operatorfunctions::TensorTrace(leftBlock->get_leftBlock(), leftOp, leftBlock, &(leftBlock->get_stateInfo()), Aop, 1.0);
    }
  else if (leftindices == 0)
    {
      //Aop.set_fermion() = false;
      Aop.set_orbs() = dotOp.get_orbs();
      Aop.set_deltaQuantum(1, dotOp.get_deltaQuantum(0));
      Aop.allocate(leftBlock->get_stateInfo());
      operatorfunctions::TensorTrace(leftBlock->get_rightBlock(), dotOp, leftBlock, &(leftBlock->get_stateInfo()), Aop, 1.0);
    }
  else
    {
      Aop.set_orbs() = leftOp.get_orbs(); copy(dotOp.get_orbs().begin(), dotOp.get_orbs().end(), back_inserter(Aop.set_orbs()));
      //Aop.set_fermion() = Aop.set_orbs().size() == 2 ? true : false;
      //Aop.set_fermion() = Aop.set_orbs().size() == 2 ? true : false;
      vector<SpinQuantum> spins = (dotOp.get_deltaQuantum(0) + leftOp.get_deltaQuantum(0));
      SpinQuantum dQ;
      for (int i=0; i< spins.size(); i++) {
	if (spins[i].get_s().getirrep() == totalspin) { dQ = spins[i]; break; }
      }
      if(dmrginp.spinAdapted())
        Aop.set_deltaQuantum(1, dQ);
      //FIXME
      // Above expression should also works for non-spinAdapted
      // I do not know why it does not work.
      else
        Aop.set_deltaQuantum(1, spins[0]);
      Aop.allocate(leftBlock->get_stateInfo());
      operatorfunctions::TensorProduct(leftBlock->get_leftBlock(), leftOp, dotOp, leftBlock, &(leftBlock->get_stateInfo()), Aop, 1.0);      
    }
} 

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

double DotProduct(const Wavefunction& w1, const Wavefunction& w2, double Sz, const SpinBlock& big)
{
  int leftOpSz = big.get_leftBlock()->get_stateInfo().quanta.size ();
  int rightOpSz = big.get_rightBlock()->get_stateInfo().quanta.size ();
  const StateInfo* rS = big.get_stateInfo().rightStateInfo, *lS = big.get_stateInfo().leftStateInfo;

  double output = 0.0;
  for (int lQ =0; lQ < leftOpSz; lQ++)
    for (int rQ = 0; rQ < rightOpSz; rQ++) {
      if (w1.allowed(lQ, rQ) && w2.allowed(lQ, rQ))
      {
	double b1b2 = MatrixDotProduct(w1(lQ, rQ), w2(lQ, rQ));
	output += b1b2;
       
      }	
    }

  return output;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

}
}

