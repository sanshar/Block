/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


//#include "twopdm.h"
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
  SpinQuantum dQ = wave1.get_deltaQuantum();
  opw2.initialise(dQ, &big, true);

  SpinBlock* leftBlock = big.get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();

  Cre AOp; //This is just an example class
  int totalspin = (&rightOp) ? rightOp.get_spin() : 0;
//MAWpout << "totalspin  " << totalspin << std::endl;

  if (Aindices != 0) {
//if (dotindices==3 && rightindices==1) {
//MAWpout << "Building AOp:\n";
//MAWpout << "dotOp:\n";
//MAWpout << dotOp;
//MAWpout << "rhsOp:\n"; // can't print Transposeview
//MAWpout << rightOp;
//}
    FormLeftOp(leftBlock, leftOp, dotOp, AOp, totalspin);
//MAWif (leftindices > 0) { pout << "leftOp:\n"; pout << leftOp; }
//MAWif (dotindices > 0) { pout << "dotOp:\n"; pout << dotOp; }
//pout << "AOp:\n";
//pout << AOp;
  }
  
  //different cases
  if (Aindices == 0 && Bindices != 0)
  {
//FIXME    assert( Bindices == 4 ); //2PDM
//FIXME    assert( Bindices == 6 ); //3PDM
    operatorfunctions::TensorMultiply(rightBlock, rightOp, &big, wave2, opw2, dQ, 1.0);
//MAW    expectations.push_back( DotProduct(wave1, opw2, dmrginp.Sz(), big) );
  }
  else if (Aindices != 0 && Bindices == 0)
  { 
//FIXME    assert( Aindices == 4 ); //2PDM
//FIXME    assert( Aindices == 6 ); //3PDM
    operatorfunctions::TensorMultiply(leftBlock, AOp, &big, wave2, opw2, dQ, 1.0);
//MAW    expectations.push_back( DotProduct(wave1, opw2, dmrginp.Sz(), big) );
  }
  else if (Aindices != 0 && Bindices != 0)
  { 
    operatorfunctions::TensorMultiply(leftBlock, AOp, rightOp, &big, wave2, opw2, dQ, 1.0);
//MAW    expectations.push_back( DotProduct(wave1, opw2, dmrginp.Sz(), big) );    
//MAW    if (doTranspose)
//MAW    {
//MAW      opw2.Clear();
//MAW      operatorfunctions::TensorMultiply(leftBlock, Transposeview(AOp), rightOp, &big, wave2, opw2, dQ, 1.0);
//MAW      expectations.push_back( DotProduct(wave1, opw2, dmrginp.Sz(), big) );    
//MAW    }
  }
  else assert(false);

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
  if (dotindices == 0)
    {
      Aop.set_fermion() = false;
      Aop.set_orbs() = leftOp.get_orbs();
      Aop.set_deltaQuantum() = leftOp.get_deltaQuantum();
      Aop.allocate(leftBlock->get_stateInfo());
      operatorfunctions::TensorTrace(leftBlock->get_leftBlock(), leftOp, leftBlock, &(leftBlock->get_stateInfo()), Aop, 1.0);
    }
  else if (leftindices == 0)
    {
      Aop.set_fermion() = false;
      Aop.set_orbs() = dotOp.get_orbs();
      Aop.set_deltaQuantum() = dotOp.get_deltaQuantum();
      Aop.allocate(leftBlock->get_stateInfo());
//MAWpout << "dotOp:\n";
//MAWpout << dotOp;
//MAWcout << "dotOp.get_sign():\n";
//MAWcout << dotOp.get_sign() << std::endl;
//MAWcout << "dotOp.get_initialised():\n";
//MAWcout << dotOp.get_initialised() << std::endl;
//MAWcout << "dotOp.get_fermion():\n";
//MAWcout << dotOp.get_fermion() << std::endl;
//MAWcout << "dotOp.get_spin():\n";
//MAWcout << dotOp.get_spin() << std::endl;
//MAWcout << "dotOp.get_symm():\n";
//MAWcout << dotOp.get_symm() << std::endl;
//MAWcout << "dotOp.get_deltaQuantum():\n";
//MAWcout << dotOp.get_deltaQuantum() << std::endl;
//MAWcout << "dotOp.get_orbs():\n";
//MAWcout << dotOp.get_orbs()[0] << dotOp.get_orbs()[1] << dotOp.get_orbs()[2] << std::endl;
//MAWcout << "dotOp.get_built():\n";
//MAWcout << dotOp.get_orbs()[0] << dotOp.get_orbs()[1] << dotOp.get_orbs()[2] << std::endl;
//MAWcout << "leftBlock->get_stateInfo():\n";
//MAWcout << leftBlock->get_stateInfo() << std::endl;
//MAWcout << "Aop.get_orbs():\n";
//MAWcout << Aop.get_orbs()[0] << Aop.get_orbs()[1] << Aop.get_orbs()[2] << std::endl;
      operatorfunctions::TensorTrace(leftBlock->get_rightBlock(), dotOp, leftBlock, &(leftBlock->get_stateInfo()), Aop, 1.0);
    }
  else
    {
      Aop.set_orbs() = leftOp.get_orbs(); copy(dotOp.get_orbs().begin(), dotOp.get_orbs().end(), back_inserter(Aop.set_orbs()));
      Aop.set_fermion() = Aop.set_orbs().size() == 2 ? true : false;
      vector<SpinQuantum> spins = (dotOp.get_deltaQuantum() + leftOp.get_deltaQuantum());
      SpinQuantum dQ;
      for (int i=0; i< spins.size(); i++) {
	if (spins[i].get_s() == totalspin) { dQ = spins[i]; break; }
      }
      Aop.set_deltaQuantum() = dQ;
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

