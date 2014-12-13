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
  Wavefunction opw2;
  vector<SpinQuantum> dQ = wave1.get_deltaQuantum();
  if(dmrginp.setStateSpecific() || !dmrginp.doimplicitTranspose()) opw2.initialisebra(dQ, &big, true);
  else opw2.initialise(dQ, &big, true);
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

  FormLeftOp(leftBlock, leftOp, dotOp, AOp, totalspin);
  SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));

  operatorfunctions::TensorMultiply(leftBlock, AOp, rightOp, &big, wave2, opw2, hq, 1.0);

  return DotProduct(wave1, opw2, big) ;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void FormLeftOp(const SpinBlock* leftBlock, const SparseMatrix& leftOp, const SparseMatrix& dotOp, SparseMatrix& Aop, int totalspin)
{
  //Cre is just a class..it is not actually cre
  int leftindices=0, dotindices=0;

  leftindices = leftOp.get_orbs().size() ;
  dotindices =  dotOp.get_orbs().size() ;
  
  int Aindices;
  Aindices = leftindices+dotindices;

  Aop.CleanUp();
  Aop.set_initialised() = true;
  Aop.set_fermion() = Aindices%2==1 ? true : false;
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
  //Aop.allocate(leftBlock->get_stateInfo());
  Aop.allocate(leftBlock->get_braStateInfo(),leftBlock->get_ketStateInfo());
  operatorfunctions::TensorProduct(leftBlock->get_leftBlock(), leftOp, dotOp, leftBlock, &(leftBlock->get_stateInfo()), Aop, 1.0);      
} 

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

double DotProduct(const Wavefunction& w1, const Wavefunction& w2, const SpinBlock& big)
{
  // After multipling ket by cd operator, it has the same basis with ket
  int leftOpSz = big.get_leftBlock()->get_braStateInfo().quanta.size ();
  int rightOpSz = big.get_rightBlock()->get_braStateInfo().quanta.size ();
  const StateInfo* rS = big.get_braStateInfo().rightStateInfo, *lS = big.get_braStateInfo().leftStateInfo;
  double output = 0.0;
  SpinQuantum Q= w1.get_deltaQuantum(0);
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

