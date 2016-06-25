
#include "sweep_gen_nevpt2.h"
#include "pario.h"
#include "density.h"
#include "nevpt2_renormalize.h"
#include "guess_wavefunction.h"
#include "nevpt2_mpi.h"
#include "nevpt2_util.h"
#include "nevpt2_info.h"
#include "nevpt2_operators.h"

#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif



namespace SpinAdapted{
  
  void nevpt2::BlockAndDecimate_(SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, 
                                 const bool &useSlater, const bool& dot_with_sys,IntegralContainer &IKJL, 
                                 IntegralContainer &IKJA, IntegralContainer &IAJB, IntegralContainer &IJKA,
                                 NEVPT2Info &Info)
{
  if (dmrginp.outputlevel() > 0) 
  mcheck("at the start of block and decimate");
  // figure out if we are going forward or backwards
  p1out << "\t\t\t Performing Blocking"<<endl;
  dmrginp.guessgenT -> start();
  bool forward = (system.get_sites() [0] == 0);
  SpinBlock systemDot;
  int systemDotStart, systemDotEnd;
  int systemDotSize = sweepParams.get_sys_add() - 1;
  if (forward)
  {
    systemDotStart = *system.get_sites().rbegin () + 1;
    systemDotEnd = systemDotStart + systemDotSize;
  }
  else
  {
    systemDotStart = system.get_sites() [0] - 1;
    systemDotEnd = systemDotStart - systemDotSize;
  }
  vector<int> spindotsites(2); 
  spindotsites[0] = systemDotStart;
  spindotsites[1] = systemDotEnd;
  systemDot = SpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), true);

  const int nexact = forward ? sweepParams.get_forward_starting_size() : sweepParams.get_backward_starting_size();

  system.addAdditionalCompOps();
  InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, sweepParams.current_root(), sweepParams.current_root(), 
                                 sweepParams.get_sys_add(), dmrginp.direct(), system.get_integralIndex(), DISTRIBUTED_STORAGE_FOR_ONEPDM, true, true);
  

  pout << "\t\t\t System  Block"<<newSystem;
  newSystem.printOperatorSummary();

  std::vector<Matrix> rotateMatrix;


  //this should be done when we actually have wavefunctions stored, otherwise not!!
  SpinBlock environment, environmentDot, newEnvironment;
  int environmentDotStart, environmentDotEnd, environmentStart, environmentEnd;
  InitBlocks::InitNewEnvironmentBlock(environment, systemDot, newEnvironment, system, systemDot,
                                      sweepParams.current_root(), sweepParams.current_root(),
                                      sweepParams.get_sys_add(), sweepParams.get_env_add(), forward, dmrginp.direct(),
                                      sweepParams.get_onedot(), nexact, useSlater, system.get_integralIndex(), true, true, true);
  
  newSystem.set_loopblock(true);
  system.set_loopblock(false);
  newEnvironment.set_loopblock(false);
  SpinBlock big;
  InitBlocks::InitBigBlock(newSystem, newEnvironment, big); 
  
  const int nroots = dmrginp.nroots();
  DiagonalMatrix e;
  std::vector<Wavefunction> solutions(nroots);

  GuessWave::guess_wavefunctions(solutions, e, big, sweepParams.get_guesstype(), true, true, 0.0);
  char msg[512];
  for (int istate=0;istate<nroots;istate++){
    solutions[istate].SaveWavefunctionInfo (big.get_stateInfo(), big.get_leftBlock()->get_sites(), istate);
  }

#ifndef SERIAL
  mpi::communicator world;
  mpi::broadcast(world, solutions, 0);
#endif

  int SweepIter = sweepParams.get_sweep_iter();
  int BlockIter = sweepParams.get_block_iter();
  DensityMatrix tracedMatrix;
  tracedMatrix.allocate(newSystem.get_stateInfo());
  tracedMatrix.makedensitymatrix(solutions, big, dmrginp.weights(SweepIter), 0.0, 0.0, false);
  ThreeIndOpArray Dummy,Dummy_;
  //AddFOISDensity(big,tracedMatrix,solutions,Info,sweepParams.get_sweep_iter(),sweepParams.get_block_iter());
  AddFOISDensity(big,tracedMatrix,solutions,Dummy,Dummy_,IKJL,IKJA,IAJB,IJKA,Info,SweepIter,BlockIter);
  Dummy.Clear();
  Dummy_.Clear();
  rotateMatrix.clear();
  if (!mpigetrank()){
    double error = makeRotateMatrix(tracedMatrix, rotateMatrix, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());    
  }
  NEVPT2_AddToRotationMat(big,rotateMatrix,newSystem,solutions,sweepParams.get_block_iter());
 
#ifndef SERIAL
  mpi::broadcast(world, rotateMatrix, 0);
#endif

  for (int istate=0;istate<nroots;istate++){
    SaveRotationMatrix (newSystem.get_sites(), rotateMatrix, istate);
  }

  p1out <<"\t\t\t Performing Renormalization "<<endl;
  newSystem.transform_operators(rotateMatrix);
  if (dmrginp.outputlevel() > 0) 
    mcheck("after rotation and transformation of block");
  p2out <<newSystem<<endl;
  newSystem.printOperatorSummary();
  //mcheck("After renorm transform");
  
  //clear up memory
  for (int p=0;p<solutions.size();p++){
    solutions[p].Clear();
  }
  solutions.clear();
  
}

double nevpt2::do_one_(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize)
{

  int IntegralIndex =0;
  
  //---------------
  //the NEVPT2 Info
  //---------------
  int OrbWin[6];
  char BaseName[512];
  char msg[512];
  NEVPT2Info Info;
  Info.ReadData();
  Info.ReadOrbOrder();
  Info.ReadOrbEnergies();
  Info.GetOrbWin(OrbWin);
  Info.getBaseName(BaseName);
  Info.CalcMaxBlockIter(sweepParams.get_n_iters());
  Info.ResetNevSweep();
  int NInternal=OrbWin[1]-OrbWin[0]+1;  
  int NActive  =OrbWin[3]-OrbWin[2]+1;
  int NExternal=OrbWin[5]-OrbWin[4]+1;
  int NOcc = NInternal + NActive;
  int OrbDim = NInternal+NActive+NExternal;
  int n;//the number of possible spin couplings for V(i) and V(a) with |psi>
  int MaxCore = Info.GetMaxCore();
  
  //----------------------
  //the integral container
  //----------------------
  Matrix h,h_eff,h_eff_;
  h.ReSize(OrbDim,OrbDim);
  h_eff.ReSize(OrbDim,OrbDim);
  h_eff_.ReSize(OrbDim,OrbDim);
  vector<int> ReOrder;
  Info.GetOrbOrder(ReOrder);
  IntegralContainer IJKL(NOcc,NOcc,NOcc,NOcc,_COULOMB_);
  IntegralContainer IKJL(NOcc,NOcc,NOcc,NOcc,_EXCHANGE_);
  IntegralContainer IAJB(NOcc,NOcc,NExternal,NExternal,_EXCHANGE_);
  IntegralContainer IJAB(NOcc,NOcc,NExternal,NExternal,_COULOMB_);
  IntegralContainer IJKA(NOcc,NOcc,NOcc,NExternal,_NO_SYMM_);
  IntegralContainer IKJA(NOcc,NOcc,NOcc,NExternal,_NO_SYMM_);
  sprintf(msg,"%s.block.ijkl.tmp",BaseName);IJKL.SetFileName(msg);
  sprintf(msg,"%s.block.ikjl.tmp",BaseName);IKJL.SetFileName(msg);
  sprintf(msg,"%s.block.iajb.tmp",BaseName);IAJB.SetFileName(msg);
  sprintf(msg,"%s.block.ijab.tmp",BaseName);IJAB.SetFileName(msg);
  sprintf(msg,"%s.block.ijka.tmp",BaseName);IJKA.SetFileName(msg);
  sprintf(msg,"%s.block.ikja.tmp",BaseName);IKJA.SetFileName(msg);
  if (dmrginp.nevpt2()) ReadIntegrals(IJKL,IKJL,IAJB,IJAB,IJKA,IKJA,h,OrbWin,BaseName,ReOrder);
    
  //-----------------------------------------------------
  //construct or read the effective one-electron matrices
  //-----------------------------------------------------
  GenerateHeff(OrbWin,h,h_eff,h_eff_,IJKL,IKJL,IJAB,IAJB,IJKA,IKJA);
  //ReadHeff(OrbWin,h_eff,h_eff_,BaseName,ReOrder);
  Info.AddH(h);Info.AddH(h_eff);Info.AddH(h_eff_);
  
  //release some memory
  IJKL.Clear();
  IJAB.Clear();
  
  SpinBlock system;
  const int nroots = dmrginp.nroots();
  std::vector<double> finalEnergy(nroots,0.);
  std::vector<double> finalEnergy_spins(nroots,0.);
  double finalError = 0.;

  sweepParams.set_sweep_parameters();
  // a new renormalisation sweep routine
  pout << ((forward) ? "\t\t\t Starting renormalisation sweep in forwards direction" : "\t\t\t Starting renormalisation sweep in backwards direction") << endl;
  pout << "\t\t\t ============================================================================ " << endl;
  
  InitBlocks::InitStartingBlock (system,forward, sweepParams.current_root(), 
                                 sweepParams.current_root(), sweepParams.get_forward_starting_size(), 
                                 sweepParams.get_backward_starting_size(), restartSize, restart, warmUp, IntegralIndex);
  if(!restart)
    sweepParams.set_block_iter() = 0;

  pout << "\t\t\t Starting block is :: " << endl << system << endl;
  if (!restart) 
    SpinBlock::store (forward, system.get_sites(), system, sweepParams.current_root(), sweepParams.current_root()); // if restart, just restoring an existing block --
  sweepParams.savestate(forward, system.get_sites().size());
  bool dot_with_sys = true;
  /*if (restart)
  {
    if (forward && system.get_complementary_sites()[0] >= dmrginp.last_site()/2)
      dot_with_sys = false;
    if (!forward && system.get_sites()[0]-1 < dmrginp.last_site()/2)
      dot_with_sys = false;
  }
*/
  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); )
    {
      pout << "\n\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ----------------------------" << endl;
      if (forward)
	{ p1out << "\t\t\t Current direction is :: Forwards " << endl; }
      else
	{ p1out << "\t\t\t Current direction is :: Backwards " << endl; }
      //if (SHOW_MORE) pout << "system block" << endl << system << endl;
  
      if (dmrginp.no_transform())
	      sweepParams.set_guesstype() = BASIC;
      else if (!warmUp && sweepParams.get_block_iter() != 0) 
  	    sweepParams.set_guesstype() = TRANSFORM;
      else if (!warmUp && sweepParams.get_block_iter() == 0 && 
                ((dmrginp.algorithm_method() == TWODOT_TO_ONEDOT && dmrginp.twodot_to_onedot_iter() != sweepParams.get_sweep_iter()) ||
                  dmrginp.algorithm_method() != TWODOT_TO_ONEDOT))
        sweepParams.set_guesstype() = TRANSPOSE;
      else
        sweepParams.set_guesstype() = BASIC;
      
	p1out << "\t\t\t Blocking and Decimating " << endl;
	  
      SpinBlock newSystem;

      nevpt2::BlockAndDecimate_(sweepParams, system, newSystem, warmUp, dot_with_sys,IKJL,IKJA,IAJB,IJKA,Info);

      
      system = newSystem;

      //system size is going to be less than environment size
      /*
      if (forward && system.get_complementary_sites()[0] >= dmrginp.last_site()/2)
	dot_with_sys = false;
      if (!forward && system.get_sites()[0]-1 < dmrginp.last_site()/2)
	dot_with_sys = false;
       */
      SpinBlock::store (forward, system.get_sites(), system, sweepParams.current_root(), sweepParams.current_root());	 	

      p1out << "\t\t\t saving state " << system.get_sites().size() << endl;
      ++sweepParams.set_block_iter();
      //if (sweepParams.get_onedot())
      //pout << "\t\t\tUsing one dot algorithm!!"<<endl; 
      sweepParams.savestate(forward, system.get_sites().size());
    }
  pout << "\t\t\t Finished Generate-Blocks Sweep. " << endl;
  pout << "\t\t\t ============================================================================ " << endl;

  //clear disk and memory
  IKJL.Clear();
  IKJA.Clear();
  IAJB.Clear();
  IJKA.Clear();

  // update the static number of iterations

  ++sweepParams.set_sweep_iter();

  return finalEnergy[0];
}

  
  
  
}




