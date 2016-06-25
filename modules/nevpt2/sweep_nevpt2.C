#include "sweep_nevpt2.h"
#include "sweep_gen_nevpt2.h"
#include "sweepgenblock.h"
#include "pario.h"
#include "density.h"
#include "nevpt2.h"
#include "nevpt2_operators.h"
#include "nevpt2_util.h"
#include "nevpt2_mpi.h"
#include "nevpt2_opconstruct.h"
#include "nevpt2_info.h"
#include "nevpt2_renormalize.h"
#include "ripdm.h"
#include "../onepdm/sweeponepdm.h"
#include "../twopdm/sweeptwopdm.h"


void dmrg(double sweep_tol);
void restart(double sweep_tol, bool reset_iter);
void ReadInput(char* conf);
void fullrestartGenblock();


namespace SpinAdapted{
  //============================================================================
  // Block and Decimate function 
  //============================================================================
  //Note: The second set of OperatorArrays will be used only in the second sweep.
  //In this case it contains the operators generated in the first sweep. For the first
  //sweep these are just dummys
  void nevpt2::BlockAndDecimate (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem,
                                 const bool &useSlater, const bool& dot_with_sys,
                                 ThreeIndOpArray &CCD, ThreeIndOpArray &CDD, IntegralContainer &IKJL,
                                 IntegralContainer &IKJA, IntegralContainer &IAJB, 
                                 IntegralContainer &IJKA, NEVPT2Info &Info)
{
  //mcheck("at the start of block and decimate");
  // figure out if we are going forward or backwards
  dmrginp.guessgenT -> start();
  bool forward = (system.get_sites() [0] == 0);
  SpinBlock systemDot;
  SpinBlock envDot;
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

  SpinBlock environment, environmentDot, newEnvironment;

  int environmentDotStart, environmentDotEnd, environmentStart, environmentEnd;

  const int nexact = forward ? sweepParams.get_forward_starting_size() : sweepParams.get_backward_starting_size();

  system.addAdditionalCompOps();
  //initialize the new system block. Note: SEMI_DISTRIBUTED_STORAGE means that 
  //all one-index operators are available for all processes
  InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, sweepParams.current_root(), sweepParams.current_root(), 
                                 sweepParams.get_sys_add(), dmrginp.direct(), system.get_integralIndex(), DISTRIBUTED_STORAGE_FOR_ONEPDM, true, true);
  InitBlocks::InitNewEnvironmentBlock(environment, systemDot, newEnvironment, system, systemDot,
                                      sweepParams.current_root(), sweepParams.current_root(),
				      sweepParams.get_sys_add(), sweepParams.get_env_add(), forward, dmrginp.direct(),
				      sweepParams.get_onedot(), nexact, useSlater, system.get_integralIndex(), true, true, true);
  
  newEnvironment.addAdditionalCompOps();
  
  SpinBlock big;
  newSystem.set_loopblock(true);
  system.set_loopblock(false);
  newEnvironment.set_loopblock(false);
  InitBlocks::InitBigBlock(newSystem, newEnvironment, big); 

  const int nroots = dmrginp.nroots();
  std::vector<Wavefunction> solutions(nroots);

  DiagonalMatrix e;
  GuessWave::guess_wavefunctions(solutions, e, big, sweepParams.get_guesstype(), true, true, 0.0); 
  
#ifndef SERIAL
  mpi::communicator world;
  mpi::broadcast(world, solutions, 0);
#endif

  

#ifdef SERIAL
  const int numprocs = 1;
#endif
#ifndef SERIAL
  const int numprocs = world.size();
#endif
  
  //----------------------------------------------------------------------------
  // perform the majority of the NEVPT2 calculations
  //----------------------------------------------------------------------------
  char msg[512];
  int iter=sweepParams.get_block_iter();
  bool FirstNevpt2Sweep = (Info.GetNevSweep()==0);
  int iSweep = Info.GetNevSweep();
  int MaxIter=Info.GetMaxBlockIter(iSweep);
  bool NonConvNevpt2 = ((dmrginp.nevpt2())&&(!dmrginp.read_higherpdm()));
  if (FirstNevpt2Sweep){
    if (sweepParams.get_block_iter()==MaxIter){
      if (!dmrginp.nevpt2()) generate_RI_density_matrices(solutions, big);
      if (dmrginp.nevpt2())  NEVPT2_Driver(big, solutions, Info);
    }//the middle of the sweep
  }//first NEVPT2 sweep

  //----------------------------------------------------------------------------
  // Construct the three-index operators required for V(i) and V(a)
  //----------------------------------------------------------------------------
  if (NonConvNevpt2){
    if (iter<=MaxIter){
      //-------------------------------------------------
      //set the current StateInfo for the Operator Arrays
      //-------------------------------------------------
      CCD.SetStateInfo(big.get_stateInfo());
      CDD.SetStateInfo(big.get_stateInfo());
      //---------------------------------------------
      //actually construct the operator contributions
      //---------------------------------------------
      ConstructCreCreDes(big,solutions, CCD, IKJL, sweepParams, Info);
      ConstructCreDesDes(big,solutions, CDD, IKJA, sweepParams, Info);
    };
  }//nevpt2?

  std::vector<Matrix> rotateMatrix;
  DensityMatrix tracedMatrix;
  tracedMatrix.allocate(newSystem.get_stateInfo());
  tracedMatrix.makedensitymatrix(solutions, big, dmrginp.weights(sweepParams.get_sweep_iter()), 0.0, 0.0, false);
  //AddFOISDensity(big,tracedMatrix,solutions,Info,sweepParams.get_sweep_iter(),sweepParams.get_block_iter());
  AddFOISDensity(big,tracedMatrix,solutions,CCD,CDD,IKJL,IKJA,IAJB,IJKA,Info,sweepParams.get_sweep_iter(),sweepParams.get_block_iter());
  rotateMatrix.clear();
  if (!mpigetrank()){
    double error = makeRotateMatrix(tracedMatrix, rotateMatrix, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());
  }
  NEVPT2_AddToRotationMat(big,rotateMatrix,newSystem,solutions,sweepParams.get_block_iter());


#ifndef SERIAL
  mpi::broadcast(world,rotateMatrix,0);
#endif
  //----------------------------------------------------------------------------
  //Save rotation matrices and transform operators
  //----------------------------------------------------------------------------
  //save the rotation matrix and wavefunction from this iteration 
  for (int iroot=0;iroot<nroots;iroot++){
    SaveRotationMatrix (newSystem.get_sites(), rotateMatrix, iroot);
    solutions[iroot].SaveWavefunctionInfo (big.get_stateInfo(), big.get_leftBlock()->get_sites(), iroot);
  }

  mpi_barrier();
  //TensorTrace all operators from previous block iterations
  if (iter<=MaxIter&&NonConvNevpt2) {
    double StartTraceTime = GetTime();
    CDD.TensorTrace(big,newSystem.get_stateInfo());
    CCD.TensorTrace(big,newSystem.get_stateInfo());
    double FinishTraceTime = GetTime();
    double TraceTime = FinishTraceTime-StartTraceTime;
    Info.AddTime(3,TraceTime);
  }

  mpi_barrier();
  //transform the operators of the new system block and update the StateInfo
  newSystem.transform_operators(rotateMatrix);

  //transform the Vtuv operators
  if (iter<=MaxIter&&NonConvNevpt2) {
    double StartTrafoTime = GetTime();
    CDD.RenormalizeTransform(big,newSystem.get_stateInfo(),rotateMatrix);
    CCD.RenormalizeTransform(big,newSystem.get_stateInfo(),rotateMatrix);
    double FinishTrafoTime = GetTime();
    double TrafoTime = FinishTrafoTime-StartTrafoTime;
    Info.AddTime(3,TrafoTime);
  }
  mpi_barrier();
  
  //clear up memory
  for (int p=0;p<solutions.size();p++){
    solutions[p].Clear();
  }
  solutions.clear();
  
}

  
  
  
  
  
  //============================================================================
  // Sweep Iteration function reproduced form the one-particl reduced density
  // module. It serves to read all operators and generate a big SpinBlock
  //============================================================================
  double nevpt2::do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize)
{

  SpinBlock system;
  const int nroots = dmrginp.nroots();
  std::vector<double> finalEnergy(nroots,0.);
  std::vector<double> finalEnergy_spins(nroots,0.);
  double finalError = 0.;
  int IntegralIndex = 0;
  int MaxM = dmrginp.sweep_state_schedule()[dmrginp.sweep_qstate_schedule().size()];

  bool NonConvNevpt2 = ((dmrginp.nevpt2())&&(!dmrginp.read_higherpdm()));
  //if (NonConvNevpt2) mpi_debug_break();
  
  //-----------------------------------------------------
  //the storage container for the V(i) and V(a) operators
  //-----------------------------------------------------
  int OrbWin[6];
  char BaseName[512];
  char msg[512];
  //---------------
  //the NEVPT2 Info
  //---------------
  NEVPT2Info Info;
  Info.ReadData();
  Info.ReadOrbOrder();
  Info.ReadOrbEnergies();
  Info.SetNRoots(nroots);
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
  //--------------------------------
  // the three-index operator arrays
  //--------------------------------
  ThreeIndOpArray CCD,CCD_,CDD,CDD_;
  sprintf(msg,"%s.CCD.tmp",BaseName);
  CCD.Initialize(NActive,msg,mpi_world_size());
  CCD.CalcBufferSize(MaxCore,MaxM);
  sprintf(msg,"%s.CCD_.tmp",BaseName);
  CCD_.Initialize(NActive,msg,mpi_world_size());
  CCD_.CalcBufferSize(MaxCore,MaxM);
  sprintf(msg,"%s.CDD.tmp",BaseName);
  CDD.Initialize(NActive,msg,mpi_world_size());
  CDD.CalcBufferSize(MaxCore,MaxM);
  sprintf(msg,"%s.CDD_.tmp",BaseName);
  CDD_.Initialize(NActive,msg,mpi_world_size());
  CDD_.CalcBufferSize(MaxCore,MaxM);
  //invoke disk storage
  CCD.SetInCore(false);
  CCD_.SetInCore(false);
  CDD.SetInCore(false);
  CDD_.SetInCore(false);
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
  
  //----------------------
  // start the first sweep
  //----------------------
  sweepParams.set_sweep_parameters();
  // a new renormalisation sweep routine
  pout << ((forward) ? "\t\t\t Starting renormalisation sweep in forwards direction" : "\t\t\t Starting renormalisation sweep in backwards direction") << endl;
  pout << "\t\t\t ============================================================================ " << endl;
  
  InitBlocks::InitStartingBlock (system,forward, sweepParams.current_root(), sweepParams.current_root(),
                                 sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(),
                                 restartSize, restart, warmUp, IntegralIndex);
  SpinBlock::store (forward, system.get_sites(), system, sweepParams.current_root(),sweepParams.current_root()); // if restart, just restoring an existing block --
  sweepParams.savestate(forward, system.get_sites().size());
  bool dot_with_sys = true;
  pout << "\t\t\t Starting block is :: " << endl << system << endl;

  sweepParams.set_block_iter() = 0;

  //----------------------------------------------------------------------------
  //here we perform iterations of the first sweep where the contributions for 
  //V(ijab), V(ija), V(iab), V(ab), V(ij) aand V(ia) are evaluated. Furthermore,
  //the perturber functions for V(i) and V(a) are generated. However, their
  //contributions are evaluated at the middle of a second sweep 
  //----------------------------------------------------------------------------
  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); ){
    
    //print information
    pout << "\n\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
    pout << "\t\t\t ----------------------------" << endl;
    if (forward)
      { p1out << "\t\t\t Current direction is :: Forwards " << endl; }
    else
      { p1out << "\t\t\t Current direction is :: Backwards " << endl; }

    //determine the way of guessing the wavefunction
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
      
    //block and decimate. Note: when we have reached the middle position, the 
    //density generator is called
    p1out << "\t\t\t Blocking and Decimating " << endl;
    SpinBlock newSystem;
    nevpt2::BlockAndDecimate (sweepParams, system, newSystem, warmUp, dot_with_sys, CCD, CDD, IKJL, IKJA, IAJB, IJKA, Info);
    
    //update the system block
    system = newSystem;
    
    //increase the block iteration counter
    ++sweepParams.set_block_iter();
    
    //save the block
    p1out << "\t\t\t saving state " << endl;
    SpinBlock::store (forward, system.get_sites(), system, sweepParams.current_root(), sweepParams.current_root());
    sweepParams.savestate(forward, system.get_sites().size());
    
    pout.precision(12);
    pout << "\t\t\t The lowest sweep energy : "<< sweepParams.get_lowest_energy()[0]<<endl;
    pout << "\t\t\t ============================================================================ " << endl;
  }//blockiter
  
  //increment the sweep counter
  Info.IncrementNevSweep();
  
  //store the three-index operators on disk
  CDD.Store(false);
  CCD.Store(false);
  //and free the memory
  if (CDD.GetInCore()) CDD.Clear();
  if (CCD.GetInCore()) CCD.Clear();
  
  //-----------------------
  // start the second sweep
  //-----------------------
  sweepParams.set_sweep_iter()++;
  sweepParams.set_sweep_parameters();
  // a new renormalisation sweep routine
  pout << ((!forward) ? "\t\t\t Starting renormalisation sweep in forwards direction" : "\t\t\t Starting renormalisation sweep in backwards direction") << endl;
  pout << "\t\t\t ============================================================================ " << endl;
  
  InitBlocks::InitStartingBlock (system,!forward, sweepParams.current_root(), sweepParams.current_root(), 
                                 sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), 
                                 restartSize, restart, warmUp, IntegralIndex);
  SpinBlock::store (!forward, system.get_sites(), system, sweepParams.current_root(), sweepParams.current_root()); // if restart, just restoring an existing block --
  sweepParams.savestate(!forward, system.get_sites().size());
  pout << "\t\t\t Starting block is :: " << endl << system << endl;

  sweepParams.set_block_iter() = 0;

  //------------------------
  //perform the second sweep
  //------------------------
  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); ){
    
    //print information
    pout << "\n\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
    pout << "\t\t\t ----------------------------" << endl;
    if (!forward)
      { p1out << "\t\t\t Current direction is :: Forwards " << endl; }
    else
      { p1out << "\t\t\t Current direction is :: Backwards " << endl; }

    //determine the way of guessing the wavefunction
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
      
    //block and decimate. Note: when we have reached the middle position, the 
    //density generator is called
    p1out << "\t\t\t Blocking and Decimating " << endl;
    SpinBlock newSystem;
    nevpt2::BlockAndDecimate (sweepParams, system, newSystem, warmUp, dot_with_sys, CCD_, CDD_, IKJL, IKJA, IAJB, IJKA, Info);
    
    //update the system block
    system = newSystem;
    
    //increase the block iteration counter
    ++sweepParams.set_block_iter();
    
    //save the block
    p1out << "\t\t\t saving state " << endl;
    SpinBlock::store (!forward, system.get_sites(), system, sweepParams.current_root(), sweepParams.current_root());
    sweepParams.savestate(!forward, system.get_sites().size());
    
    pout.precision(12);
    pout << "\t\t\t The lowest sweep energy : "<< sweepParams.get_lowest_energy()[0]<<endl;
    pout << "\t\t\t ============================================================================ " << endl;
  }//blockiter
  //increment the sweep counter
  Info.IncrementNevSweep();
  
  //-----------------------
  //Print the NEVPT2 Output
  //-----------------------
  if (dmrginp.nevpt2()&&(!Info.ConvNevPT2())) GiveNEVPT2Output(Info);
  
  //------------------
  //clean up and leave
  //------------------
  CCD.Clear();
  CDD.Clear();
  CCD_.Clear();
  CDD_.Clear();
  IKJL.Clear();
  IKJA.Clear();
  IAJB.Clear();
  IJKA.Clear();
  
  return sweepParams.get_lowest_energy()[0];

}

  
  void nevpt2::nevpt2(){
    double sweep_tol = dmrginp.get_sweep_tol();
    SweepParams sweepParams;
    SweepParams sweep_copy;
    bool direction,direction_copy;
    int restartsize,restartsize_copy;
    sweepParams.current_root() = -1;

    if(sym != "c1") {
      pout << "NEVPT2 is only implemented with c1 symmetry"<<endl;
      abort();
    }

    if (dmrginp.algorithm_method() == TWODOT) {
      pout << "ripdm not allowed with twodot algorithm" << endl;
      abort();
    }
    if (RESTART && !FULLRESTART)
      restart(sweep_tol, reset_iter);
    else if (FULLRESTART) {
      fullrestartGenblock();
      reset_iter = true;
      sweepParams.restorestate(direction, restartsize);
      sweepParams.calc_niter();
      sweepParams.savestate(direction, restartsize);
      restart(sweep_tol, reset_iter);
    }
    else {
      dmrg(sweep_tol);
    }
    //dmrginp.screen_tol() = 0.0; //need to turn screening off for ripdm
    //dmrginp.do_cd() = true;
    dmrginp.do_npdm_ops() = true;
    dmrginp.do_pdm() = true;
    dmrginp.oneindex_screen_tol() = 0.0; 
    dmrginp.twoindex_screen_tol() = 0.0;
    dmrginp.Sz() = dmrginp.total_spin_number().getirrep();

    //first generate the 1- and 2-pdm
    sweep_copy.restorestate(direction_copy, restartsize_copy);
    sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
    dmrginp.set_fullrestart()=true;

    for (int istate=0;istate<dmrginp.nroots();istate++){
      SweepOnepdm::do_one(sweepParams, false, !direction, false, 0, istate);
    }
    
    sweep_copy.restorestate(direction_copy, restartsize_copy);
    sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
    for (int istate=0;istate<dmrginp.nroots();istate++){
      SweepTwopdm::do_one(sweepParams, false, direction,false, 0, istate);
    }
    
    //Do a sweep that generates the operators and the rotation matrices
    sweep_copy.restorestate(direction_copy, restartsize_copy);
    dmrginp.set_fullrestart() = true;
    sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
    nevpt2::do_one_(sweepParams, false, !direction, false, 0);

    // Do the NEVPT2 calculation
    dmrginp.set_fullrestart() = false;    
    nevpt2::do_one(sweepParams,false,direction,false,0);
    sweep_copy.savestate(direction_copy, restartsize_copy);
    
  }
  
  void nevpt2::nevpt2_restart(){
    
    double sweep_tol = dmrginp.get_sweep_tol();
    SweepParams sweepParams;
    SweepParams sweep_copy;
    bool direction,direction_copy;
    int restartsize,restartsize_copy;
    sweepParams.current_root() = -1;
    
    if(sym != "c1") {
      pout << "NEVPT2 is only implemented with c1 symmetry"<<endl;
      abort();
    }
    sweepParams.restorestate(direction, restartsize);
    if(!sweepParams.get_onedot() || dmrginp.algorithm_method() == TWODOT) {
      pout << "ripdm only runs for the onedot algorithm" << endl;
      abort();
    }
    
    dmrginp.do_npdm_ops() = true;
    dmrginp.do_pdm() = true;
    dmrginp.oneindex_screen_tol() = 0.0; 
    dmrginp.twoindex_screen_tol() = 0.0;
    dmrginp.Sz() = dmrginp.total_spin_number().getirrep();

    //first generate the 1- and 2-pdm
    sweep_copy.restorestate(direction_copy, restartsize_copy);
    sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
    dmrginp.set_fullrestart()=true;

    for (int istate=0;istate<dmrginp.nroots();istate++){
      SweepOnepdm::do_one(sweepParams, false, !direction, true, restartsize, istate);
    }
    
    sweep_copy.restorestate(direction_copy, restartsize_copy);
    sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
    for (int istate=0;istate<dmrginp.nroots();istate++){
      SweepTwopdm::do_one(sweepParams, false, direction,true, restartsize, istate);
    }

    //Do a sweep that generates the operators and the rotation matrices
    sweep_copy.restorestate(direction_copy, restartsize_copy);
    dmrginp.set_fullrestart() = true;
    nevpt2::do_one_(sweepParams, false, !direction, false, 0); //this will generate the cd operators and extended roation matrices
    //SweepTwopdm::do_one(sweepParams, false, direction, false, 0, 0); //this will generate the cd operators and extended roation matrices
    dmrginp.set_fullrestart() = false;    

    // Do the NEVPT2 calculation
    nevpt2::do_one(sweepParams,false,direction,true,1);
    sweep_copy.savestate(direction_copy, restartsize_copy);
  }

  
}















