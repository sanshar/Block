/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "IntegralMatrix.h"
#include <fstream>
#include "input.h"
#include "pario.h"
#include "global.h"
#include "orbstring.h"
#include "least_squares.h"
#include <include/communicate.h>
#include "sweepgenblock.h"
#include "npdm.h"
#include "ds1_sweeponepdm.h"  //EL
#include "ds0_sweeponepdm.h"  //EL


#ifdef _OPENMP
#include <omp.h>
#endif

//the following can be removed later
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/export.hpp>
#include <boost/format.hpp>
#include "spinblock.h"
#include "StateInfo.h"
#include "operatorfunctions.h"
#include "wavefunction.h"
#include "solver.h"
#include "davidson.h"
#include "guess_wavefunction.h"
#include "rotationmat.h"
#include "density.h"
#include "sweep.h"
#include "sweepCompress.h"
#include "sweepResponse.h"
#include "BaseOperator.h"
#include "dmrg_wrapper.h"

#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include "pario.h"
#include "modules/nevpt2/sweep_nevpt2.h"
#include "mps_nevpt.h"


#ifdef USE_BTAS
void calculateOverlap();
#endif
void dmrg(double sweep_tol);
void compress(double sweep_tol, int targetState, int baseState);
void responseSweep(double sweep_tol, int targetState, vector<int>& correctionVector, vector<int>& baseState);
void restart(double sweep_tol, bool reset_iter);
void dmrg_stateSpecific(double sweep_tol, int targetState);
void ReadInput(char* conf);
void fullrestartGenblock();
void license() {
#ifndef MOLPRO
  pout << "Block  Copyright (C) 2012  Garnet K.-L. Chan"<<endl;
  pout << "This program comes with ABSOLUTELY NO WARRANTY; for details see license file."<<endl;
  pout << "This is free software, and you are welcome to redistribute it"<<endl;
  pout << "under certain conditions; see license file for details."<<endl;
#endif
}


namespace SpinAdapted{
  Timer globaltimer(false);
  double tcpu,twall,ecpu,ewall;
  bool DEBUGWAIT = false;
  bool DEBUG_MEMORY = false;
  bool restartwarm = false;
  double NUMERICAL_ZERO = 1e-15;
  double BWPTenergy = 0.0;
  std::vector<OneElectronArray> v_1;
  std::vector<TwoElectronArray> v_2;
  std::map<TwoPerturbType,PerturbTwoElectronArray> vpt_2;
  OneElectronArray vpt_1;
//  std::map<OnePerturbType,OnePerturbArray> vpt_1;
//  std::map<TwoPerturbType,TwoPerturbArray> vpt_2;
  OneElectronArray fock;
  std::vector<double> coreEnergy;
  PairArray v_cc;
  CCCCArray v_cccc;
  CCCDArray v_cccd;
  Input dmrginp;
  int MAX_THRD = 1;
  bool FULLRESTART;
  bool RESTART;
  bool BACKWARD;
  bool reset_iter;
  std::vector<int> NPROP;
  int PROPBITLEN=1;
}

using namespace SpinAdapted;

int calldmrg(char* input, char* output)
{
  //sleep(15);
  streambuf *backup;
  backup = cout.rdbuf();
  ofstream file;
  if (output != 0) {
    file.open(output);
    pout.rdbuf(file.rdbuf());
  }
  license();
  ReadInput(input);
  MAX_THRD = dmrginp.thrds_per_node()[mpigetrank()];
#ifdef _OPENMP
  omp_set_num_threads(MAX_THRD);
#endif
  pout.precision (12);

   //Initializing timer calls
  dmrginp.initCumulTimer();

  double sweep_tol = 1e-7;
  sweep_tol = dmrginp.get_sweep_tol();
  bool direction;
  int restartsize;
  SweepParams sweepParams;

  SweepParams sweep_copy;
  bool direction_copy; int restartsize_copy;
  Matrix O, H;
  
  switch(dmrginp.calc_type()) {

  case (COMPRESS):
  {
    bool direction; int restartsize;
    //sweepParams.restorestate(direction, restartsize);
    //sweepParams.set_sweep_iter() = 0;
    restartsize = 0;

    int targetState, baseState, correctionVector, firstorderstate;
    {
      direction = true;

      //base state is always defined
      baseState = dmrginp.baseStates()[0];

      //if targetstate is given use it otherwise use basestate+1
      if(dmrginp.targetState() == -1)
	targetState = dmrginp.baseStates()[0]+1;
      else
	targetState = dmrginp.targetState();

      algorithmTypes atype = dmrginp.algorithm_method();
      dmrginp.set_algorithm_method() = ONEDOT;
      //initialize state info and canonicalize wavefunction is always done using onedot algorithm
      if (mpigetrank()==0) {
	Sweep::InitializeStateInfo(sweepParams, direction, baseState);
	Sweep::InitializeStateInfo(sweepParams, !direction, baseState);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, baseState);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, baseState);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, baseState);
      }
      dmrginp.set_algorithm_method() = atype;
    }

    //this genblock is required to generate all the nontranspose operators
    dmrginp.setimplicitTranspose() = false;
    SweepGenblock::do_one(sweepParams, false, false, false, restartsize, baseState, baseState);


    compress(sweep_tol, targetState, baseState);

    break;
  }
  case (RESPONSEBW):
  {
    //compressing the V|\Psi_0>, here \Psi_0 is the basestate and 
    //its product with V will have a larger bond dimension and is being compressed
    //it is called the target state
    dmrginp.setimplicitTranspose() = false;


    sweepParams.restorestate(direction, restartsize);
    algorithmTypes atype = dmrginp.algorithm_method();
    dmrginp.set_algorithm_method() = ONEDOT;
    if (mpigetrank()==0 && !RESTART && !FULLRESTART) {
      for (int l=0; l<dmrginp.projectorStates().size(); l++) {
	Sweep::InitializeStateInfo(sweepParams, direction, dmrginp.projectorStates()[l]);
	Sweep::InitializeStateInfo(sweepParams, !direction, dmrginp.projectorStates()[l]);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, dmrginp.projectorStates()[l]);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, dmrginp.projectorStates()[l]);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, dmrginp.projectorStates()[l]);
      }
      for (int l=0; l<dmrginp.baseStates().size(); l++) {
	Sweep::InitializeStateInfo(sweepParams, direction, dmrginp.baseStates()[l]);
	Sweep::InitializeStateInfo(sweepParams, !direction, dmrginp.baseStates()[l]);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, dmrginp.baseStates()[l]);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, dmrginp.baseStates()[l]);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, dmrginp.baseStates()[l]);
      }
    }
    dmrginp.set_algorithm_method() = atype;

    
    pout << "DONE COMPRESSING THE CORRECTION VECTOR"<<endl;
    pout << "NOW WE WILL OPTIMIZE THE RESPONSE WAVEFUNCTION"<<endl;
    //finally now calculate the response state
    responseSweep(sweep_tol, dmrginp.targetState(), dmrginp.projectorStates(), dmrginp.baseStates());

    break;
  }
  case (RESPONSE):
  {
    //compressing the V|\Psi_0>, here \Psi_0 is the basestate and 
    //its product with V will have a larger bond dimension and is being compressed
    //it is called the target state
    dmrginp.setimplicitTranspose() = false;


    sweepParams.restorestate(direction, restartsize);
    algorithmTypes atype = dmrginp.algorithm_method();
    dmrginp.set_algorithm_method() = ONEDOT;
    if (mpigetrank()==0 && !RESTART && !FULLRESTART) {
      for (int l=0; l<dmrginp.projectorStates().size(); l++) {
	Sweep::InitializeStateInfo(sweepParams, direction, dmrginp.projectorStates()[l]);
	Sweep::InitializeStateInfo(sweepParams, !direction, dmrginp.projectorStates()[l]);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, dmrginp.projectorStates()[l]);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, dmrginp.projectorStates()[l]);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, dmrginp.projectorStates()[l]);
      }
      for (int l=0; l<dmrginp.baseStates().size(); l++) {
	Sweep::InitializeStateInfo(sweepParams, direction, dmrginp.baseStates()[l]);
	Sweep::InitializeStateInfo(sweepParams, !direction, dmrginp.baseStates()[l]);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, dmrginp.baseStates()[l]);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, dmrginp.baseStates()[l]);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, dmrginp.baseStates()[l]);
      }
    }
    dmrginp.set_algorithm_method() = atype;

    
    pout << "DONE COMPRESSING THE CORRECTION VECTOR"<<endl;
    pout << "NOW WE WILL OPTIMIZE THE RESPONSE WAVEFUNCTION"<<endl;
    //finally now calculate the response state
    responseSweep(sweep_tol, dmrginp.targetState(), dmrginp.projectorStates(), dmrginp.baseStates());

    break;
  }
  case (CALCOVERLAP):
  {
    pout.precision(12);
    if (mpigetrank() == 0) {
      for (int istate = 0; istate<dmrginp.nroots(); istate++) {
	bool direction;
	int restartsize;
	sweepParams.restorestate(direction, restartsize);
	Sweep::InitializeStateInfo(sweepParams, !direction, istate);
	Sweep::InitializeStateInfo(sweepParams, direction, istate);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, istate);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, istate);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, istate);
      }
      for (int istate = 0; istate<dmrginp.nroots(); istate++) 
	for (int j=istate; j<dmrginp.nroots() ; j++) {
	  int integralIndex = 0;
	  Sweep::InitializeOverlapSpinBlocks(sweepParams, !direction, j, istate, integralIndex);
	  Sweep::InitializeOverlapSpinBlocks(sweepParams, direction, j, istate, integralIndex);
	}
      //Sweep::calculateAllOverlap(O);
    }
    break;
  }
  case (CALCHAMILTONIAN):
  {
    pout.precision(12);

    for (int istate = 0; istate<dmrginp.nroots(); istate++) {
      bool direction;
      int restartsize;
      sweepParams.restorestate(direction, restartsize);
      
      if (mpigetrank() == 0) {
	Sweep::InitializeStateInfo(sweepParams, !direction, istate);
	Sweep::InitializeStateInfo(sweepParams, direction, istate);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, istate);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, istate);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, istate);
      }
    }
    
    //Sweep::calculateHMatrixElements(H);
    pout << "overlap "<<endl<<O<<endl;
    pout << "hamiltonian "<<endl<<H<<endl;
    break;
  }
  case (DMRG):
  {
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
    else if (BACKWARD) {
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
    break;
  }
  case (FCI):
    Sweep::fullci(sweep_tol);
    break;
    
  case (TINYCALC):
    Sweep::tiny(sweep_tol);
    break;
  case (ONEPDM):
    Npdm::npdm(NPDM_ONEPDM);
    if (dmrginp.hamiltonian() == BCS) {
      Npdm::npdm(NPDM_PAIRMATRIX,true);
    }
    break;

  case (TWOPDM):
    Npdm::npdm(NPDM_TWOPDM);
    break;

  case (THREEPDM):
    Npdm::npdm(NPDM_THREEPDM);
    break;

  case (FOURPDM):
    Npdm::npdm(NPDM_FOURPDM);
    break;

  case (NEVPT2PDM):
    Npdm::npdm(NPDM_NEVPT2);
    break;

  case(NEVPT2):
    nevpt2::nevpt2();
    break;

  case(MPS_NEVPT):
    mps_nevpt::mps_nevpt(sweep_tol);
    break;
    
  case(RESTART_MPS_NEVPT):
    mps_nevpt::mps_nevpt(sweep_tol);
    break;
    
  case (RESTART_ONEPDM):
    Npdm::npdm(NPDM_ONEPDM,true);
    if (dmrginp.hamiltonian() == BCS) {
      Npdm::npdm(NPDM_PAIRMATRIX,true);
    }
    break;

  case (RESTART_TWOPDM):
    Npdm::npdm(NPDM_TWOPDM,true);
    break;
  case (RESTART_THREEPDM):
    Npdm::npdm(NPDM_THREEPDM,true);
    break;
  case (RESTART_FOURPDM):
    Npdm::npdm(NPDM_FOURPDM,true);
    break;
  case (RESTART_NEVPT2PDM):
    Npdm::npdm(NPDM_NEVPT2,true);
    break;
  case (TRANSITION_ONEPDM):
    Npdm::npdm(NPDM_ONEPDM,false,true);
    if (dmrginp.hamiltonian() == BCS) {
      Npdm::npdm(NPDM_PAIRMATRIX,true,true);      
    }
    break;
  case (TRANSITION_TWOPDM):
    Npdm::npdm(NPDM_TWOPDM,false,true);
    break;
  case (TRANSITION_THREEPDM):
    Npdm::npdm(NPDM_THREEPDM,false,true);
    break;
  case (RESTART_T_ONEPDM):
    Npdm::npdm(NPDM_ONEPDM,true,true);
    if (dmrginp.hamiltonian() == BCS) {
      Npdm::npdm(NPDM_PAIRMATRIX,true,true);      
    }
    break;
  case (RESTART_T_TWOPDM):
    Npdm::npdm(NPDM_TWOPDM,true,true);
    break;
  case (RESTART_T_THREEPDM):
    Npdm::npdm(NPDM_THREEPDM,true,true);
    break;
  case(RESTART_NEVPT2):
    nevpt2::nevpt2_restart();
    break;
//EL
   case (DS1_ONEPDM):
     Npdm::npdm(NPDM_DS1,false,true,true);
    break;
   case (RESTART_DS1_ONEPDM):
     Npdm::npdm(NPDM_DS1,true,true,true);
    break;
   case (DS0_ONEPDM):
     Npdm::npdm(NPDM_DS0,false,true,true);
    break;
   case (RESTART_DS0_ONEPDM):
     Npdm::npdm(NPDM_DS0,true,true,true);
    break;
//EL
  default:
    pout << "Invalid calculation types" << endl; abort();
    
  }

  cout.rdbuf(backup);

  tcpu=globaltimer.totalcputime();twall=globaltimer.totalwalltime();
  pout << setprecision(3) <<"\n\n\t\t\t BLOCK CPU  Time (seconds): " << tcpu << endl;
  pout << setprecision(3) <<"\t\t\t BLOCK Wall Time (seconds): " << twall << endl;

  return 0;
}


void calldmrg_(char* input, char* output) {
   int a;
   //a=calldmrg("dmrg.inp",0);//, output);
   a=calldmrg(input,0);//, output);
}


void fullrestartGenblock() {
  SweepParams sweepParams, sweepParamsTmp;
  bool direction; int restartsize;
//Temporary fix to restore sweep direction
//FIXME: NN wrote: please let me know if this makes some erroneous behaviors
  sweepParamsTmp.restorestate(direction, restartsize);
  sweepParams.set_sweep_iter() = 0;
  sweepParams.current_root() = -1;
//direction = true;
  restartsize = 0;

  SweepGenblock::do_one(sweepParams, false, !direction, RESTART, restartsize, -1, -1);
  
  sweepParams.restorestate(direction, restartsize);
  sweepParams.set_sweep_iter()=0;
  sweepParams.set_block_iter() = 0;
  
  sweepParams.savestate(direction, restartsize);
}  


void restart(double sweep_tol, bool reset_iter)
{
  double last_fe = 100.;
  double last_be = 100.;
  double old_fe = 0.;
  double old_be = 0.;
  bool direction;
  int restartsize;
  SweepParams sweepParams;
  bool dodiis = false;

  int domoreIter = 2;

  sweepParams.restorestate(direction, restartsize);

  if (!dmrginp.setStateSpecific()) {
    if(reset_iter) { //this is when you restart from the start of the sweep
      sweepParams.set_sweep_iter() = 0;
      sweepParams.set_restart_iter() = 0;
    }
    
    if (restartwarm)
      last_fe = Sweep::do_one(sweepParams, true, direction, true, restartsize);
    else
      last_fe = Sweep::do_one(sweepParams, false, direction, true, restartsize);
    
    
    while ((fabs(last_fe - old_fe) > sweep_tol) || (fabs(last_be - old_be) > sweep_tol) || 
	   (dmrginp.algorithm_method() == TWODOT_TO_ONEDOT && dmrginp.twodot_to_onedot_iter()+1 >= sweepParams.get_sweep_iter()) )
      {
	
	old_fe = last_fe;
	old_be = last_be;
	if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	  break;
	last_be = Sweep::do_one(sweepParams, false, !direction, false, 0);
	
	
	if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	  break;
	last_fe = Sweep::do_one(sweepParams, false, direction, false, 0);	
      }
  }
  else { //this is state specific calculation  
    const int nroots = dmrginp.nroots();

    bool direction;
    int restartsize;
    sweepParams.restorestate(direction, restartsize);

    //initialize state and canonicalize all wavefunctions
    int currentRoot = sweepParams.current_root();
    for (int i=0; i<nroots; i++) {
      sweepParams.current_root() = i;
      if (mpigetrank()==0) {
	Sweep::InitializeStateInfo(sweepParams, direction, i);
	Sweep::InitializeStateInfo(sweepParams, !direction, i);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, i);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, i);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, i);
      }
    }

    //now generate overlaps with all the previous wavefunctions
    for (int i=0; i<currentRoot; i++) {
      sweepParams.current_root() = i;
      if (mpigetrank()==0) {
	for (int j=0; j<i; j++) {
	  int integralIndex = 0;
	  Sweep::InitializeOverlapSpinBlocks(sweepParams, !direction, i, j, integralIndex);
	  Sweep::InitializeOverlapSpinBlocks(sweepParams, direction, i, j, integralIndex);
	}
      }
    }
    sweepParams.current_root() = currentRoot;

    if (sweepParams.current_root() <0) {
      p1out << "This is most likely not a restart calculation and should be done without the restart command!!"<<endl;
      p1out << "Aborting!!"<<endl;
      exit(0);
    }
    pout << "RESTARTING STATE SPECIFIC CALCULATION OF STATE "<<sweepParams.current_root()<<" AT SWEEP ITERATION  "<<sweepParams.get_sweep_iter()<<endl;

    //this is so that the iteration is not one ahead after generate block for restart
    --sweepParams.set_sweep_iter(); sweepParams.savestate(direction, restartsize);
    for (int i=sweepParams.current_root(); i<nroots; i++) {
      sweepParams.current_root() = i;

      p1out << "RUNNING GENERATE BLOCKS FOR STATE "<<i<<endl;

      if (mpigetrank()==0) {
	Sweep::InitializeStateInfo(sweepParams, direction, i);
	Sweep::InitializeStateInfo(sweepParams, !direction, i);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, i);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, i);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, i);
	for (int j=0; j<i ; j++) {
	  int integralIndex = 0;
	  Sweep::InitializeOverlapSpinBlocks(sweepParams, direction, i, j, integralIndex);
	  Sweep::InitializeOverlapSpinBlocks(sweepParams, !direction, i, j, integralIndex);
	}
      }
      SweepGenblock::do_one(sweepParams, false, !direction, false, 0, i, i);
      
      
      p1out << "STATE SPECIFIC CALCULATION FOR STATE: "<<i<<endl;
      dmrg_stateSpecific(sweep_tol, i);
      p1out << "STATE SPECIFIC CALCULATION FOR STATE: "<<i<<" FINSIHED"<<endl;

      sweepParams.set_sweep_iter() = 0;
      sweepParams.set_restart_iter() = 0;
      sweepParams.savestate(!direction, restartsize);
    }

    p1out << "ALL STATE SPECIFIC CALCUALTIONS FINISHED"<<endl;
  }


  if(dmrginp.max_iter() <= sweepParams.get_sweep_iter()){
    pout << "\n\t\t\t Maximum sweep iterations achieved " << std::endl;
  }

}

void dmrg(double sweep_tol)
{
  double last_fe = 10.e6;
  double last_be = 10.e6;
  double old_fe = 0.;
  double old_be = 0.;
  SweepParams sweepParams;

  int old_states=sweepParams.get_keep_states();
  int new_states;
  double old_error=0.0;
  double old_energy=0.0;
  // warm up sweep ...
  bool dodiis = false;

  int domoreIter = 0;
  bool direction;

  //this is regular dmrg calculation
  if(!dmrginp.setStateSpecific()) {
    sweepParams.current_root() = -1;
    last_fe = Sweep::do_one(sweepParams, true, true, false, 0);
    direction = false;
    while ((fabs(last_fe - old_fe) > sweep_tol) || (fabs(last_be - old_be) > sweep_tol) || 
	   (dmrginp.algorithm_method() == TWODOT_TO_ONEDOT && dmrginp.twodot_to_onedot_iter()+1 >= sweepParams.get_sweep_iter()) )
    {
      old_fe = last_fe;
      old_be = last_be;
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	break;
      last_be = Sweep::do_one(sweepParams, false, false, false, 0);
      direction = true;
      pout << "\t\t\t Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;
      
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	break;
      
      //For obtaining the extrapolated energy
      old_states=sweepParams.get_keep_states();
      new_states=sweepParams.get_keep_states_ls();
      
      last_fe = Sweep::do_one(sweepParams, false, true, false, 0);
      direction = false;
      
      new_states=sweepParams.get_keep_states();
      
      
      pout << "\t\t\t Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;
      if (domoreIter == 2) {
	dodiis = true;
	break;
      }
      
    }
  }
  else { //this is state specific calculation  
    const int nroots = dmrginp.nroots();

    bool direction=true;
    int restartsize;
    //sweepParams.restorestate(direction, restartsize);
    //sweepParams.set_sweep_iter() = 0;
    //sweepParams.set_restart_iter() = 0;

    algorithmTypes atype;
    pout << "STARTING STATE SPECIFIC CALCULATION "<<endl;
    for (int i=0; i<nroots; i++) {
      atype = dmrginp.algorithm_method();
      dmrginp.set_algorithm_method() = ONEDOT;
      sweepParams.current_root() = i;

      p1out << "RUNNING GENERATE BLOCKS FOR STATE "<<i<<endl;

      if (mpigetrank()==0) {
	Sweep::InitializeStateInfo(sweepParams, direction, i);
	Sweep::InitializeStateInfo(sweepParams, !direction, i);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, i);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, i);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, i);
	Sweep::InitializeStateInfo(sweepParams, direction, i);
	Sweep::InitializeStateInfo(sweepParams, !direction, i);

      }

      for (int j=0; j<i ; j++) {
	int integralIndex = 0;
	Sweep::InitializeOverlapSpinBlocks(sweepParams, direction, i, j, integralIndex);
	Sweep::InitializeOverlapSpinBlocks(sweepParams, !direction, i, j, integralIndex);
      }
      dmrginp.set_algorithm_method() = atype;

      p1out << "RUNNING GENERATE BLOCKS FOR STATE "<<i<<endl;

      SweepGenblock::do_one(sweepParams, false, !direction, false, 0, i, i);
      sweepParams.set_sweep_iter() = 0;
      sweepParams.set_restart_iter() = 0;
      sweepParams.savestate(!direction, restartsize);

      
      pout << "STATE SPECIFIC CALCULATION FOR STATE: "<<i<<endl;
      dmrg_stateSpecific(sweep_tol, i);
      pout << "STATE SPECIFIC CALCULATION FOR STATE: "<<i<<" FINSIHED"<<endl;
    }

    pout << "ALL STATE SPECIFIC CALCUALTIONS FINISHED"<<endl;
  }
}

void responseSweep(double sweep_tol, int targetState, vector<int>& projectors, vector<int>& baseStates)
{
  double last_fe = 1.e6;
  double last_be = 1.e6;
  double old_fe = 0.;
  double old_be = 0.;
  SweepParams sweepParams;

  bool direction, warmUp, restart;
  int restartSize=0;
  direction = true; //forward
  warmUp = true; //startup sweep
  restart = false; //not a restart

  sweepParams.current_root() = -1;

  algorithmTypes atype = dmrginp.algorithm_method();
  dmrginp.set_algorithm_method() = ONEDOT;

  //the baseState is the initial guess for the targetState
  if (FULLRESTART) {
    sweepParams.restorestate(direction, restartSize);
    direction = !direction;
    dmrginp.setGuessState() = targetState;
    last_fe = SweepResponse::do_one(sweepParams, warmUp, direction, restart, restartSize, targetState, projectors, baseStates);
    bool tempdirection;
    sweepParams.restorestate(tempdirection, restartSize);
    sweepParams.calc_niter();
    sweepParams.set_sweep_iter() = 0;
    sweepParams.set_restart_iter() = 0;
    sweepParams.savestate(tempdirection, restartSize);
  }
  else if (RESTART) {
    dmrginp.set_algorithm_method() = atype;
    warmUp = false;
    restart = true;
    sweepParams.restorestate(direction, restartSize);
    last_fe = SweepResponse::do_one(sweepParams, warmUp, direction, restart, restartSize, targetState, projectors, baseStates);
  }
  else 
    last_fe = SweepResponse::do_one(sweepParams, warmUp, direction, restart, restartSize, targetState, projectors, baseStates);

  dmrginp.set_algorithm_method() = atype;
  restart = false;
  restartSize = 0;
  warmUp = false;
  while ( true)
    {
      old_fe = last_fe;
      old_be = last_be;
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	break;

      last_be = SweepResponse::do_one(sweepParams, warmUp, !direction, restart, restartSize, targetState, projectors, baseStates);
      p1out << "\t\t\t Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;
      
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	break;
      
      last_fe = SweepResponse::do_one(sweepParams, warmUp, direction, restart, restartSize, targetState, projectors, baseStates);

      
      pout << "\t\t\t Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;
      
    }
  
}

/*
void restartResponseSweep(double sweep_tol, int targetState, int correctionVector, int baseState)
{
  double last_fe = 10.e6;
  double last_be = 10.e6;
  double old_fe = 0.;
  double old_be = 0.;
  SweepParams sweepParams;
  bool direction, warmUp=false, restart=true;
  int restartSize=0;

  sweepParams.restorestate(direction, restartSize);

  sweepParams.current_root() = -1;

  last_fe = SweepResponse::do_one(sweepParams, warmUp, direction, restart, restartSize, targetState, correctionVector, baseState);

  warmUp = false;
  restart = false;
  while ( true)
    {
      old_fe = last_fe;
      old_be = last_be;
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	break;
      last_be = SweepResponse::do_one(sweepParams, warmUp, !direction, restart, restartSize, targetState, correctionVector, baseState);
      p1out << "\t\t\t Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;
      
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	break;
      
      direction = true;
      last_fe = SweepResponse::do_one(sweepParams, warmUp, direction, restart, restartSize, targetState, correctionVector, baseState);

      
      p1out << "\t\t\t Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;
      
    }
  
}
*/


void compress(double sweep_tol, int targetState, int baseState)
{
  double last_fe = 10.e6;
  double last_be = 10.e6;
  double old_fe = 0.;
  double old_be = 0.;
  SweepParams sweepParams;
  bool direction;

  sweepParams.current_root() = -1;
  //this is the warmup sweep, the baseState is used as the initial guess for the targetState
  last_fe = SweepCompress::do_one(sweepParams, true, true, false, 0, targetState, baseState);

  direction = false;
  while ( true)
    {
      old_fe = last_fe;
      old_be = last_be;
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	break;
      last_be = SweepCompress::do_one(sweepParams, false, false, false, 0, targetState, baseState);
      direction = true;
      pout << "\t\t\t Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;
      
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	break;
      
      last_fe = SweepCompress::do_one(sweepParams, false, true, false, 0, targetState, baseState);
      direction = false;
      
      
      pout << "\t\t\t Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;
      
    }

  //we finally canonicalize the targetState
  //one has to canonicalize the wavefunction with atleast 3 sweeps, this is a quirk of the way 
  //we transform wavefunction
  if (mpigetrank()==0) {
    Sweep::InitializeStateInfo(sweepParams, !direction, targetState);
    Sweep::InitializeStateInfo(sweepParams, direction, targetState);
    Sweep::CanonicalizeWavefunction(sweepParams, !direction, targetState);
    Sweep::CanonicalizeWavefunction(sweepParams, direction, targetState);
    Sweep::CanonicalizeWavefunction(sweepParams, !direction, targetState);
    Sweep::InitializeStateInfo(sweepParams, !direction, targetState);
    Sweep::InitializeStateInfo(sweepParams, direction, targetState);
    
  }
  
}

void dmrg_stateSpecific(double sweep_tol, int targetState)
{
  double last_fe = 10.e6;
  double last_be = 10.e6;
  double old_fe = 0.;
  double old_be = 0.;
  int ls_count=0;
  SweepParams sweepParams;
  int old_states=sweepParams.get_keep_states();
  int new_states;
  double old_error=0.0;
  double old_energy=0.0;
  // warm up sweep ...

  bool direction;
  int restartsize;
  sweepParams.restorestate(direction, restartsize);

  //initialize array of size m_maxiter or dmrginp.max_iter() for dw and energy
  sweepParams.current_root() = targetState;

  last_fe = Sweep::do_one(sweepParams, false, direction, true, restartsize);

  while ((fabs(last_fe - old_fe) > sweep_tol) || (fabs(last_be - old_be) > sweep_tol)  )
    {
      old_fe = last_fe;
      old_be = last_be;
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter()) 
	break;

      last_be = Sweep::do_one(sweepParams, false, !direction, false, 0);
      pout << "\t\t\t Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;

      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	break;


      last_fe = Sweep::do_one(sweepParams, false, direction, false, 0);

      new_states=sweepParams.get_keep_states();


      pout << "\t\t\t Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;

    }
  pout << "Converged Energy  " << sweepParams.get_lowest_energy()[0]<< std::endl;
  if(dmrginp.max_iter() <= sweepParams.get_sweep_iter()) {
    
    pout << "Maximum sweep iterations achieved " << std::endl;
  }

  //one has to canonicalize the wavefunction with atleast 3 sweeps, this is a quirk of the way 
  //we transform wavefunction
  if (mpigetrank()==0) {
    Sweep::InitializeStateInfo(sweepParams, !direction, targetState);
    Sweep::InitializeStateInfo(sweepParams, direction, targetState);
    Sweep::CanonicalizeWavefunction(sweepParams, !direction, targetState);
    Sweep::CanonicalizeWavefunction(sweepParams, direction, targetState);
    Sweep::CanonicalizeWavefunction(sweepParams, !direction, targetState);
    Sweep::InitializeStateInfo(sweepParams, !direction, targetState);
    Sweep::InitializeStateInfo(sweepParams, direction, targetState);
    
  }

}


