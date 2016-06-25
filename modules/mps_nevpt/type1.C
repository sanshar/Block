#include "guess_wavefunction.h"
#include "global.h"
#include "solver.h"
#include "initblocks.h"
#include "MatrixBLAS.h"
#include <boost/format.hpp>
#ifndef SERIAL
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include "rotationmat.h"
#include "density.h"
#include "pario.h"
#include "davidson.h"
#include "sweep.h"
//#include "mps_mpi.h"
#include "type1.h"
#include "operatorfunctions.h"
#include "StateInfo.h"
#include <stdio.h>

using namespace boost;
using namespace std;


double SpinAdapted::mps_nevpt::type1::do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize, perturber& pb, int baseState)
{
  int integralIndex = 0;
  SpinBlock system;
  system.nonactive_orb() = pb.orb();
  const int nroots = dmrginp.nroots(sweepParams.get_sweep_iter());

  std::vector<double> finalEnergy(nroots,-1.0e10);
  std::vector<double> finalEnergy_spins(nroots,0.);
  double finalError = 0.;

  sweepParams.set_sweep_parameters();
  // a new renormalisation sweep routine
  if (forward)
    {
      if (dmrginp.outputlevel() > 0)
	{ pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in forwards direction"<<endl; }
    }
  else
    {
      if (dmrginp.outputlevel() > 0)
	{
	  pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in backwards direction" << endl;
	  pout << "\t\t\t ============================================================================ " << endl;
	}
    }

  InitBlocks::InitStartingBlock (system,forward, baseState, pb.wavenumber(), sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), restartSize, restart, warmUp, integralIndex, pb.braquanta, pb.ketquanta);
  if(!restart)
    sweepParams.set_block_iter() = 0;

 
  if (dmrginp.outputlevel() > 0)
    pout << "\t\t\t Starting block is :: " << endl << system << endl;

  SpinBlock::store (forward, system.get_sites(), system, pb.wavenumber(), baseState); // if restart, just restoring an existing block --
  sweepParams.savestate(forward, system.get_sites().size());
  bool dot_with_sys = true;
  vector<int> syssites = system.get_sites();

  if (restart)
  {
    if (forward && system.get_complementary_sites()[0] >= dmrginp.last_site()/2)
      dot_with_sys = false;
    if (!forward && system.get_sites()[0]-1 < dmrginp.last_site()/2)
      dot_with_sys = false;
  }
  if (dmrginp.outputlevel() > 0)
    mcheck("at the very start of sweep");  // just timer

  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); ) // get_n_iters() returns the number of blocking iterations needed in one sweep
    {
      if (dmrginp.outputlevel() > 0)
      {
        pout << "\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
        pout << "\t\t\t ----------------------------" << endl;
      }
      if (dmrginp.outputlevel() > 0) {
	    if (forward)
      {
	      pout << "\t\t\t Current direction is :: Forwards " << endl;
      }
	    else
      {
	      pout << "\t\t\t Current direction is :: Backwards " << endl;
      }
      }

      if (sweepParams.get_block_iter() != 0) 
	sweepParams.set_guesstype() = TRANSFORM;
      else
        sweepParams.set_guesstype() = TRANSPOSE;


      
      if (dmrginp.outputlevel() > 0)
         pout << "\t\t\t Blocking and Decimating " << endl;
	  
      SpinBlock newSystem; // new system after blocking and decimating
      newSystem.nonactive_orb() = pb.orb();

      //Need to substitute by:
     // if (warmUp )
     //   Startup(sweepParams, system, newSystem, dot_with_sys, pb.wavenumber(), baseState);
     // else {
     //   BlockDecimateAndCompress (sweepParams, system, newSystem, false, dot_with_sys, pb.wavenumber(), baseState);
     // }
      
        BlockDecimateAndCompress (sweepParams, system, newSystem, warmUp, dot_with_sys,pb, baseState);
      //Need to substitute by?


      system = newSystem;
      if (dmrginp.outputlevel() > 0){
	    pout << system<<endl;
	    pout << system.get_braStateInfo()<<endl;
	    system.printOperatorSummary();
      }
      
      //system size is going to be less than environment size
      if (forward && system.get_complementary_sites()[0] >= dmrginp.last_site()/2)
	    dot_with_sys = false;
      if (!forward && system.get_sites()[0]-1 < dmrginp.last_site()/2)
	    dot_with_sys = false;

      SpinBlock::store (forward, system.get_sites(), system, pb.wavenumber(), baseState);	 	
      syssites = system.get_sites();
      if (dmrginp.outputlevel() > 0)
	      pout << "\t\t\t saving state " << syssites.size() << endl;
      ++sweepParams.set_block_iter();
      
#ifndef SERIAL
      mpi::communicator world;
      world.barrier();
#endif
      sweepParams.savestate(forward, syssites.size());
      if (dmrginp.outputlevel() > 0)
         mcheck("at the end of sweep iteration");
    }

  //FIXME
  //It does not seem necessary.

  //when we are doing twodot, we still need to do the last sweep to make sure that the
  //correctionVector and base wavefunction are propogated correctly across sweeps
//  //especially when we switch from twodot to onedot algorithm
//  if (!sweepParams.get_onedot() && !warmUp) {
//      pout << "\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
//      pout << "\t\t\t ----------------------------" << endl;
//      if (dmrginp.outputlevel() > 0) {
//	    if (forward)
//	      pout << "\t\t\t Current direction is :: Forwards " << endl;
//	    else
//	      pout << "\t\t\t Current direction is :: Backwards " << endl;
//      }
//    sweepParams.set_onedot() = true;
//    sweepParams.set_env_add() = 0;
//    bool dot_with_sys = true;
//    WavefunctionCanonicalize(sweepParams, system, warmUp, dot_with_sys, targetState, baseState);
//    sweepParams.set_onedot() = false;
//    sweepParams.set_env_add() = 1;
//  }
//

  pout << "\t\t\t Largest Error for Sweep with " << sweepParams.get_keep_states() << " states is " << finalError << endl;
  pout << "\t\t\t Largest overlap for Sweep with " << sweepParams.get_keep_states() << " states is " << finalEnergy[0] << endl;
  sweepParams.set_largest_dw() = finalError;
  

  pout << "\t\t\t ============================================================================ " << endl;

  // update the static number of iterations

  ++sweepParams.set_sweep_iter();

  return finalError;
}

void SpinAdapted::mps_nevpt::type1::BlockDecimateAndCompress (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys, perturber& pb, int baseState)
{
  int sweepiter = sweepParams.get_sweep_iter();
  if (dmrginp.outputlevel() > 0) {
    mcheck("at the start of block and decimate");
    pout << "\t\t\t dot with system "<<dot_with_sys<<endl;
    pout <<endl<< "\t\t\t Performing Blocking"<<endl;
  }
  // figure out if we are going forward or backwards
  dmrginp.guessgenT -> start();
  bool forward = (system.get_sites() [0] == 0);
  SpinBlock systemDot;
  SpinBlock environment, environmentDot, newEnvironment;
  SpinBlock big;
  environment.nonactive_orb() = pb.orb();
  newEnvironment.nonactive_orb() = pb.orb();
  int systemDotStart, systemDotEnd;
  int environmentDotStart, environmentDotEnd, environmentStart, environmentEnd;
  int systemDotSize = sweepParams.get_sys_add() - 1;
  int environmentDotSize = sweepParams.get_env_add() -1;
  if (forward)
  {
    systemDotStart = dmrginp.spinAdapted() ? *system.get_sites().rbegin () + 1 : (*system.get_sites().rbegin ())/2 + 1 ;
    systemDotEnd = systemDotStart + systemDotSize;
    environmentDotStart = systemDotEnd + 1;
    environmentDotEnd = environmentDotStart + environmentDotSize;
  }
  else
  {
    systemDotStart = dmrginp.spinAdapted() ? system.get_sites()[0] - 1 : (system.get_sites()[0])/2 - 1 ;
    systemDotEnd = systemDotStart - systemDotSize;
    environmentDotStart = systemDotEnd - 1;
    environmentDotEnd = environmentDotStart - environmentDotSize;
  }
  systemDot = SpinBlock(systemDotStart, systemDotEnd, pb.orb());
  environmentDot = SpinBlock(environmentDotStart, environmentDotEnd, pb.orb());

  Sweep::makeSystemEnvironmentBigBlocks(system, systemDot, newSystem, environment, environmentDot, newEnvironment, big, sweepParams, dot_with_sys, useSlater, system.get_integralIndex(), pb.wavenumber(), baseState,pb.braquanta,pb.ketquanta);


  //analyse_operator_distribution(big);
  dmrginp.guessgenT -> stop();
  dmrginp.multiplierT -> start();
  std::vector<Matrix> rotatematrix;

  if (dmrginp.outputlevel() > 0)
    mcheck(""); 
  if (dmrginp.outputlevel() > 0) {
    if (!dot_with_sys && sweepParams.get_onedot())  { pout << "\t\t\t System  Block"<<system;    }
    else pout << "\t\t\t System  Block"<<newSystem;
    pout << "\t\t\t Environment Block"<<newEnvironment<<endl;
    pout << "\t\t\t Solving wavefunction "<<endl;
  }

  std::vector<Wavefunction> solution; solution.resize(1);
  std::vector<Wavefunction> outputState; outputState.resize(1);

  DiagonalMatrix e;


  //read the 0th wavefunction which we keep on the ket side because by default the ket stateinfo is used to initialize wavefunction
  //also when you use spinblock operators to multiply a state, it does so from the ket side i.e.  H|ket>
  //GuessWave::guess_wavefunctions(solution, e, big, sweepParams.set_guesstype(), sweepParams.get_onedot(), dot_with_sys, 0.0, baseState); 
  GuessWave::guess_wavefunctions(solution[0], e, big, sweepParams.set_guesstype(), sweepParams.get_onedot(), baseState, dot_with_sys, 0.0); 

#ifndef SERIAL
  mpi::communicator world;
  broadcast(world, solution, 0);
#endif
  
  outputState[0].AllowQuantaFor(big.get_leftBlock()->get_braStateInfo(), big.get_rightBlock()->get_braStateInfo(),pb.braquanta);
  outputState[0].set_onedot(sweepParams.get_onedot());
  outputState[0].Clear();
  if (pb.type() == TwoPerturbType::Va)
    big.multiplyCDD_sum(solution[0],&(outputState[0]),MAX_THRD);
  if (pb.type() == TwoPerturbType::Vi)
    big.multiplyCCD_sum(solution[0],&(outputState[0]),MAX_THRD);

  //davidson_f(solution[0], outputState[0]);
  SpinBlock newbig;

  if (sweepParams.get_onedot() && !dot_with_sys)
  {
    InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, baseState, pb.wavenumber(), systemDot.size(), dmrginp.direct(), system.get_integralIndex(), DISTRIBUTED_STORAGE, false, true,NO_PARTICLE_SPIN_NUMBER_CONSTRAINT,pb.braquanta,pb.ketquanta);
    InitBlocks::InitBigBlock(newSystem, environment, newbig,pb.braquanta,pb.ketquanta); 

    Wavefunction tempwave = outputState[0];
    GuessWave::onedot_shufflesysdot(big.get_braStateInfo(), newbig.get_braStateInfo(), outputState[0], tempwave);  
    outputState[0] = tempwave;

    tempwave = solution[0];
    GuessWave::onedot_shufflesysdot(big.get_ketStateInfo(), newbig.get_ketStateInfo(), solution[0], tempwave);  
    solution[0] = tempwave;

    big.get_rightBlock()->clear();
    big.clear();
  }
  else
    newbig = big;
  
  DensityMatrix bratracedMatrix(newSystem.get_braStateInfo());
  bratracedMatrix.allocate(newSystem.get_braStateInfo());

  //bratracedMatrix.makedensitymatrix(outputState, newbig, dmrginp.weights(sweepiter), 0.0, 0.0, true);
  bratracedMatrix.makedensitymatrix(outputState, newbig, std::vector<double>(1,1.0), 0.0, 0.0, true);
  if (sweepParams.get_noise() > NUMERICAL_ZERO) {
    pout << "adding noise  "<<trace(bratracedMatrix)<<"  "<<sweepiter<<"  "<<dmrginp.weights(sweepiter)[0]<<endl;
    bratracedMatrix.add_onedot_noise_forCompression(solution[0], newbig, sweepParams.get_noise()*max(1.0,trace(bratracedMatrix)));
    if (trace(bratracedMatrix) <1e-14) 
      bratracedMatrix.SymmetricRandomise();
      
    pout << "after noise  "<<trace(bratracedMatrix)<<"  "<<sweepParams.get_noise()<<endl;
  }
  environment.clear();
  newEnvironment.clear();


  std::vector<Matrix> brarotateMatrix, ketrotateMatrix;
  LoadRotationMatrix (newSystem.get_sites(), ketrotateMatrix, baseState);

  double braerror;
  if (!mpigetrank()) {
    braerror = makeRotateMatrix(bratracedMatrix, brarotateMatrix, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());
  }

#ifndef SERIAL
  broadcast(world, ketrotateMatrix, 0);
  broadcast(world, brarotateMatrix, 0);
#endif

  if (dmrginp.outputlevel() > 0)
    pout << "\t\t\t Total bra discarded weight "<<braerror<<endl<<endl;

  sweepParams.set_lowest_error() = braerror;

  SaveRotationMatrix (newbig.get_leftBlock()->get_sites(), brarotateMatrix, pb.wavenumber());
  //FIXME
  //It is neccessary for twodot algorithm to save baseState wavefuntion.
  //I do not know why. 
  solution[0].SaveWavefunctionInfo (newbig.get_ketStateInfo(), newbig.get_leftBlock()->get_sites(), baseState);
  outputState[0].SaveWavefunctionInfo (newbig.get_braStateInfo(), newbig.get_leftBlock()->get_sites(), pb.wavenumber());
  //TODO 
  //Why do I need this?
  //They should have been consistent.
//  solution[0].SaveWavefunctionInfo (newbig.get_ketStateInfo(), newbig.get_leftBlock()->get_sites(), baseState);
//  SaveRotationMatrix (newbig.get_leftBlock()->get_sites(), ketrotateMatrix, baseState);

  if (dmrginp.outputlevel() > 0)
    pout <<"\t\t\t Performing Renormalization "<<endl;
  newSystem.transform_operators(brarotateMatrix, ketrotateMatrix);

  if (dmrginp.outputlevel() > 0)
    mcheck("after rotation and transformation of block");

  if (dmrginp.outputlevel() > 0){
    pout << *dmrginp.guessgenT<<" "<<*dmrginp.multiplierT<<" "<<*dmrginp.operrotT<< "  "<<globaltimer.totalwalltime()<<" timer "<<endl;
    pout << *dmrginp.makeopsT<<" makeops "<<endl;
    pout << *dmrginp.datatransfer<<" datatransfer "<<endl;
    pout <<"oneindexopmult   twoindexopmult   Hc  couplingcoeff"<<endl;  
    pout << *dmrginp.oneelecT<<" "<<*dmrginp.twoelecT<<" "<<*dmrginp.hmultiply<<" "<<*dmrginp.couplingcoeff<<" hmult"<<endl;
    pout << *dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build block"<<endl;
    pout << "addnoise  S_0_opxop  S_1_opxop   S_2_opxop"<<endl;
    pout << *dmrginp.addnoise<<" "<<*dmrginp.s0time<<" "<<*dmrginp.s1time<<" "<<*dmrginp.s2time<<endl;
  }

}

void SpinAdapted::mps_nevpt::type1::cleanup(int baseState, const perturber& pb, int cleanlevel)
{
  for (int site=0; site < dmrginp.last_site(); site++)
  {
    std::string file = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%s") % dmrginp.save_prefix() % "/SpinBlock-forward-"% 0 % "-" % site % "." %pb.wavenumber() % "." % baseState % "." %0 % "." % mpigetrank() % ".tmp" );
    if (boost::filesystem::exists(file)) remove(file.c_str());
    file = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%s") % dmrginp.save_prefix() % "/SpinBlock-backward-"% 0 % "-" % site % "." % pb.wavenumber() % "." % baseState % "." %0 % "." % mpigetrank() % ".tmp" );
    if (boost::filesystem::exists(file)) remove(file.c_str());
    file = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%s") % dmrginp.save_prefix() % "/SpinBlock-forward-"% (site+1) % "-" % (dmrginp.last_site()-1) % "." % pb.wavenumber() % "." % baseState % "." %0 % "." % mpigetrank() % ".tmp" );
    if (boost::filesystem::exists(file)) remove(file.c_str());
    file = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%s") % dmrginp.save_prefix() % "/SpinBlock-backward-"% (site+1) % "-" % (dmrginp.last_site()-1) % "." % pb.wavenumber() % "." % baseState % "." %0 % "." % mpigetrank() % ".tmp" );
    if (boost::filesystem::exists(file)) remove(file.c_str());

  }
}


void SpinAdapted::mps_nevpt::type1::subspace_Vi(int baseState)
{
  
  double energy=0;
  double overlap=0;
  ViPerturber pb;
  MPS::siteBlocks.clear();
  int coresize = dmrginp.spinAdapted()? dmrginp.core_size():dmrginp.core_size()*2;
  int coreshift = dmrginp.spinAdapted()? dmrginp.act_size(): dmrginp.act_size()*2;
  //int virtsize = dmrginp.virt_size();
  //int virtshift = dmrginp.core_size()+dmrginp.act_size();
  for(int i=0; i< coresize; i++){
    double perturberEnergy=0;
    dmrginp.calc_type() = MPS_NEVPT;
    pb.init(i+coreshift);
    pout << "Begin Vi subspace with i = " << pb.orb(0)<<endl;
    SweepParams sweepParams;
    sweepParams.set_sweep_parameters();
    sweepParams.current_root() = baseState;
    //sweepParams.current_root() = -1;
    //double last_fe = Startup(sweepParams, true, true, false, 0, pb, baseState);
    Timer timer;
    Startup(sweepParams, true, pb, baseState);
    pout <<"Start up time :" << timer.elapsedwalltime();
    //sweepParams.current_root() = baseState;
    timer.start();
    while(true)
    {
      do_one(sweepParams, false, false, false, 0, pb, baseState);
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	      break;
      do_one(sweepParams, false, true, false, 0, pb, baseState);
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
      {
        cleanup(baseState, pb);
	      break;
      }
    }
    pout <<"Sweep time :" << timer.elapsedwalltime();
//
//    if (mpigetrank()==0) {
//      bool direction = true;
//      Sweep::InitializeStateInfo(sweepParams, !direction, pb.wavenumber(),pb.braquanta );
//      Sweep::InitializeStateInfo(sweepParams, direction, pb.wavenumber(),pb.braquanta );
//      Sweep::CanonicalizeWavefunction(sweepParams, !direction, pb.wavenumber(),pb.braquanta );
//      Sweep::CanonicalizeWavefunction(sweepParams, direction, pb.wavenumber(),pb.braquanta );
//      Sweep::CanonicalizeWavefunction(sweepParams, !direction, pb.wavenumber(),pb.braquanta );
//      Sweep::InitializeStateInfo(sweepParams, !direction, pb.wavenumber(),pb.braquanta );
//      Sweep::InitializeStateInfo(sweepParams, direction, pb.wavenumber(),pb.braquanta );
//      
//    }
//

    MPS pbmps(pb.wavenumber());
    double o, h;
    dmrginp.calc_type() = DMRG;

    timer.start();
    calcHamiltonianAndOverlap(pbmps, h, o,pb);

    pout <<"Calculate Expectation time :" << timer.elapsedwalltime();



    if(!dmrginp.spinAdapted())
    {
      //In nonspinAdapted, alpha and beta have the results. Only one is neccessary. 
      o*=2;
      h*=2;
      i++;
    }
    if(o> NUMERICAL_ZERO){
      double fock =dmrginp.spinAdapted()? v_1[0](2*(i+coreshift),2*(i+coreshift)): v_1[0](i+coreshift,i+coreshift);
      //perturberEnergy = h/o+fock+perturber::CoreEnergy[0];
      perturberEnergy = h/o - fock;
      energy += o/(mps_nevpt::ZeroEnergy[baseState]- perturberEnergy) ;
      //overlap +=o;
      overlap += sqrt(o)/(mps_nevpt::ZeroEnergy[baseState]- perturberEnergy);
      if (dmrginp.outputlevel() > 0)
      {
        pout << "Zero Energy: " << mps_nevpt::ZeroEnergy[baseState]<<endl;
        pout << "Amplitude : " << sqrt(o)/(mps_nevpt::ZeroEnergy[baseState]- perturberEnergy) <<endl;
        pout << "Ener(only CAS part) : " << h/o<<endl;
        pout << "Energy : " << perturberEnergy<<endl;
        pout << "Correction Energy: "<< o/(mps_nevpt::ZeroEnergy[baseState]- perturberEnergy)<<endl; 
      }
    }
    else{
      if (dmrginp.outputlevel() > 0){
        pout << "Amplitude : " << 0.0 <<endl;
        pout << "Energy : " << 0.0<<endl;
      }
    }



  }
  pout << "Nevpt2 correction to the energy for state 0 in subspace Vi is " << energy<<endl;;
  pout << "Nevpt2 Vi subspace perturber Amplitude : " << overlap<<endl;;
  //pout << "Core Energy of nevpt2 " <<perturber::CoreEnergy[0]<<endl;
  std::string file = str(boost::format("%s%s%d") % dmrginp.load_prefix() % "/Vi_" % baseState);
  std::fstream f(file,std::fstream::out);
  f << energy <<endl;
  f << overlap <<endl;
  f.close();
}

void SpinAdapted::mps_nevpt::type1::subspace_Va(int baseState)
{
  
  double energy=0;
  double overlap=0;
  VaPerturber pb;
  MPS::siteBlocks.clear();
  int virtsize = dmrginp.spinAdapted()? dmrginp.virt_size():dmrginp.virt_size()*2;
  int virtshift = dmrginp.spinAdapted()? dmrginp.core_size()+dmrginp.act_size(): (dmrginp.core_size()+dmrginp.act_size())*2;
  //int virtsize = dmrginp.virt_size();
  //int virtshift = dmrginp.core_size()+dmrginp.act_size();
  for(int i=0; i< virtsize; i++){
    double perturberEnergy=0;
    dmrginp.calc_type() = MPS_NEVPT;
    pb.init(i+virtshift);
    pout << "Begin Va subspace with a = " << pb.orb(0)<<endl;
    SweepParams sweepParams;
    sweepParams.set_sweep_parameters();
    sweepParams.current_root() = baseState;
    //sweepParams.current_root() = -1;
    //double last_fe = Startup(sweepParams, true, true, false, 0, pb, baseState);
    Timer timer;
    Startup(sweepParams, true, pb, baseState);
    pout <<"Start up time :" << timer.elapsedwalltime();
    //sweepParams.current_root() = baseState;
    timer.start();
    while(true)
    {
      do_one(sweepParams, false, false, false, 0, pb, baseState);
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	      break;
      do_one(sweepParams, false, true, false, 0, pb, baseState);
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
      {
        cleanup(baseState, pb);
	      break;
      }
    }
    pout <<"Sweep time :" << timer.elapsedwalltime();
//    while ( true)
//      {
//        old_fe = last_fe;
//        old_be = last_be;
//        if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
//	        break;
//        last_be = do_one(sweepParams, false, false, false, 0, pb, baseState);
//        if (dmrginp.outputlevel() > 0) 
//	        pout << "Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;
//        
//        if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
//	        break;
//        
//        last_fe = do_one(sweepParams, false, true, false, 0, pb, baseState);
//        
//        if (dmrginp.outputlevel() > 0)
//	        pout << "Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;
//        
//      }
//

 //   if (mpigetrank()==0) {
 //     bool direction = true;
 //     Sweep::InitializeStateInfo(sweepParams, !direction, pb.wavenumber(),pb.braquanta );
 //     Sweep::InitializeStateInfo(sweepParams, direction, pb.wavenumber(),pb.braquanta );
 //     Sweep::CanonicalizeWavefunction(sweepParams, !direction, pb.wavenumber(),pb.braquanta );
 //     Sweep::CanonicalizeWavefunction(sweepParams, direction, pb.wavenumber(),pb.braquanta );
 //     Sweep::CanonicalizeWavefunction(sweepParams, !direction, pb.wavenumber(),pb.braquanta );
 //     Sweep::InitializeStateInfo(sweepParams, !direction, pb.wavenumber(),pb.braquanta );
 //     Sweep::InitializeStateInfo(sweepParams, direction, pb.wavenumber(),pb.braquanta );
 //     
 //   }
   
    MPS pbmps(pb.wavenumber());
    double o, h;
    dmrginp.calc_type() = DMRG;

    timer.start();

    calcHamiltonianAndOverlap(pbmps, h, o,pb);

    pout <<"Calculate Expectation time :" << timer.elapsedwalltime();

    if(!dmrginp.spinAdapted())
    {
      //In nonspinAdapted, alpha and beta have the results. Only one is neccessary. 
      o*=2;
      h*=2;
      i++;
    }
    if(o> NUMERICAL_ZERO){
      double fock =dmrginp.spinAdapted()? v_1[0](2*(i+virtshift),2*(i+virtshift)): v_1[0](i+virtshift,i+virtshift);
      //perturberEnergy = h/o+fock+perturber::CoreEnergy[0];
      perturberEnergy = h/o+fock;
      energy += o/(mps_nevpt::ZeroEnergy[baseState]- perturberEnergy) ;
      //overlap +=o;
      overlap += sqrt(o)/(mps_nevpt::ZeroEnergy[baseState]- perturberEnergy);
      if (dmrginp.outputlevel() > 0){
        pout << "Amplitude : " << sqrt(o)/(mps_nevpt::ZeroEnergy[baseState]- perturberEnergy) <<endl;
        pout << "Ener(only CAS part) : " << h/o<<endl;
        pout << "Energy : " << perturberEnergy<<endl;
        pout << "Correction Energy: "<< o/(mps_nevpt::ZeroEnergy[baseState]- perturberEnergy)<<endl; 
      }
    }
    else{
      if (dmrginp.outputlevel() > 0){
        pout << "Amplitude : " << 0.0 <<endl;
        pout << "Energy : " << 0.0<<endl;
      }
    }



  }
  pout << "Nevpt2 correction to the energy for state 0 in subspace Va is " << energy<<endl;;
  pout << "Nevpt2 Va subspace perturber Amplitude : " << overlap<<endl;;
  //pout << "Core Energy of nevpt2 " <<perturber::CoreEnergy[0]<<endl;
  std::string file = str(boost::format("%s%s%d") % dmrginp.load_prefix() % "/Va_" % baseState);
  std::fstream f(file,std::fstream::out);
  f << energy <<endl;
  f << overlap <<endl;
  f.close();
}

void SpinAdapted::mps_nevpt::type1::calcHamiltonianAndOverlap(const MPS& statea, double& h, double& o, perturber& pb) {

#ifndef SERIAL
  mpi::communicator world;
#endif


  SpinBlock system, siteblock;
  bool forward = true, restart=false, warmUp = false;
  int leftState=0, rightState=0, forward_starting_size=1, backward_starting_size=1, restartSize =0;
  InitBlocks::InitStartingBlock(system, forward, leftState, rightState, forward_starting_size, backward_starting_size, restartSize, restart, warmUp, 0,statea.getw().get_deltaQuantum(), statea.getw().get_deltaQuantum()); 

  if (dmrginp.outputlevel() > 0)
    pout << system<<endl;
  system.transform_operators(const_cast<std::vector<Matrix>&>(statea.getSiteTensors(0)), 
      		       const_cast<std::vector<Matrix>&>(statea.getSiteTensors(0)), false, false );

  int sys_add = true; bool direct = true; 

  for (int i=0; i<mps_nevpt::sweepIters-1; i++) {
    SpinBlock newSystem;
    system.addAdditionalCompOps();

    InitBlocks::InitNewSystemBlock(system, siteBlocks_noDES[i+1], newSystem, leftState, rightState, sys_add, direct, 0, DISTRIBUTED_STORAGE, false, true,NO_PARTICLE_SPIN_NUMBER_CONSTRAINT,statea.getw().get_deltaQuantum(),statea.getw().get_deltaQuantum());

    newSystem.transform_operators(const_cast<std::vector<Matrix>&>(statea.getSiteTensors(i+1)), 
      			    const_cast<std::vector<Matrix>&>(statea.getSiteTensors(i+1)), false );

    system = newSystem;
  }

  SpinBlock newSystem, big;
  system.addAdditionalCompOps();
  //system.printOperatorSummary();
  //To set implicit_transpose to true.
  //The last site spinblock should have implicit_transpose true.
  InitBlocks::InitNewSystemBlock(system, siteBlocks_noDES[mps_nevpt::sweepIters], newSystem,  leftState, rightState, sys_add, direct, 0, DISTRIBUTED_STORAGE, false, true,NO_PARTICLE_SPIN_NUMBER_CONSTRAINT,statea.getw().get_deltaQuantum(),statea.getw().get_deltaQuantum());
  
  newSystem.set_loopblock(false); system.set_loopblock(false);
  newSystem.addAdditionalCompOps();
  //siteBlocks_noDES[mps_nevpt::sweepIters+1].set_loopblock(false);
  InitBlocks::InitBigBlock(newSystem, siteBlocks_noDES[mps_nevpt::sweepIters+1], big,statea.getw().get_deltaQuantum(),statea.getw().get_deltaQuantum()); 
  

  //FIXME
  //Assume statea.getw() has only one deltaquantum.
  //Spin Embeding for zero order wavefunction is needed.

  Wavefunction temp = statea.getw();
  temp.Clear();
  big.multiplyH(const_cast<Wavefunction&>(statea.getw()), &temp, 1);

  if (mpigetrank() == 0)
    h = DotProduct(statea.getw(), temp);

  temp.Clear();
  big.multiplyOverlap(const_cast<Wavefunction&>(statea.getw()), &temp, 1);
  if (mpigetrank() == 0)
    o = DotProduct(statea.getw(), temp);
  if(dmrginp.spinAdapted())
  {

    double cg= 0.0;
    //TODO
    //Assume the zero order wavefunction has a spin zero. 
    //Spin Embeding must be used.
    SpinQuantum wQ= statea.getw().get_deltaQuantum(0);
    //cg*= dmrginp.get_ninej()(wQ.get_s().getirrep(), 1, dmrginp.effective_molecule_quantum().get_s().getirrep(), 
    //      					   pb.delta.get_s().getirrep(), pb.delta.get_s().getirrep(), 0,
    //                          dmrginp.effective_molecule_quantum().get_s().getirrep(), 0, dmrginp.effective_molecule_quantum().get_s().getirrep());
    //cg*= Symmetry::spatial_ninej(wQ.get_symm().getirrep(), -pb.delta.get_symm().getirrep(), dmrginp.effective_molecule_quantum().get_symm().getirrep(), 
    //      					   pb.delta.get_symm().getirrep(), -pb.delta.get_symm().getirrep(), 0,
    //                          dmrginp.effective_molecule_quantum().get_symm().getirrep(), 0, dmrginp.effective_molecule_quantum().get_symm().getirrep());
    cg += pow(clebsch(wQ.get_s().getirrep(),-1,1,1,0,0),2);
    cg += pow(clebsch(wQ.get_s().getirrep(),1,1,-1,0,0),2);
    //cout << "cg coefficient: " <<cg<<endl;
    h*= cg*cg;
    o*= cg*cg;
  }


#ifndef SERIAL
  mpi::broadcast(world, h, 0);
  mpi::broadcast(world, o, 0);
#endif

  return;
}

//
//void SpinAdapted::mps_nevpt::type1::Startup(const SweepParams &sweepParams, const bool &forward, perturber& pb, int baseState) {
//
//#ifndef SERIAL
//  mpi::communicator world;
//#endif
//  assert(forward);
//  SpinBlock system;
//  system.nonactive_orb() =pb.orb();
//  bool restart=false, warmUp = false;
//  int forward_starting_size=1, backward_starting_size=0, restartSize =0;
//  InitBlocks::InitStartingBlock(system, forward, pb.wavenumber(), baseState, forward_starting_size, backward_starting_size, restartSize, restart, warmUp, 0,pb.braquanta, pb.ketquanta); 
//
//  SpinBlock::store (forward, system.get_sites(), system, pb.wavenumber(), baseState); // if restart, just restoring an existing block --
//
//  for (int i=0; i<mps_nevpt::sweepIters; i++) {
//    SpinBlock newSystem;
//    SpinBlock dotSystem(i+1,i+1,pb.orb(),false);
//
//    system.addAdditionalCompOps();
//    //newSystem.default_op_components(true, system, dotSystem, true, true, false);
//    newSystem.perturb_op_components(false, system, dotSystem, pb);
//    newSystem.setstoragetype(DISTRIBUTED_STORAGE);
//    newSystem.BuildSumBlock(LessThanQ, system, dotSystem, pb.ketquanta, pb.ketquanta);
//    //newSystem.BuildSumBlock(LessThanQ, system, dotSystem, pb.braquanta, pb.ketquanta);
//    newSystem.printOperatorSummary();
//    //SpinBlock Environment, big;
//    //SpinBlock::restore (!forward, newSystem.get_complementary_sites() , Environment, baseState, baseState);
//    //TODO
//    //SpinBlock::restore (!forward, newSystem.get_complementary_sites() , Environment,sweepParams.current_root(),sweepParams.current_root());
//
//    //big.BuildSumBlock(PARTICLE_SPIN_NUMBER_CONSTRAINT, newSystem, Environment, pb.braquanta, pb.ketquanta);
//
//    //StateInfo envStateInfo;
//    StateInfo ketStateInfo;
//    StateInfo braStateInfo;
//    StateInfo halfbraStateInfo;// It has the same left and right StateInfo as braStateInfo. However, its total quanta is pb.ketquanta.
//    // It is used to project solution into to braStateInfo.
//
//    //TensorProduct (newSystem.get_braStateInfo(), *(ketStateInfo.rightStateInfo), pb.braquanta[0], EqualQ, braStateInfo);
//    //TODO
//    //TensorProduct do not support const StateInfo&
//
//    //StateInfo::restore(forward, environmentsites, envStateInfo, baseState);
//
//    //DiagonalMatrix e;
//    //if(i == 0)
//    //  GuessWave::guess_wavefunctions(solution, e, big, TRANSPOSE, true, true, 0.0, baseState); 
//    //else
//    //  GuessWave::guess_wavefunctions(solution, e, big, TRANSFORM, true, true, 0.0, baseState); 
//
//
//    //SpinAdapted::operatorfunctions::Product(&newSystem, ccd, solution[0], &ketStateInfo, stateb.getw(), temp, SpinQuantum(0, SpinSpace(0), IrrepSpace(0)), true, 1.0);
//
//    
//
//    boost::shared_ptr<SparseMatrix> O;
//    if (pb.type() == TwoPerturbType::Va)
//      O = newSystem.get_op_array(CDD_SUM).get_local_element(0)[0]->getworkingrepresentation(&newSystem);
//    if (pb.type() == TwoPerturbType::Vi)
//      O = newSystem.get_op_array(CCD_SUM).get_local_element(0)[0]->getworkingrepresentation(&newSystem);
//    boost::shared_ptr<SparseMatrix> overlap = newSystem.get_op_array(OVERLAP).get_local_element(0)[0]->getworkingrepresentation(&newSystem);
//    std::vector<Matrix> brarotateMatrix, ketrotateMatrix;
//    LoadRotationMatrix (newSystem.get_sites(), ketrotateMatrix, baseState);
//
//    newSystem.transform_operators(ketrotateMatrix,ketrotateMatrix);
//    SpinBlock::store (forward, newSystem.get_sites(), newSystem, pb.wavenumber(), baseState); // if restart, just restoring an existing block --
//    system=newSystem;
//  }
//  //TODO
//  //It seems that there is no need to do Last Step of Sweep.
//}
//

void SpinAdapted::mps_nevpt::type1::Startup(const SweepParams &sweepParams, const bool &forward, perturber& pb, int baseState) {

#ifndef SERIAL
  mpi::communicator world;
#endif
  assert(forward);
  SpinBlock system;
  system.nonactive_orb() =pb.orb();
  bool restart=false, warmUp = false;
  int forward_starting_size=1, backward_starting_size=0, restartSize =0;
  InitBlocks::InitStartingBlock(system, forward, pb.wavenumber(), baseState, forward_starting_size, backward_starting_size, restartSize, restart, warmUp, 0,pb.braquanta, pb.ketquanta); 

  SpinBlock::store (forward, system.get_sites(), system, pb.wavenumber(), baseState); // if restart, just restoring an existing block --

  for (int i=0; i<mps_nevpt::sweepIters; i++) {
    SpinBlock newSystem;
    SpinBlock dotSystem(i+1,i+1,pb.orb(),false);

    system.addAdditionalCompOps();
    //newSystem.default_op_components(true, system, dotSystem, true, true, false);
    newSystem.perturb_op_components(false, system, dotSystem, pb);
    newSystem.setstoragetype(DISTRIBUTED_STORAGE);
    newSystem.BuildSumBlock(LessThanQ, system, dotSystem, pb.braquanta, pb.ketquanta);
    newSystem.printOperatorSummary();
    //SpinBlock Environment, big;
    //SpinBlock::restore (!forward, newSystem.get_complementary_sites() , Environment, baseState, baseState);
    //TODO
    //SpinBlock::restore (!forward, newSystem.get_complementary_sites() , Environment,sweepParams.current_root(),sweepParams.current_root());

    //big.BuildSumBlock(PARTICLE_SPIN_NUMBER_CONSTRAINT, newSystem, Environment, pb.braquanta, pb.ketquanta);

    //StateInfo envStateInfo;
    StateInfo ketStateInfo;
    StateInfo braStateInfo;
    StateInfo halfbraStateInfo;// It has the same left and right StateInfo as braStateInfo. However, its total quanta is pb.ketquanta.
    // It is used to project solution into to braStateInfo.

    std::vector<Wavefunction> solution; solution.resize(1);
    std::vector<Wavefunction> outputState; outputState.resize(1);
    std::vector<Wavefunction> solutionprojector; solutionprojector.resize(1);
    solution[0].LoadWavefunctionInfo(ketStateInfo, newSystem.get_sites(), baseState);
    #ifndef SERIAL
      broadcast(world, ketStateInfo, 0);
      broadcast(world, solution, 0);
    #endif
    outputState[0].AllowQuantaFor(newSystem.get_braStateInfo(), *(ketStateInfo.rightStateInfo), pb.braquanta);
    outputState[0].set_onedot(solution[0].get_onedot());
    outputState[0].Clear();
    solutionprojector[0].AllowQuantaFor(newSystem.get_braStateInfo(), *(ketStateInfo.rightStateInfo), pb.ketquanta);
    solutionprojector[0].set_onedot(solution[0].get_onedot());
    solutionprojector[0].Clear();
    //TensorProduct (newSystem.get_braStateInfo(), *(ketStateInfo.rightStateInfo), pb.braquanta[0], EqualQ, braStateInfo);
    //TODO
    //TensorProduct do not support const StateInfo&
    TensorProduct (newSystem.set_braStateInfo(), *(ketStateInfo.rightStateInfo), pb.braquanta[0], EqualQ, braStateInfo);
    TensorProduct (newSystem.set_braStateInfo(), *(ketStateInfo.rightStateInfo), pb.ketquanta[0], EqualQ, halfbraStateInfo);

    //StateInfo::restore(forward, environmentsites, envStateInfo, baseState);

    //DiagonalMatrix e;
    //if(i == 0)
    //  GuessWave::guess_wavefunctions(solution, e, big, TRANSPOSE, true, true, 0.0, baseState); 
    //else
    //  GuessWave::guess_wavefunctions(solution, e, big, TRANSFORM, true, true, 0.0, baseState); 


    //SpinAdapted::operatorfunctions::Product(&newSystem, ccd, solution[0], &ketStateInfo, stateb.getw(), temp, SpinQuantum(0, SpinSpace(0), IrrepSpace(0)), true, 1.0);

    

    boost::shared_ptr<SparseMatrix> O;
    if (pb.type() == TwoPerturbType::Va)
      O = newSystem.get_op_array(CDD_SUM).get_local_element(0)[0]->getworkingrepresentation(&newSystem);
    if (pb.type() == TwoPerturbType::Vi)
      O = newSystem.get_op_array(CCD_SUM).get_local_element(0)[0]->getworkingrepresentation(&newSystem);
    boost::shared_ptr<SparseMatrix> overlap = newSystem.get_op_array(OVERLAP).get_local_element(0)[0]->getworkingrepresentation(&newSystem);
    SpinAdapted::operatorfunctions::TensorMultiply(*O, &braStateInfo, &ketStateInfo , solution[0], outputState[0], pb.delta, true, 1.0);
    SpinAdapted::operatorfunctions::TensorMultiply(*overlap, &halfbraStateInfo, &ketStateInfo , solution[0], solutionprojector[0], overlap->get_deltaQuantum(0), true, 1.0);
    DensityMatrix bratracedMatrix(newSystem.get_braStateInfo());
    bratracedMatrix.allocate(newSystem.get_braStateInfo());
    double norm = DotProduct(outputState[0], outputState[0]);
    if(norm > NUMERICAL_ZERO)
      SpinAdapted::operatorfunctions::MultiplyProduct(outputState[0], Transpose(const_cast<Wavefunction&> (outputState[0])), bratracedMatrix, 0.5/norm);
    SpinAdapted::operatorfunctions::MultiplyProduct(solutionprojector[0], Transpose(const_cast<Wavefunction&> (solutionprojector[0])), bratracedMatrix, 0.5);
    std::vector<Matrix> brarotateMatrix, ketrotateMatrix;
    LoadRotationMatrix (newSystem.get_sites(), ketrotateMatrix, baseState);
    double error;
    if (!mpigetrank())
      error = makeRotateMatrix(bratracedMatrix, brarotateMatrix, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());
    #ifndef SERIAL
      broadcast(world, ketrotateMatrix, 0);
      broadcast(world, brarotateMatrix, 0);
    #endif

    SaveRotationMatrix (newSystem.get_sites(), brarotateMatrix, pb.wavenumber());
    newSystem.transform_operators(brarotateMatrix,ketrotateMatrix);
    SpinBlock::store (forward, newSystem.get_sites(), newSystem, pb.wavenumber(), baseState); // if restart, just restoring an existing block --
    system=newSystem;
  }
  //TODO
  //It seems that there is no need to do Last Step of Sweep.
}


