#include "type1.h"
#include "mps_nevpt.h"
#include "sweep.h"
#include "npdm.h"

void dmrg(double sweep_tol);
vector<double> perturber::ZeroEnergy;
vector<double> perturber::CoreEnergy;

double readZeroEnergy(){
  perturber::ZeroEnergy.resize(dmrginp.nroots());
  perturber::CoreEnergy.resize(dmrginp.nroots(),0.0);

  std::string efile;
  efile = str(boost::format("%s%s") % dmrginp.load_prefix() % "/dmrg.e" );
  FILE* f = fopen(efile.c_str(), "rb");      
  for(int j=0;j<dmrginp.nroots();++j) {
    fread( &(perturber::ZeroEnergy[j]), 1, sizeof(double), f);
  }
  fclose(f);

  for(int i=0; i< dmrginp.core_size(); i++)
    perturber::CoreEnergy[0] += 2*v_1[0](2*(i+dmrginp.act_size()),2*(i+dmrginp.act_size()));
  //perturber::ZeroEnergy[0] +=perturber::CoreEnergy[0];
  pout << "Zero order energy for state 0 is " << perturber::ZeroEnergy[0]<<endl;;
}

void SpinAdapted::mps_nevpt::mps_nevpt(int baseState)
{
  //if(!restartpdm){
  //  if (RESTART && !FULLRESTART)
  //    restart(sweep_tol, reset_iter);
  //  else if (FULLRESTART) {
  //    fullrestartGenblock();
  //    reset_iter = true;
  //    sweepParams.restorestate(direction, restartsize);
  //    sweepParams.calc_niter();
  //    sweepParams.savestate(direction, restartsize);
  //    restart(sweep_tol, reset_iter);
  //  }
  //  else {
  //    dmrg(sweep_tol);
  //  }
  //}

  if(dmrginp.calc_type() == MPS_NEVPT){
  //  double sweep_tol = 1e-7;
    dmrginp.calc_type() = DMRG;
    Npdm::npdm(NPDM_ONEPDM);
    dmrginp.do_pdm() = false;
  //  dmrg(sweep_tol);
    dmrginp.calc_type() = MPS_NEVPT;
  }
  dmrginp.calc_type() = DMRG;
  if(!dmrginp.spinAdapted())
  {
    MPS::sweepIters = dmrginp.last_site()/2-2;
    for (int j=0; j<dmrginp.last_site(); j+=2){
      MPS::siteBlocks_noDES.push_back(SpinBlock(j/2, j/2, 0, true)); 
//    MPS::siteBlocks.push_back(SpinBlock(j, j, 0, false)); 
    }
  }
  else
  {
    MPS::sweepIters = dmrginp.last_site()-2;
    for (int j=0; j<dmrginp.last_site(); j++){
      MPS::siteBlocks_noDES.push_back(SpinBlock(j, j, 0, true)); 
//      MPS::siteBlocks.push_back(SpinBlock(j, j, 0, false)); 
    }
  }
  dmrginp.calc_type() = MPS_NEVPT;
  SweepParams sweepParams;
  bool direction = true;
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
  readZeroEnergy();
  SpinAdapted::mps_nevpt::type1::subspace_Vi(baseState);
  SpinAdapted::mps_nevpt::type1::subspace_Va(baseState);
}
