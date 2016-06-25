#include "type1.h"
#include "mps_nevpt.h"
#include "sweep.h"
#include "npdm.h"

void dmrg(double sweep_tol);

namespace SpinAdapted{
namespace mps_nevpt{
  std::vector<SpinBlock> siteBlocks_noDES;
  vector<double> ZeroEnergy;
  int sweepIters ;
  void readZeroEnergy(){
    ZeroEnergy.resize(dmrginp.nroots());
  //  CoreEnergy.resize(dmrginp.nroots(),0.0);
  
    std::string efile;
    efile = str(boost::format("%s%s") % dmrginp.load_prefix() % "/dmrg.e" );
    FILE* f = fopen(efile.c_str(), "rb");      
    for(int j=0;j<dmrginp.nroots();++j) {
      fread( &(ZeroEnergy[j]), 1, sizeof(double), f);
    }
    fclose(f);
  
  //  for(int i=0; i< dmrginp.core_size(); i++)
  //    CoreEnergy[0] += 2*v_1[0](2*(i+dmrginp.act_size()),2*(i+dmrginp.act_size()));
  //  //ZeroEnergy[0] +=CoreEnergy[0];
  //  pout << "Zero order energy for state 0 is " << ZeroEnergy[0]<<endl;;
  }
}
}

using namespace SpinAdapted;


void SpinAdapted::mps_nevpt::mps_nevpt(double sweep_tol)
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
//    Npdm::npdm(NPDM_ONEPDM);
//    dmrginp.do_pdm() = false;
    dmrg(sweep_tol);
    dmrginp.calc_type() = MPS_NEVPT;
  }
  dmrginp.calc_type() = DMRG;
  if(!dmrginp.spinAdapted())
  {
    mps_nevpt::sweepIters = dmrginp.last_site()/2-2;
    MPS::sweepIters = dmrginp.last_site()/2-2;
    for (int j=0; j<dmrginp.last_site(); j+=2){
      siteBlocks_noDES.push_back(SpinBlock(j/2, j/2, 0, true)); 
//    MPS::siteBlocks.push_back(SpinBlock(j, j, 0, false)); 
    }
  }
  else
  {
    mps_nevpt::sweepIters = dmrginp.last_site()-2;
    MPS::sweepIters = dmrginp.last_site()-2;
    for (int j=0; j<dmrginp.last_site(); j++){
      siteBlocks_noDES.push_back(SpinBlock(j, j, 0, true)); 
//      MPS::siteBlocks.push_back(SpinBlock(j, j, 0, false)); 
    }
  }
  dmrginp.calc_type() = MPS_NEVPT;
  SweepParams sweepParams;
  bool direction = true;
  algorithmTypes atype = dmrginp.algorithm_method();
  dmrginp.set_algorithm_method() = ONEDOT;
  //initialize state info and canonicalize wavefunction is always done using onedot algorithm
  int baseState = dmrginp.nevpt_state_num();
  if (mpigetrank()==0) {
    Sweep::InitializeStateInfo(sweepParams, direction, baseState);
    Sweep::InitializeStateInfo(sweepParams, !direction, baseState);
    Sweep::CanonicalizeWavefunction(sweepParams, direction, baseState);
    Sweep::CanonicalizeWavefunction(sweepParams, !direction, baseState);
    Sweep::CanonicalizeWavefunction(sweepParams, direction, baseState);
  }
  dmrginp.set_algorithm_method() = atype;
  SpinAdapted::mps_nevpt::readZeroEnergy();
  SpinAdapted::mps_nevpt::type1::subspace_Vi(baseState);
  SpinAdapted::mps_nevpt::type1::subspace_Va(baseState);
}
