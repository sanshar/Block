#include "fciqmchelper.h"
#include "rotationmat.h"
#include "operatorfunctions.h"
#include "spinblock.h"
#include "initblocks.h"

namespace SpinAdapted{

  //initializing the static variables
  int MPS::sweepIters ;
  bool MPS::spinAdapted ;
  std::vector<SpinBlock> MPS::siteBlocks;

  
  //assumes that the state has been canonicalized in the left canonical form already and stored on disk
  MPS::MPS(int stateindex) {

    std::vector<int> rotSites(2,0);
    if (!dmrginp.spinAdapted()) {
      rotSites[0] = 0; rotSites[1] = 3; //this should be 1 for spinadapted
    }
    else {
      rotSites[0] = 0; rotSites[1] = 1; //this should be 1 for spinadapted
    }

    Matrix m(1,1); m=0;
    std::vector<Matrix> rotMat;
    if (!dmrginp.spinAdapted()) {
      rotMat = std::vector<Matrix>(4,m);
      rotMat[0](1,1) = 1.0;rotMat[1](1,1) = 1.0;rotMat[2](1,1) = 1.0;rotMat[3](1,1) = 1.0;
    }
    else {
      rotMat = std::vector<Matrix>(3,m);
      rotMat[0](1,1) = 1.0;rotMat[1](1,1) = 1.0;rotMat[2](1,1) = 1.0;
    }

    SiteTensors.push_back(rotMat);

    for (int i=0; i<MPS::sweepIters; i++) {
      LoadRotationMatrix(rotSites, rotMat, stateindex);
      SiteTensors.push_back(rotMat);
      if (!dmrginp.spinAdapted())
	rotSites[1] += 2;
      else
	rotSites[1] += 1;
    }

    //loading the final wavefunction
    if (!dmrginp.spinAdapted())
      rotSites[1] -=2;
    else
      rotSites[1] -=1;

    StateInfo s;
    w.LoadWavefunctionInfo(s, rotSites, stateindex);
  }


  //this is a helper function to make a MPS from a occupation number representation of determinant
  void MPS::Init(std::vector<bool>& occnum)
  {
    assert(occnum.size() == dmrginp.last_site());
    Matrix m(1,1); m=1.0;
    Matrix dummy;

    //first rotation matrix 
    std::vector<Matrix> rotMat; rotMat.resize(4, dummy);
    int index = occnum[0]*2+occnum[1];
    rotMat[index]=m;
    SiteTensors.push_back(rotMat);
    SpinQuantum sTotal = MPS::siteBlocks[0].get_stateInfo().quanta[index];

    for (int i=0; i<MPS::sweepIters-1; i++) {
      //stateinfo of in incoming bond of dimension 1
      SpinQuantum sq[] = {sTotal}; int qs[] = {1}; int n = 1;
      StateInfo stateTotal(n, sq, qs), currentState;

      TensorProduct(stateTotal, const_cast<StateInfo&>(MPS::siteBlocks[i+1].get_stateInfo()), currentState, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);

      std::vector<Matrix> rotMat; rotMat.resize(currentState.quanta.size(), dummy);
      int index = occnum[2*i+2]*2+occnum[2*i+3];
      index = currentState.quantaMap(0, index)[0];
      rotMat[index]=m;

      sTotal = currentState.quanta[index];
      SiteTensors.push_back(rotMat);
    }

    //stateinfo of in incoming bond of dimension 1
    SpinQuantum sq[] = {sTotal}; int qs[] = {1}; int n = 1;
    StateInfo stateTotal(n, sq, qs), secondLastState;

    //the incoming bond x k-1 site stateinfo
    TensorProduct(stateTotal, const_cast<StateInfo&>(MPS::siteBlocks[MPS::sweepIters].get_stateInfo()), 
		  secondLastState, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
    int index1 =  occnum[2*MPS::sweepIters]*2+occnum[2*MPS::sweepIters+1];
    index1 = secondLastState.quantaMap(0, index1)[0];

    //now make wavefunction with the big state A
    w.AllowQuantaFor(secondLastState, MPS::siteBlocks[MPS::sweepIters+1].get_stateInfo(), 
		     std::vector<SpinQuantum>(1,dmrginp.effective_molecule_quantum()));

    int index2 = occnum[2*MPS::sweepIters+2]*2+occnum[2*MPS::sweepIters+3];

    w(index1, index2)(1,1) = 1.0;
  }


  //this is a helper function to make a MPS from a occupation number representation of determinant
  //this representation is slightly different than the usual occupation, here each integer
  //element is a spatial orbital which can have a value 0, -1, 1, or 2.
  MPS::MPS(long *occnum, int length)
  {
    assert(length*64 >= dmrginp.last_site());

    //convert the int array into a vector<bool>
    std::vector<bool> occ(dmrginp.last_site(), 0);

    ulong temp = 1;
    int index = 0;
    for (int i=0; i <length ; i++) 
      for (int j=63; j>=0; j--) {
	if (index >=dmrginp.last_site()) break;

	occ[index] = occnum[i] & ( temp << j ) ;
	index++;
      }

    Init(occ);
  }


  MPS::MPS(std::vector<bool>& occ) {
    Init(occ);
  }


  double calculateOverlap(const MPS& statea, const MPS& stateb) {

    Overlap siteOverlap;
    siteOverlap.makeIdentity(MPS::siteBlocks[0].get_stateInfo());

    StateInfo stateA=MPS::siteBlocks[0].get_stateInfo(), stateB=MPS::siteBlocks[0].get_stateInfo();
    Overlap o = siteOverlap;

    for (int i=0; i<MPS::sweepIters; i++) {

      StateInfo renormA, renormB;
      SpinAdapted::StateInfo::transform_state(statea.getSiteTensors(i), stateA, renormA);
      SpinAdapted::StateInfo::transform_state(stateb.getSiteTensors(i), stateB, renormB);
      o.renormalise_transform(statea.getSiteTensors(i), &renormA, stateb.getSiteTensors(i), &renormB);

      //make overlap and state info for the current site
      siteOverlap.makeIdentity(MPS::siteBlocks[i+1].get_stateInfo());
      
      //take tensor product of stateinfo
      StateInfo A, B;
      TensorProduct(renormA, const_cast<StateInfo&>(MPS::siteBlocks[i+1].get_stateInfo()), 
		    A, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
      TensorProduct(renormB, const_cast<StateInfo&>(MPS::siteBlocks[i+1].get_stateInfo()), 
		    B, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
      A.CollectQuanta(); B.CollectQuanta();
      stateA = A; stateB = B;

      //build the new Overlap matrix
      Overlap onew;
      onew.set_deltaQuantum(1, SpinQuantum(0, SpinSpace(0), IrrepSpace(0)));
      onew.set_built()=true; onew.allocate(stateA, stateB); 
      onew.set_initialised() = true;
      SpinAdapted::operatorfunctions::TensorProduct(o, siteOverlap, &stateA, &stateB, onew, 1.0, true);
      o = onew;
    }
    Wavefunction temp = statea.getw();
    temp.Clear();

    siteOverlap.makeIdentity(MPS::siteBlocks[MPS::sweepIters+1].get_stateInfo());

    StateInfo braStateInfo, ketStateInfo;
    TensorProduct(stateA, const_cast<StateInfo&>(MPS::siteBlocks[MPS::sweepIters+1].get_stateInfo()), 
		  braStateInfo, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
    TensorProduct(stateB, const_cast<StateInfo&>(MPS::siteBlocks[MPS::sweepIters+1].get_stateInfo())
		  , ketStateInfo, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);

    SpinAdapted::operatorfunctions::TensorMultiply(o, siteOverlap, &braStateInfo, &ketStateInfo, stateb.getw(), temp, SpinQuantum(0, SpinSpace(0), IrrepSpace(0)), true, 1.0);
    double overlap = DotProduct(statea.getw(), temp);
    return overlap;
  }


  void calcHamiltonianAndOverlap(const MPS& statea, const MPS& stateb, double& h, double& o) {

    SpinBlock system = SpinBlock(0,0,false), siteblock;
    SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));

    system.transform_operators(const_cast<std::vector<Matrix>&>(statea.getSiteTensors(0)), 
			       const_cast<std::vector<Matrix>&>(stateb.getSiteTensors(0)) );

    int sys_add = true; bool direct = true; 

    for (int i=0; i<MPS::sweepIters-1; i++) {
      SpinBlock newSystem;
      InitBlocks::InitNewSystemBlock(system, MPS::siteBlocks[i+1], newSystem, 0, 1, sys_add, direct, DISTRIBUTED_STORAGE, false, true);
      newSystem.transform_operators(const_cast<std::vector<Matrix>&>(statea.getSiteTensors(i+1)), 
				    const_cast<std::vector<Matrix>&>(stateb.getSiteTensors(i+1)), false );

      system = newSystem;
    }

    SpinBlock newSystem, big;
    InitBlocks::InitNewSystemBlock(system, MPS::siteBlocks[MPS::sweepIters], newSystem, 0, 1, sys_add, direct, DISTRIBUTED_STORAGE, false, true);
    
    InitBlocks::InitBigBlock(newSystem, MPS::siteBlocks[MPS::sweepIters+1], big); 
    
    Wavefunction temp = statea.getw();
    temp.Clear();

    big.multiplyH(const_cast<Wavefunction&>(stateb.getw()), &temp, 1);
    h = DotProduct(statea.getw(), temp);

    temp.Clear();
    big.multiplyOverlap(const_cast<Wavefunction&>(stateb.getw()), &temp, 1);
    o = DotProduct(statea.getw(), temp);
    return;
  }


}

      
      


