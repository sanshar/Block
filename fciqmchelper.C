#include "fciqmchelper.h"
#include "rotationmat.h"
#include "operatorfunctions.h"

namespace SpinAdapted{
  
  //assumes that the state has been canonicalized in the left canonical form already and stored on disk
  MPS::MPS(int stateindex) {

    std::vector<int> rotSites(2,0);
    rotSites[0] = 0; rotSites[1] = 3; //this should be 1 for spinadapted


    Matrix m(1,1); m=0;
    std::vector<Matrix> rotMat(4,m);
    rotMat[0](1,1) = 1.0;rotMat[1](1,1) = 1.0;rotMat[2](1,1) = 1.0;rotMat[3](1,1) = 1.0;

    SiteTensors.push_back(rotMat);

    for (int i=0; i<MPS::sweepIters; i++) {
      LoadRotationMatrix(rotSites, rotMat, stateindex);
      SiteTensors.push_back(rotMat);

      Wavefunction w1;StateInfo s;
      w1.LoadWavefunctionInfo(s, rotSites, stateindex);
      cout << w1<<endl;

      rotSites[1] += 2;
    }
    //loading the final wavefunction
    rotSites[1] -=2;

    StateInfo s;
    w.LoadWavefunctionInfo(s, rotSites, stateindex);
  }



  //this is a helper function to make a MPS from a occupation number representation of determinant
  MPS::MPS(std::vector<bool>& occnum)
  {
    assert(occnum.size() == dmrginp.last_site());
    Matrix m(1,1); m=1.0;
    Matrix dummy;

    StateInfo siteInfo;
    makeStateInfo(siteInfo, 0);

    //first rotation matrix which does not get added to the rotSites
    std::vector<Matrix> rotMat; rotMat.resize(4, dummy);
    int index = occnum[0]*2+occnum[1];
    rotMat[index]=m;
    SiteTensors.push_back(rotMat);
    SpinQuantum sTotal = siteInfo.quanta[index];

    for (int i=0; i<MPS::sweepIters-1; i++) {
      //stateinfo of in incoming bond of dimension 1
      SpinQuantum sq[] = {sTotal}; int qs[] = {1}; int n = 1;
      StateInfo stateTotal(n, sq, qs), currentState;

      makeStateInfo(siteInfo, i+1);
      TensorProduct(stateTotal, siteInfo, currentState, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);

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

    //second to last site stateinfo
    makeStateInfo(siteInfo, MPS::sweepIters);

    //the incoming bond x k-1 site stateinfo
    TensorProduct(stateTotal, siteInfo, secondLastState, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
    int index1 =  occnum[2*MPS::sweepIters]*2+occnum[2*MPS::sweepIters+1];
    index1 = secondLastState.quantaMap(0, index1)[0];

    //last site state info
    StateInfo finalSite;
    makeStateInfo(finalSite, MPS::sweepIters+1);
    StateInfo A;

    //A will be the stateinfo the big block, with secondLastState as left and final state as right
    TensorProduct(secondLastState, finalSite, A, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);

    //now make wavefunction with the big state A
    w.AllowQuantaFor(secondLastState, finalSite, std::vector<SpinQuantum>(1,dmrginp.effective_molecule_quantum()));

    int index2 = occnum[2*MPS::sweepIters+2]*2+occnum[2*MPS::sweepIters+3];

    w(index1, index2)(1,1) = 1.0;
  }


  double calculateOverlap(const MPS& statea, const MPS& stateb) {

    Overlap siteOverlap;
    StateInfo siteInfo;
    makeStateInfo(siteInfo, 0);
    siteOverlap.makeIdentity(siteInfo);

    StateInfo stateA=siteInfo, stateB=siteInfo;
    Overlap o = siteOverlap;

    for (int i=0; i<MPS::sweepIters; i++) {

      StateInfo renormA, renormB;
      cout << stateA<<endl;
      cout << o <<endl;
      SpinAdapted::StateInfo::transform_state(statea.getSiteTensors(i), stateA, renormA);
      SpinAdapted::StateInfo::transform_state(stateb.getSiteTensors(i), stateB, renormB);
      o.renormalise_transform(statea.getSiteTensors(i), &renormA, stateb.getSiteTensors(i), &renormB);

      //make overlap and state info for the current site
      makeStateInfo(siteInfo, i+1);
      siteOverlap.makeIdentity(siteInfo);
      
      //take tensor product of stateinfo
      StateInfo A, B;
      TensorProduct(renormA, siteInfo, A, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
      TensorProduct(renormB, siteInfo, B, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
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

    makeStateInfo(siteInfo, MPS::sweepIters+1);
    siteOverlap.makeIdentity(siteInfo);

    StateInfo braStateInfo, ketStateInfo;
    TensorProduct(stateA, siteInfo, braStateInfo, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
    TensorProduct(stateB, siteInfo, ketStateInfo, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);

    SpinAdapted::operatorfunctions::TensorMultiply(o, siteOverlap, &braStateInfo, &ketStateInfo, stateb.getw(), temp, SpinQuantum(0, SpinSpace(0), IrrepSpace(0)), true, 1.0);
    double overlap = DotProduct(statea.getw(), temp);
    return overlap;
  }
}

      
      


