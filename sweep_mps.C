/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "guess_wavefunction.h"
#include "sweepgenblock.h"
#include "sweep.h"
#ifdef USE_BTAS
#include "btas/SPARSE/STArray.h"
#include "btas/SPARSE/SDcontract.h"
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
#include <btas/TVector.h>
#include "overlaptensor.h"

using namespace boost;
using namespace std;

//This take the sites and just makes a stateinfo containing all the states in it
void makeStateInfo(StateInfo& s, int site)
{
  std::vector< std::vector<Csf> > ladders;
  std::vector< Csf > dets; 
  std::vector<int> new_sites(1, site), sites;
  if (dmrginp.spinAdapted()) {
    sites = new_sites;
    dets = CSFUTIL::spinfockstrings(new_sites, ladders);
  }
  else {
    for (int i=0; i<new_sites.size(); i++) {
      sites.push_back( dmrginp.spatial_to_spin()[new_sites[i]]   );
      sites.push_back( dmrginp.spatial_to_spin()[new_sites[i]]+1 );
    }
    dets = CSFUTIL::spinfockstrings(new_sites);
    for (int j=0; j<dets.size(); j++)
      ladders.push_back(std::vector<Csf>(1,dets[j]));
  }
  s = StateInfo(dets);

  if (dmrginp.add_noninteracting_orbs() && dmrginp.molecule_quantum().get_s().getirrep() != 0 && dmrginp.spinAdapted() && site == 0)
  {
    SpinQuantum sq = dmrginp.molecule_quantum();
    sq = SpinQuantum(sq.get_s().getirrep(), sq.get_s(), IrrepSpace(0));
    int qs = 1, ns = 1;
    StateInfo addstate(ns, &sq, &qs), newstate; 

    TensorProduct(s, addstate, newstate, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
    s = newstate;
  }

}



void getComplementarySites(std::vector<int> &sites, std::vector<int> &complementarySites) 
{
  for (int i=0; i<dmrginp.last_site(); ++i)
    if (find(sites.begin(), sites.end(), i) == sites.end())
      complementarySites.push_back(i);
  return;
}

//reads the overlap matrix 
void readOverlap(const std::vector<int>& sites, btas::STArray<double, 2>& Overlap, int iState, int currentState)
{
  LoadOverlapTensor(sites, Overlap, iState, currentState);
}

void SpinAdapted::Sweep::makeDMRGOverlapFromBTASOverlap(SparseMatrix &Overlap, const std::vector<int>& sites, int left, int right)
{
  btas::STArray<double, 2> OverlapBtas;
  LoadOverlapTensor(sites, OverlapBtas, left, right);

  for (int i=0; i<Overlap.nrows(); i++)
    for (int j=0; j<Overlap.ncols(); j++)
      if (Overlap.allowed(i,j) && OverlapBtas.find(btas::make_array(i,j)) != OverlapBtas.end()) 
	DCOPY(Overlap(i,j).Storage(), &OverlapBtas.find(btas::make_array(i,j))->second->operator()(0,0), 1, Overlap(i,j).Store(), 1);
}


//makes the Btas overlap from dmrg overlap. Requires the btas overlap to have the block sparse structure initialized
void SpinAdapted::Sweep::makeBTASOverlapFromDMRGOverlap(const SparseMatrix &Overlap, btas::STArray<double, 2>& OverlapBtas)
{
  for (int i=0; i<Overlap.nrows(); i++)
    for (int j=0; j<Overlap.ncols(); j++)
      if (Overlap.allowed(i,j)) {
	OverlapBtas.reserve(btas::make_array(i,j));
	DCOPY(Overlap(i,j).Storage(), Overlap(i,j).Store(), 1, &OverlapBtas.find(btas::make_array(i,j))->second->operator()(0,0), 1);
      }
}


//takes an lower energy wavefunction and then contracts the right and left overlap matrices to give back a (tensor/wavefunction which can be contracted with the current wavefunction to get a overlap
void SpinAdapted::Sweep::getLowerStatesBlockRow(int currentState, const std::vector<int>& sites, const std::vector<int>& complementSites, std::vector<Wavefunction>& lowerStates, const StateInfo& leftState, const StateInfo& rightState, const vector<StateInfo>& stateInfoi)
{
  lowerStates.resize(currentState);
  for (int istate = 0; istate<currentState; istate++) {

    btas::STArray<double,2> fOverlap, rOverlap; 		     
    std::vector<int> prevSites;
    if (sites[0] == 0) {
      prevSites = sites; prevSites.pop_back();
    }
    else 
      prevSites = std::vector<int>((++sites.begin()), sites.end());

    readOverlap(prevSites, fOverlap, istate, currentState);
    readOverlap(complementSites, rOverlap, istate, currentState);

    Wavefunction w1 = lowerStates[istate];

    btas::TVector<btas::Dshapes,3> w1shape, w2shape, inter1;
    w1shape[0] = stateInfoi[istate].leftStateInfo->leftStateInfo->quantaStates;w1shape[1] = stateInfoi[istate].leftStateInfo->rightStateInfo->quantaStates;w1shape[2] = stateInfoi[istate].rightStateInfo->quantaStates;
    w2shape[0] = fOverlap.dshape()[1];w2shape[1] = stateInfoi[istate].leftStateInfo->rightStateInfo->quantaStates;w2shape[2] = rOverlap.dshape()[1];
    inter1[0] = fOverlap.dshape()[1];inter1[1] = stateInfoi[istate].leftStateInfo->rightStateInfo->quantaStates;inter1[2] = stateInfoi[istate].rightStateInfo->quantaStates;
    btas::STArray<double, 3> w1Tensor(w1shape, false), intermediate(inter1, false), w2Tensor(w2shape, false);

    w1.UnCollectQuantaAlongRows(*stateInfoi[istate].leftStateInfo, *stateInfoi[istate].rightStateInfo, w1Tensor);
    btas::SDcontract(1.0, fOverlap, btas::shape(0), w1Tensor, btas::shape(0), 0.0, intermediate);
    btas::SDcontract(1.0, intermediate, btas::shape(2), rOverlap, btas::shape(0), 0.0, w2Tensor);

    lowerStates[istate].Clear();
    lowerStates[istate].CollectQuantaAlongRows(*leftState.unCollectedStateInfo, rightState, w2Tensor, w1.get_deltaQuantum());
    //Normalise(lowerStates[istate]);
    
  }

}

//takes an lower energy wavefunction and then contracts the right and left overlap matrices to give back a (tensor/wavefunction which can be contracted with the current wavefunction to get a overlap
void SpinAdapted::Sweep::getLowerStatesBlockCol(int currentState, const std::vector<int>& sites, const std::vector<int>& complementSites, std::vector<Wavefunction>& lowerStates, const StateInfo& leftState, const StateInfo& rightState, const vector<StateInfo>& stateInfoi)
{
  lowerStates.resize(currentState);
  for (int istate = 0; istate<currentState; istate++) {

    btas::STArray<double,2> fOverlap, rOverlap; 		     
    std::vector<int> prevSites;
    if (complementSites[0] == 0) {
      prevSites = complementSites; prevSites.pop_back();
    }
    else 
      prevSites = std::vector<int>((++complementSites.begin()), complementSites.end());

    readOverlap(sites, fOverlap, istate, currentState);
    readOverlap(prevSites, rOverlap, istate, currentState);

    Wavefunction w1 = lowerStates[istate];

    btas::TVector<btas::Dshapes,3> w1shape, w2shape, inter1;
    w1shape[0] = stateInfoi[istate].leftStateInfo->leftStateInfo->quantaStates;w1shape[1] = stateInfoi[istate].leftStateInfo->rightStateInfo->quantaStates;w1shape[2] = stateInfoi[istate].rightStateInfo->quantaStates;
    w2shape[0] = fOverlap.dshape()[1];w2shape[1] = stateInfoi[istate].leftStateInfo->rightStateInfo->quantaStates;w2shape[2] = rOverlap.dshape()[1];
    inter1[0] = fOverlap.dshape()[1];inter1[1] = stateInfoi[istate].leftStateInfo->rightStateInfo->quantaStates;inter1[2] = stateInfoi[istate].rightStateInfo->quantaStates;
    btas::STArray<double, 3> w1Tensor(w1shape, false), intermediate(inter1, false), w2Tensor(w2shape, false);


    w1.UnCollectQuantaAlongRows(*stateInfoi[istate].leftStateInfo, *stateInfoi[istate].rightStateInfo, w1Tensor);
    btas::SDcontract(1.0, fOverlap, btas::shape(0), w1Tensor, btas::shape(0), 0.0, intermediate);
    btas::SDcontract(1.0, intermediate, btas::shape(2), rOverlap, btas::shape(0), 0.0, w2Tensor);

    lowerStates[istate].Clear();
    lowerStates[istate].CollectQuantaAlongColumns(leftState, *rightState.unCollectedStateInfo, w2Tensor, w1.get_deltaQuantum());

  }

}

//during the opimization of a wavefunction you will continuously keep updating the overlap of the current state with all the lower lying states
void SpinAdapted::Sweep::saveUpdatedLocalOverlapMatrix(int currentState, const std::vector<int>& sites, StateInfo& leftState, StateInfo& rightState)
{
  for (int istate = 0; istate<currentState; istate++) {
    
    btas::STArray<double, 2> Overlap;

    if (sites[0] == 0) {
      std::vector<int> prevSites = sites; prevSites.pop_back();
      LoadOverlapTensor(prevSites, Overlap, istate, currentState);
    }
    else {
      std::vector<int> prevSites((++sites.begin()), sites.end());
      LoadOverlapTensor(prevSites, Overlap, istate, currentState);
    }


    std::vector<Matrix> rotation2, rotation1;
    LoadRotationMatrix(sites, rotation2, currentState);
    LoadRotationMatrix(sites, rotation1, istate);

    StateInfo stateInfoi;
    {
      Wavefunction w1;
      w1.LoadWavefunctionInfo(stateInfoi, sites, istate);
    }
    StateInfo leftStatei = *stateInfoi.leftStateInfo;
    StateInfo renormState1, renormState2;
    SpinAdapted::StateInfo::transform_state(rotation1, *stateInfoi.leftStateInfo, renormState1);
    SpinAdapted::StateInfo::transform_state(rotation2, leftState, renormState2);

    btas::TVector<btas::Dshapes, 2> overlapShape;
    btas::TVector<btas::Dshapes, 3> shape1, shape2, interShape;
    shape2[0] = leftState.leftStateInfo->quantaStates; shape2[1] = leftState.rightStateInfo->quantaStates; shape2[2] = renormState2.quantaStates;
    shape1[0] = leftStatei.leftStateInfo->quantaStates; shape1[1] = leftStatei.rightStateInfo->quantaStates; shape1[2] = renormState1.quantaStates;

    overlapShape[0] = renormState1.quantaStates; overlapShape[1] = renormState2.quantaStates;
    interShape[0] = leftState.leftStateInfo->quantaStates; interShape[1] = leftStatei.rightStateInfo->quantaStates; interShape[2] = renormState1.quantaStates;


    btas::STArray<double, 3> SiteTensor1(shape1,false), 
      SiteTensor2(shape2,false), 
      intermediate(interShape, false);
    btas::STArray<double, 2>  newOverlap(overlapShape,false); //O(ni, mi) matrix

    UnCollectQuantaAlongRows(leftState, rightState, rotation2, SiteTensor2);
    UnCollectQuantaAlongRows(*stateInfoi.leftStateInfo, renormState1, rotation1, SiteTensor1);
    
    btas::SDcontract(1.0, Overlap, btas::shape(0), SiteTensor1, btas::shape(0), 0.0, intermediate);
    btas::SDcontract(1.0, intermediate, btas::shape(0,1), SiteTensor2, btas::shape(0,1), 0.0, newOverlap);

    for (int i=0; i<renormState1.quanta.size(); i++)
      for (int k=0; k<renormState2.quanta.size(); k++) 
	if (renormState2.quanta[k] != renormState1.quanta[i]) 
	  newOverlap.erase(btas::make_array(i,k));

    SaveOverlapTensor(sites, newOverlap, istate, currentState);

  }    
}


//spinBlocks also have overlap stored in them. This overlap is updated everytime the spinblock is 
//renormalized. The overlap is usually an identity matrix unless the bra and the ket states are
//different
boost::shared_ptr<SparseMatrix> SpinAdapted::Sweep::updateLocalOverlapMatrix(const SparseMatrix &Overlap, const std::vector<Matrix>& leftrotateMatrix, const std::vector<Matrix>& rightrotateMatrix,  const StateInfo& braState, const StateInfo& ketState)
{
  btas::TVector<btas::Dshapes, 2> overlapShape, newOverlapShape;
  overlapShape[0] = braState.leftStateInfo->quantaStates; overlapShape[1] = ketState.leftStateInfo->quantaStates;
  btas::STArray<double, 2> OverlapBtas(overlapShape, false);

  makeBTASOverlapFromDMRGOverlap(Overlap, OverlapBtas);

  
  StateInfo braRenormState, ketRenormState;
  SpinAdapted::StateInfo::transform_state(leftrotateMatrix, braState, braRenormState);
  SpinAdapted::StateInfo::transform_state(rightrotateMatrix, ketState, ketRenormState);
  
  btas::TVector<btas::Dshapes, 3> shapeBra, shapeKet, interShape;
  shapeBra[0] = braState.leftStateInfo->quantaStates; shapeBra[1] = braState.rightStateInfo->quantaStates; shapeBra[2] = braRenormState.quantaStates;
  shapeKet[0] = ketState.leftStateInfo->quantaStates; shapeKet[1] = ketState.rightStateInfo->quantaStates; shapeKet[2] = ketRenormState.quantaStates;
  interShape[0] = ketState.leftStateInfo->quantaStates; interShape[1] = ketState.rightStateInfo->quantaStates; interShape[2] = braRenormState.quantaStates;
  newOverlapShape[0] = braRenormState.quantaStates; newOverlapShape[1] = ketRenormState.quantaStates;


  btas::STArray<double, 3> BraSiteTensor(shapeBra,false), 
    KetSiteTensor(shapeKet,false), 
    intermediate(interShape, false);
  btas::STArray<double, 2>  newOverlap(newOverlapShape,false); //O(ni, mi) matrix

  UnCollectQuantaAlongRows(ketState, ketRenormState, rightrotateMatrix, KetSiteTensor);
  UnCollectQuantaAlongRows(braState, braRenormState, leftrotateMatrix, BraSiteTensor);
  
  btas::SDcontract(1.0, OverlapBtas, btas::shape(0), BraSiteTensor, btas::shape(0), 0.0, intermediate);
  btas::SDcontract(1.0, intermediate, btas::shape(0,1), KetSiteTensor, btas::shape(0,1), 0.0, newOverlap);

  boost::shared_ptr<SparseMatrix> newOverlapOp = boost::shared_ptr<SparseMatrix>(new Ham);
  newOverlapOp->set_initialised() = true;
  newOverlapOp->set_fermion() = false;
  newOverlapOp->set_built() = true;
  newOverlapOp->set_deltaQuantum(1, SpinQuantum(0,SpinSpace(0), IrrepSpace(0)));
  newOverlapOp->allocate(braRenormState, ketRenormState);

  for (int i=0; i<braRenormState.quanta.size(); i++)
    for (int k=0; k<ketRenormState.quanta.size(); k++) 
      if (braRenormState.quanta[i] == ketRenormState.quanta[k]) 
	DCOPY(newOverlapOp->operator_element(i,k).Storage(), &newOverlap.find(btas::make_array(i,k))->second->operator()(0, 0), 1, newOverlapOp->operator_element(i,k).Store(), 1);
  return newOverlapOp;

}



//Canonicalize wavefunction, takes the wavefunction and does a sweep to update all the roation matrices so that we get a consistent wavefunction along the whole sweep
void SpinAdapted::Sweep::CanonicalizeWavefunction(SweepParams &sweepParams, const bool &forward, int currentstate)
{

  sweepParams.set_sweep_parameters();
  sweepParams.set_block_iter() = 0;

  std::vector<int> sites;
  int new_site, wave_site;
  if (forward) {
    pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in forwards direction"<<endl;
    new_site = 0;
  }
  else {
    pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in backwards direction" << endl;
    new_site = dmrginp.last_site()-1;
  }
  pout << "\t\t\t ============================================================================ " << 
    endl;
  pout << new_site<<endl;
  sites.push_back(new_site);
    
    
  //only need statinfos
  StateInfo stateInfo1; makeStateInfo(stateInfo1, new_site);
  
  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); ) {
      
    pout << "\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
    pout << "\t\t\t ----------------------------" << endl;
    
    if (forward) {
      new_site++;
      wave_site = new_site+1;
      pout << "\t\t\t Current direction is :: Forwards " << endl;
    }
    else {
      new_site--;
      wave_site = new_site-1;
      pout << "\t\t\t Current direction is :: Backwards " << endl;
    }
    std::vector<int> complementarySites, spindotsites(1, new_site), oldsites = sites, oldcomplement;
    sites.push_back(new_site);
    getComplementarySites(sites, complementarySites);
    getComplementarySites(oldsites, oldcomplement);
    
    StateInfo siteState, newState1, bigstate, envstate; 
    makeStateInfo(siteState, new_site);
    TensorProduct(stateInfo1, siteState, newState1, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
    newState1.CollectQuanta();

    Wavefunction w; w.set_deltaQuantum() = dmrginp.effective_molecule_quantum_vec();
    w.set_onedot(true);
    StateInfo::restore(!forward, complementarySites, envstate, currentstate);

    TensorProduct(newState1, envstate, bigstate, PARTICLE_SPIN_NUMBER_CONSTRAINT);

    if (sweepParams.get_block_iter() == 0) 
      GuessWave::transpose_previous_wavefunction(w, bigstate, complementarySites, spindotsites, currentstate, true, true);
    else 
      GuessWave::transform_previous_wavefunction(w, bigstate, oldsites, oldcomplement, currentstate, true, true);
    
    w.SaveWavefunctionInfo(bigstate, sites, currentstate);

      
    //make the newstate
    std::vector<Matrix> rotation1; 
      
      
    DensityMatrix tracedMatrix;
    tracedMatrix.allocate(*bigstate.leftStateInfo);
    operatorfunctions::MultiplyProduct(w, Transpose(const_cast<Wavefunction&> (w)), tracedMatrix, 1.0);
    if (!mpigetrank())
      double error = makeRotateMatrix(tracedMatrix, rotation1, envstate.totalStates, sweepParams.get_keep_qstates());
    SaveRotationMatrix (sites, rotation1, currentstate);
    
    StateInfo renormState1;
    SpinAdapted::StateInfo::transform_state(rotation1, newState1, renormState1);
    StateInfo::store(forward, sites, renormState1, currentstate);
    stateInfo1 = renormState1;
    ++sweepParams.set_block_iter();
  }
  
}


//Initialize stateinfo using the rotation matrices
void SpinAdapted::Sweep::InitializeStateInfo(SweepParams &sweepParams, const bool &forward, int currentstate)
{
  sweepParams.set_sweep_parameters();
  sweepParams.set_block_iter() = 0;

  std::vector<int> sites;
  int new_site, wave_site;
  if (forward) 
    new_site = 0;
  else 
    new_site = dmrginp.last_site()-1;
  
  sites.push_back(new_site);
    
    
  //only need statinfos
  StateInfo stateInfo1; makeStateInfo(stateInfo1, new_site);
  StateInfo::store(forward, sites, stateInfo1, currentstate);
  
  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); ) {
      
    if (forward) 
      new_site++;
    else 
      new_site--;

    sites.push_back(new_site);
    
    StateInfo siteState, newState1; 
    makeStateInfo(siteState, new_site);
    TensorProduct(stateInfo1, siteState, newState1, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
    newState1.CollectQuanta();

      
    //make the newstate
    std::vector<Matrix> rotation1; 
      
    LoadRotationMatrix (sites, rotation1, currentstate);
    
    StateInfo renormState1;
    SpinAdapted::StateInfo::transform_state(rotation1, newState1, renormState1);
    StateInfo::store(forward, sites, renormState1, currentstate);
    stateInfo1 = renormState1;
    ++sweepParams.set_block_iter();
  }
  
}


//before you start optimizing each state you want to initalize all the overlap matrices
void Sweep::InitializeAllOverlaps(SweepParams &sweepParams, const bool &forward, int stateA, int stateB)
{
  int lowerState = stateA < stateB ? stateA : stateB ;
  int higherState = stateA < stateB ? stateB : stateA ;

  sweepParams.set_sweep_parameters();
  sweepParams.set_block_iter() = 0;
  
  std::vector<int> sites;
  int new_site, wave_site;
  if (forward) {
    pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in forwards direction"<<endl;
    new_site = 0;
  }
  else {
    pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in backwards direction" << endl;
    new_site = dmrginp.last_site()-1;
  }
  pout << "\t\t\t ============================================================================ " << 
    endl;
  sites.push_back(new_site);
  
  
  //only need statinfos
  StateInfo stateInfo1; makeStateInfo(stateInfo1, new_site);
  StateInfo stateInfo2; makeStateInfo(stateInfo2, new_site);
  
  
  
  btas::TVector<btas::Dshapes, 2> overlapShape; overlapShape[0]=stateInfo1.quantaStates; overlapShape[1] = stateInfo2.quantaStates;
  btas::STArray<double, 2> Overlap(overlapShape, false); //O(ni, mi) matrix
  
  for (int i=0; i<stateInfo1.quanta.size(); i++) {
    Overlap.reserve(btas::make_array(i,i));
    Overlap.find(btas::make_array(i,i))->second->operator()(0,0) = 1;
  }
  SaveOverlapTensor(sites, Overlap, lowerState, higherState);
  
  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); ) {
    
    pout << "\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
    pout << "\t\t\t ----------------------------" << endl;
    
    if (forward) {
      new_site++;
      wave_site = new_site+1;
      pout << "\t\t\t Current direction is :: Forwards " << endl;
    }
    else {
      new_site--;
      wave_site = new_site-1;
      pout << "\t\t\t Current direction is :: Backwards " << endl;
    }
    sites.push_back(new_site);
    StateInfo siteState, newState1, newState2; makeStateInfo(siteState, new_site);
    
    //make the newstate
    TensorProduct(stateInfo1, siteState, newState1, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
    newState1.CollectQuanta();
    
    TensorProduct(stateInfo2, siteState, newState2, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
    newState2.CollectQuanta();
    
    
    std::vector<Matrix> rotation1, rotation2; 
    LoadRotationMatrix(sites, rotation1, lowerState);
    LoadRotationMatrix(sites, rotation2, higherState);
    
    StateInfo renormState1, renormState2;
    SpinAdapted::StateInfo::transform_state(rotation1, newState1, renormState1);
    SpinAdapted::StateInfo::transform_state(rotation2, newState2, renormState2);
    
    btas::TVector<btas::Dshapes, 3> shape1, shape2, interShape;
    shape1[0] = stateInfo1.quantaStates; shape1[1] = siteState.quantaStates; shape1[2] = renormState1.quantaStates;
    shape2[0] = stateInfo2.quantaStates; shape2[1] = siteState.quantaStates; shape2[2] = renormState2.quantaStates;
    overlapShape[0] = renormState1.quantaStates; overlapShape[1] = renormState2.quantaStates;
    interShape[0] = stateInfo2.quantaStates; interShape[1] = siteState.quantaStates; interShape[2] = renormState1.quantaStates;
    
    
    btas::STArray<double, 3> SiteTensor1(shape1,false), 
      SiteTensor2(shape2,false), 
      intermediate(interShape, false);
    btas::STArray<double, 2>  newOverlap(overlapShape,false); //O(ni, mi) matrix
    
    UnCollectQuantaAlongRows(newState1, renormState1, rotation1, SiteTensor1);
    UnCollectQuantaAlongRows(newState2, renormState2, rotation2, SiteTensor2);
    
    btas::SDcontract(1.0, Overlap, btas::shape(0), SiteTensor1, btas::shape(0), 0.0, intermediate);
    btas::SDcontract(1.0, intermediate, btas::shape(0,1), SiteTensor2, btas::shape(0,1), 0.0, newOverlap);
    
    
    for (int i=0; i<renormState1.quanta.size(); i++)
      for (int k=0; k<renormState2.quanta.size(); k++) 
	if (renormState2.quanta[k] != renormState1.quanta[i]) 
	  newOverlap.erase(btas::make_array(i,k));
    
    SaveOverlapTensor(sites, newOverlap, lowerState, higherState);
    
    
    stateInfo1 = renormState1;
    stateInfo2 = renormState2;
    Overlap.copy(newOverlap);
    
    
    ++sweepParams.set_block_iter();
  }
  
  
}

//generate all the overlap elements between all the mps wavefunction stored on disk. 
void SpinAdapted::Sweep::calculateAllOverlap(Matrix& O)
{
  bool direction;
  int restartsize;
  SweepParams sweepParams;
  sweepParams.restorestate(direction, restartsize);
  sweepParams.current_root() = 0;

  std::vector<int> sites(dmrginp.last_site()-2,0), wavesites(dmrginp.last_site()-1,0), quanta(3,1);
  for (int i=0; i<dmrginp.last_site()-2; i++)
    sites[i] = i;
  for (int i=0; i<dmrginp.last_site()-1; i++)
    wavesites[i] = i;

  Wavefunction w1, w2;
  StateInfo statew1, statew2;
  btas::TVector<btas::Dshapes, 2> overlapShape; overlapShape[0]=quanta; overlapShape[1] = quanta;
  btas::STArray<double, 2>  output(overlapShape, false); //O(ni, mi) matrix
  O.ReSize(dmrginp.nroots(), dmrginp.nroots()); O = 0.0;

  for (int i=0; i<dmrginp.nroots(); i++) {
    for (int j=i; j<dmrginp.nroots(); j++) { 
      Sweep::InitializeAllOverlaps(sweepParams, true, i, j);

      w1.LoadWavefunctionInfo(statew1, wavesites, i);
      w2.LoadWavefunctionInfo(statew2, wavesites, j);

      btas::TVector<btas::Dshapes,3> w1shape, w2shape, inter;
      w1shape[0] = statew1.leftStateInfo->leftStateInfo->quantaStates;w1shape[1] = statew1.leftStateInfo->rightStateInfo->quantaStates;w1shape[2] = statew1.rightStateInfo->quantaStates;
      w2shape[0] = statew2.leftStateInfo->leftStateInfo->quantaStates;w2shape[1] = statew2.leftStateInfo->rightStateInfo->quantaStates;w2shape[2] = statew2.rightStateInfo->quantaStates;
      inter[0] = statew2.leftStateInfo->leftStateInfo->quantaStates;inter[1] = statew1.leftStateInfo->rightStateInfo->quantaStates;inter[2] = statew1.rightStateInfo->quantaStates;
      overlapShape[0] = statew1.leftStateInfo->leftStateInfo->quantaStates;overlapShape[1] = statew2.leftStateInfo->leftStateInfo->quantaStates;

      btas::STArray<double, 3> w1Tensor(w1shape, false), w2Tensor(w2shape, false), intermediate(inter, false);
      btas::STArray<double, 2> Overlap(overlapShape, false);
      LoadOverlapTensor(sites, Overlap, i, j);

      w1.UnCollectQuantaAlongRows(*statew1.leftStateInfo, *statew1.rightStateInfo, w1Tensor);
      w2.UnCollectQuantaAlongRows(*statew2.leftStateInfo, *statew2.rightStateInfo, w2Tensor);
      
      btas::SDcontract(1.0, Overlap, btas::shape(0), w1Tensor, btas::shape(0), 0.0, intermediate);
      btas::SDcontract(1.0, intermediate, btas::shape(0,1), w2Tensor, btas::shape(0,1), 0.0, output);

      //cout << output<<endl;
      double o = 0.0;
      for (int k=0; k<3; k++)
	o += output.find(btas::make_array(k,k))->second->operator()(0,0);
      cout << "overlap between "<<i<<"  and "<<j<<"  =  "<<o<<endl;
      O(i+1, j+1) = o;

    }
  }
}

void SpinAdapted::Sweep::calculateHMatrixElements(Matrix& H)
{
  bool direction, restoredirection;
  int restartsize;
  SweepParams sweepParams;
  sweepParams.restorestate(restoredirection, restartsize);
  sweepParams.current_root() = 0;

  direction = restoredirection;

  int sysdotsite, envsite, envsize=2;
  std::vector<int> sites(dmrginp.last_site()-envsize-1,0), wavesites(dmrginp.last_site()-envsize,0), envsites(envsize,0);
  if (direction) {
    sysdotsite = dmrginp.last_site()-envsize-1;
    envsite = dmrginp.last_site()-envsize;
    for (int i=0; i<dmrginp.last_site()-envsize-1; i++)
      sites[i] = i;
    for (int i=0; i<dmrginp.last_site()-envsize; i++)
      wavesites[i] = i;
    for (int i=dmrginp.last_site()-envsize; i<dmrginp.last_site(); i++) 
      envsites[-dmrginp.last_site()+envsize+i] = i;
  }
  else {
    sysdotsite = envsize;
    envsite = 0;
    for (int i=0; i<dmrginp.last_site()-envsize-1; i++)
      sites[i] = envsize+1+i;
    for (int i=0; i<dmrginp.last_site()-envsize; i++)
      wavesites[i] = envsize+i;
    for (int i=0; i<envsize; i++)
      envsites[i] = i;
  }

  pout.precision(12);
  Wavefunction w1, w2;
  StateInfo statew1, statew2;
  H.ReSize(dmrginp.nroots(), dmrginp.nroots()); H = 0.0;

  for (int i=0; i<dmrginp.nroots(); i++) {
    for (int j=i; j<dmrginp.nroots(); j++) { 
      SpinAdapted::SweepGenblock::do_one(sweepParams, direction, i, j);

      w1.LoadWavefunctionInfo(statew1, wavesites, i);
      w2.LoadWavefunctionInfo(statew2, wavesites, j);


      SpinBlock newSystem, system, systemDot, env, big;
      SpinBlock::restore(direction, sites, system, i, j);
      systemDot = SpinBlock(sysdotsite, sysdotsite, i==j);
      system.addAdditionalCompOps();

      InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, i, j, 1, direction, 
				     DISTRIBUTED_STORAGE, false, true);


      if (envsize == 1) {
	int dotsize = 1, restartsize = 0;
	int restart = false, warmup = false;
	InitBlocks::InitStartingBlock(env, !direction, i, j, dotsize, dotsize, restartsize, restart, warmup);
      }
      else 
	SpinBlock::restore(!direction, envsites, env, i, j);

      pout << env<<endl;
      pout << env.get_ketStateInfo()<<endl;

      newSystem.set_loopblock(false);
      system.set_loopblock(false);
      env.set_loopblock(true);
      
      InitBlocks::InitBigBlock(newSystem, env, big); 

      w1.LoadWavefunctionInfo(statew1, wavesites, i);
      w2.LoadWavefunctionInfo(statew2, wavesites, j);

      pout << w1<<endl;
      pout << w2<<endl;

#ifndef SERIAL
      mpi::communicator world;
      broadcast(world, w1, 0);
      broadcast(world, w2, 0);
#endif
      Wavefunction Hw2 = w1; Hw2.Clear();

      big.multiplyH(w2, &Hw2, MAX_THRD);
      double o = DotProduct(w1, Hw2);

      cout <<"matrix element "<<o<<endl;
      if (i==j) o += dmrginp.get_coreenergy();
      H(i+1, j+1) = o;
      H(j+1, i+1) = o;
      pout << H <<endl; 
    }
  }
}

/*
void SpinAdapted::Sweep::calculateHMatrixElements(Matrix& H)
{
  bool direction;
  int restartsize;
  SweepParams sweepParams;
  sweepParams.restorestate(direction, restartsize);
  sweepParams.current_root() = 0;

  /*
  std::vector<int> sites(dmrginp.last_site()-2,0), wavesites(dmrginp.last_site()-1,0), quanta(3,1);
  for (int i=0; i<dmrginp.last_site()-2; i++)
    sites[i] = i;
  for (int i=0; i<dmrginp.last_site()-1; i++)
    wavesites[i] = i;

  pout.precision(12);
  Wavefunction w1, w2;
  StateInfo statew1, statew2;
  H.ReSize(dmrginp.nroots(), dmrginp.nroots()); H = 0.0;
  for (int i=0; i<dmrginp.nroots(); i++) {
    for (int j=i; j<dmrginp.nroots(); j++) { 
      SpinAdapted::SweepGenblock::do_one(sweepParams, true, i, j);
      SpinAdapted::SweepGenblock::do_one(sweepParams, false, i, j);

      std::vector<int> sites(1,0), wavesites(1,0); wavesites.push_back(1);
      for (int k=0; k<5; k++) {
	w1.LoadWavefunctionInfo(statew1, wavesites, i);
	w2.LoadWavefunctionInfo(statew2, wavesites, j);

	Wavefunction Hw1 = w2; Hw1.Clear();
	SpinBlock newSystem, system, systemDot, env, big;
	SpinBlock::restore(true, sites, system, i, j);
	systemDot = SpinBlock(k+1, k+1, i==j);
	
	InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, i, j, 1, true, 
				       DISTRIBUTED_STORAGE, true, true);
	
	SpinBlock::restore(false, newSystem.get_complementary_sites(), env, i, j);
	//cout << env<<endl;
	//env.printOperatorSummary();

	newSystem.set_loopblock(true);
	system.set_loopblock(false);
	env.set_loopblock(false);
	InitBlocks::InitBigBlock(newSystem, env, big); 
	
	//cout << *env.getOverlap()<<endl;
	cout << DotProduct(w1, w2)<<endl;
	cout << DotProduct(w1, w1)<<endl;
	cout << DotProduct(w2, w2)<<endl;
	cout << DotProduct(Hw1, Hw1)<<endl;
	
	big.multiplyH(w1, &Hw1, MAX_THRD);
	cout <<i<<"  "<<j<<" ****************MAT ELEMENT "<< DotProduct(w2, Hw1)<<endl;
	cout << DotProduct(Hw1, Hw1)<<endl;
	double o = DotProduct(w2, Hw1);
	
	sites.push_back(k+1);
	wavesites.push_back(k+2);
	H(i+1, j+1) = 0;
	H(j+1, i+1) = 0;
      }
      cout << "Matrix element between "<<i<<" and  "<<j<<" = "<<H(i+1, j+1)<<endl;
    }
  }
}
*/

void SpinAdapted::Wavefunction::CollectQuantaAlongRows (const StateInfo& sRow, const StateInfo& sCol, btas::STArray<double, 3>& siteWave, const vector<SpinQuantum> dQ)
{
  //mdebugcheck("before collectquantaalongrows");
  try
    {
      Timer ctimer;
      deltaQuantum = dQ;
      //first make the unpacked wavefunction
      AllowQuantaFor (sRow, sCol, deltaQuantum);

      for (int i=0; i<sRow.quanta.size(); i++)
	    for (int k=0; k<sCol.quanta.size(); k++) {
	      int ileft = sRow.leftUnMapQuanta[i], iright = sRow.rightUnMapQuanta[i];
          for (int j = 0; j < deltaQuantum.size(); ++j) {
	        if (deltaQuantum[j].allow(sCol.quanta[k], sRow.quanta[i])) {
	          if( siteWave.find(btas::make_array(ileft,iright,k)) == siteWave.end()) continue;
	          DCOPY(operator_element(i,k).Storage(), &siteWave.find(btas::make_array(ileft,iright,k))->second->operator()(0, 0, 0), 1, operator_element(i,k).Store(), 1);
	        }
          }
	    }
      
      
      //now collect the wavefunction
      StateInfo tmpState = sRow;
      tmpState.CollectQuanta ();
      Wavefunction tmpOper;

      tmpOper.AllowQuantaFor (tmpState, sCol, deltaQuantum);
      ObjectMatrix<Matrix*> matRef;
      for (int i = 0; i < tmpState.quanta.size (); ++i)
	for (int j = 0; j < sCol.quanta.size (); ++j)
	  {
	    std::vector<int> dum (1); dum [0] = j;
	    if (tmpOper.allowed(i, j))
	      {
		OperatorMatrixReference (matRef, tmpState.oldToNewState [i], dum);
		CatenateProduct (matRef, tmpOper.operator_element(i,j), true);
	      }
	  }
      allowedQuantaMatrix = tmpOper.allowedQuantaMatrix;
      copy(tmpOper.operatorMatrix, operatorMatrix);
    }
  catch (Exception)
    {
      Exception::what ();
      abort ();
    }
  //mdebugcheck("after collectquantaalongrows");
}


void SpinAdapted::Wavefunction::UnCollectQuantaAlongRows (const StateInfo& sRow, const StateInfo& sCol, btas::STArray<double, 3>& siteWave)
{
  try
    {
      Wavefunction tmpOper;
      tmpOper.AllowQuantaFor (*sRow.unCollectedStateInfo, sCol, deltaQuantum);
      for (int i = 0; i < sRow.quanta.size (); ++i)
	{
	  const std::vector<int>& oldToNewStateI = sRow.oldToNewState [i];
	  int firstRow = 0;
	  for (int iSub = 0; iSub < oldToNewStateI.size (); ++iSub)
	    {
	      int unCollectedI = oldToNewStateI [iSub];
	      int lastRowSize = sRow.unCollectedStateInfo->quantaStates [unCollectedI];
	      for (int j = 0; j < sCol.quanta.size (); ++j)
		if (tmpOper.allowedQuantaMatrix (unCollectedI, j))
		  tmpOper.operator_element(unCollectedI, j) = operator_element(i, j).Rows (firstRow + 1, firstRow + lastRowSize);
	      firstRow += lastRowSize;
	    }
	}

      StateInfo* lState =sRow.leftStateInfo, *rState=sRow.rightStateInfo; 
      boost::shared_ptr<StateInfo> unCollected = sRow.unCollectedStateInfo;
      
      for (int i=0; i<unCollected->quanta.size(); i++)
	for (int k=0; k<sCol.quanta.size(); k++) {
	  int ileft = sRow.leftUnMapQuanta[i], iright = sRow.rightUnMapQuanta[i];
	  for (int j = 0; j < deltaQuantum.size(); ++j) {
        if (deltaQuantum[j].allow(sCol.quanta[k], unCollected->quanta[i])) {
	      siteWave.reserve(btas::make_array(ileft,iright,k));
	    
	      DCOPY(tmpOper.operator_element(i,k).Storage(), tmpOper.operator_element(i,k).Store(), 1, &siteWave.find(btas::make_array(ileft,iright,k))->second->operator()(0, 0, 0), 1);
	    }
      }
	}
      
      
    }
  catch (Exception)
    {
      Exception::what ();
      abort ();
    }
}

void SpinAdapted::Wavefunction::CollectQuantaAlongColumns (const StateInfo& sRow, const StateInfo& sCol, btas::STArray<double, 3>& siteWave, const vector<SpinQuantum> dQ)
{
  //mdebugcheck("before collectquantaalongrows");
  try
    {
      Timer ctimer;
      deltaQuantum = dQ;
      //first make the unpacked wavefunction
      AllowQuantaFor (sRow, sCol, deltaQuantum);

      for (int i=0; i<sRow.quanta.size(); i++)
	for (int k=0; k<sCol.quanta.size(); k++) {
	  int kleft = sCol.leftUnMapQuanta[k], kright = sCol.rightUnMapQuanta[k];
	  for (int j = 0; j < deltaQuantum.size(); ++j) {
        if (deltaQuantum[j].allow(sCol.quanta[k], sRow.quanta[i])) {


	      if( siteWave.find(btas::make_array(i,kright,kleft)) == siteWave.end()) continue;
	      //s.e --> se.
	      SpinQuantum Bq = sCol.leftStateInfo->quanta[kleft];
	      SpinQuantum Cq = sCol.rightStateInfo->quanta[kright];
	      SpinQuantum CBq = sCol.quanta[k];
	      int parity1 = getCommuteParity(Cq, Bq, CBq);		
	      DSCAL(operator_element(i,k).Storage(), parity1*1.0, &siteWave.find(btas::make_array(i,kright,kleft))->second->operator()(0, 0, 0), 1); 
	      
	      DCOPY(operator_element(i,k).Storage(), &siteWave.find(btas::make_array(i,kright,kleft))->second->operator()(0, 0, 0), 1, operator_element(i,k).Store(), 1);

	    }
      }
	}
      
      
      //now collect the wavefunction
      StateInfo tmpState = sCol;
      tmpState.CollectQuanta ();
      Wavefunction tmpOper;

      tmpOper.AllowQuantaFor (sRow, tmpState, deltaQuantum);
      ObjectMatrix<Matrix*> matRef;
      for (int i = 0; i < sRow.quanta.size (); ++i)
	for (int j = 0; j < tmpState.quanta.size (); ++j)
	  {
	    std::vector<int> dum (1); dum [0] = i;
	    if (tmpOper.allowed(i, j))
	      {
		OperatorMatrixReference (matRef, dum, tmpState.oldToNewState [j]);
		CatenateProduct (matRef, tmpOper.operator_element(i,j));
	      }
	  }
      allowedQuantaMatrix = tmpOper.allowedQuantaMatrix;
      copy(tmpOper.operatorMatrix, operatorMatrix);
    }
  catch (Exception)
    {
      Exception::what ();
      abort ();
    }
  //mdebugcheck("after collectquantaalongrows");
}


void SpinAdapted::Wavefunction::UnCollectQuantaAlongColumns (const StateInfo& sRow, const StateInfo& sCol, btas::STArray<double, 3>& siteWave)
{
  try
    {
      Wavefunction tmpOper;
      tmpOper.AllowQuantaFor (sRow, *sCol.unCollectedStateInfo, deltaQuantum);
      for (int i = 0; i < sCol.quanta.size (); ++i)
	{
	  const std::vector<int>& oldToNewStateI = sCol.oldToNewState [i];
	  int firstCol = 0;
	  for (int iSub = 0; iSub < oldToNewStateI.size (); ++iSub)
	    {
	      int unCollectedI = oldToNewStateI [iSub];
	      int lastColSize = sCol.unCollectedStateInfo->quantaStates [unCollectedI];
	      for (int j = 0; j < sRow.quanta.size (); ++j)
		if (tmpOper.allowed(j, unCollectedI))
		  tmpOper.operator_element(j, unCollectedI) = operator_element(i, j).Columns (firstCol + 1, firstCol + lastColSize);
	      firstCol += lastColSize;
	    }
	}

      StateInfo* lState =sCol.leftStateInfo, *rState=sCol.rightStateInfo; 
      boost::shared_ptr<StateInfo> unCollected = sCol.unCollectedStateInfo;
      
      for (int i=0; i<sRow.quanta.size(); i++) {
	for (int k=0; k<unCollected->quanta.size(); k++) {
	  int kleft = sCol.leftUnMapQuanta[i], kright = sCol.rightUnMapQuanta[i];
	  for (int j = 0; j < deltaQuantum.size(); ++j) {
        if (deltaQuantum[j].allow(sRow.quanta[i], unCollected->quanta[k])) {
	      siteWave.reserve(btas::make_array(i, kleft, kright));
	      
	      DCOPY(tmpOper.operator_element(i,k).Storage(), tmpOper.operator_element(i,k).Store(), 1, &siteWave.find(btas::make_array(i, kleft, kright))->second->operator()(0, 0, 0), 1);
	    }
      }
	}
      }
      
    }
  catch (Exception)
    {
      Exception::what ();
      abort ();
    }
}

void SpinAdapted::UnCollectQuantaAlongRows(const StateInfo& sRow, const StateInfo& sCol, const std::vector<Matrix> &inRotation, btas::STArray<double, 3>& SiteTensor)
{
  try
  {
    std::vector<Matrix> outMatrix;
    allocate(*sRow.unCollectedStateInfo, sCol, outMatrix);

    for (int i = 0; i < sRow.quanta.size (); ++i)
    {
      const std::vector<int>& oldToNewStateI = sRow.oldToNewState [i];
      int firstRow = 0;
      for (int iSub = 0; iSub < oldToNewStateI.size (); ++iSub)
      {
	int unCollectedI = oldToNewStateI [iSub];
	int lastRowSize = sRow.unCollectedStateInfo->quantaStates [unCollectedI];
	for (int j = 0; j < sCol.quanta.size (); ++j)
	  if (sRow.unCollectedStateInfo->quanta[unCollectedI] == sCol.quanta[j])
	    outMatrix[unCollectedI] = inRotation[i].Rows (firstRow + 1, firstRow + lastRowSize);
	firstRow += lastRowSize;
      }
    }

    StateInfo* lState =sRow.leftStateInfo, *rState=sRow.rightStateInfo; 
    boost::shared_ptr<StateInfo> unCollected = sRow.unCollectedStateInfo;

    for (int i=0; i<unCollected->quanta.size(); i++)
      for (int k=0; k<sCol.quanta.size(); k++) {
      int ileft = sRow.leftUnMapQuanta[i], iright = sRow.rightUnMapQuanta[i];
      if (sCol.quanta[k] == unCollected->quanta[i]) {
	SiteTensor.reserve(btas::make_array(ileft,iright,k));

	DCOPY(outMatrix[i].Storage(), outMatrix[i].Store(), 1, &SiteTensor.find(btas::make_array(ileft,iright,k))->second->operator()(0, 0, 0), 1);
      }
    }

  }
  catch (Exception)
  {
    Exception::what ();
    abort ();
  }
  
}

#else

void Sweep::InitializeAllOverlaps(SweepParams &sweepParams, const bool &forward, int stateA, int stateB) {return;}
void SpinAdapted::Sweep::CanonicalizeWavefunction(SweepParams &sweepParams, const bool &forward, int currentstate) {return;}
void SpinAdapted::Sweep::InitializeStateInfo(SweepParams &sweepParams, const bool &forward, int currentstate) {return;}
void SpinAdapted::Sweep::calculateAllOverlap(Matrix& O){return;}
void SpinAdapted::Sweep::calculateHMatrixElements(Matrix& H) {return;}

#endif
