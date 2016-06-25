/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "spinblock.h"
#include "op_components.h"
#include "operatorfunctions.h"
#include "opxop.h"
#include "wavefunction.h"
#include <boost/format.hpp>
#include "distribute.h"

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
#include "pario.h"

namespace SpinAdapted{
using namespace operatorfunctions;

void SpinBlock::printOperatorSummary()
{
#ifndef SERIAL
  mpi::communicator world;

  if (mpigetrank() != 0) {
    for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::const_iterator it = ops.begin(); it != ops.end(); ++it)
      sendobject(it->second->get_size(), 0);
  }
  else {
    for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::const_iterator it = ops.begin(); it != ops.end(); ++it)
    {
      if(it->second->is_core())
      {
         p2out << "\t\t\t " << it->second->size()<<" :  "<<it->second->get_op_string()<<"  Core Operators  ";
      }
      else
      {
         p2out << "\t\t\t " << it->second->size()<<" :  "<<it->second->get_op_string()<<"  Virtual Operators  ";      
      }
      
      vector<int> numops(world.size(), 0);
      for (int proc = 0; proc <world.size(); proc++) {
         if (proc != 0) 
            receiveobject(numops[proc],proc);
         else 
            numops[proc] = it->second->get_size();
         p2out << " " << numops[proc]<<"  ";
      }
      p2out << endl;
      /*
      if(it->second->is_core()) { 
        for (int i = 0; i < it->second->size(); ++i) {
           std::vector<boost::shared_ptr<SparseMatrix> > global_element = it->second->get_global_element(i);
           p2out << "\t\t\t Element " << i  << " has " << global_element.size() << " operators" << endl;
           for (int j = 0; j < global_element.size(); ++j) {
             p2out << "\t\t\t Operator " << j << endl; 
             p2out << "\t\t\t " << *(global_element[j]) << endl;
           }
        }
      }
      */
    }
  }
#else
  for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::const_iterator it = ops.begin(); it != ops.end(); ++it)
  {
    if(it->second->is_core()) 
      { p2out << "\t\t\t " << it->second->size()<<" :  "<<it->second->get_op_string()<<"  Core Operators  "; }
    else
      { p2out << "\t\t\t " << it->second->size()<<" :  "<<it->second->get_op_string()<<"  Virtual Operators  "; }
    p2out << endl;
  }
#endif
  
}
ostream& operator<< (ostream& os, const SpinBlock& b)
{
  os << "\t\t\t Sites ::  ";
  for (int i = 0; i < b.sites.size(); ++i) { os << b.sites[i] << " "; } 
  
  if (dmrginp.outputlevel() > 1) {
    os << endl;
    os << b.braStateInfo;
    os << b.ketStateInfo;
  }
  else {
    os <<"    # states: "<<b.braStateInfo.totalStates;
    os <<"    # states: "<<b.ketStateInfo.totalStates<<endl;
  }
  return os;
}


SpinBlock::SpinBlock () : 
  localstorage(false),
  name (rand()), 
  integralIndex(0),
  hasMemoryAllocated (false),
  loopblock(false),
  direct(false), complementary(false), normal(true), leftBlock(0), rightBlock(0) { }

SpinBlock::SpinBlock(int start, int finish, int p_integralIndex, bool implicitTranspose, bool is_complement) :  
  name (rand()), 
  hasMemoryAllocated (false), 
  integralIndex(p_integralIndex),
  direct(false), leftBlock(0), rightBlock(0),
  nonactive_orbs(0)
{
  complementary = is_complement;
  normal = !is_complement;

  //this is used to make dot block and we make the 
  //additional operators by default because they are cheap
  default_op_components(is_complement, implicitTranspose);

  std::vector<int> sites; 
  if (dmrginp.use_partial_two_integrals()) {
    if (start != finish) {
      pout << "Cannot use partial two electron integrals, when making spin block with more than two orbitals"<<endl;
      abort();
    }
    std::vector<int> o;
    for (int i=dmrginp.spatial_to_spin()[start]; i<dmrginp.spatial_to_spin()[start+1]; i+=2)
      o.push_back(i/2);
    twoInt = boost::shared_ptr<PartialTwoElectronArray> (new PartialTwoElectronArray(o));
    twoInt->Load(dmrginp.load_prefix(), integralIndex);
#ifndef SERIAL
    mpi::communicator world;
    PartialTwoElectronArray& ar = dynamic_cast<PartialTwoElectronArray&>(*twoInt.get());
    mpi::broadcast(world, ar, 0);
#endif
  }
  else
    twoInt = boost::shared_ptr<TwoElectronArray>( &v_2[integralIndex], boostutils::null_deleter());

  int lower = min(start, finish);
  int higher = max(start, finish);
  sites.resize(higher - lower + 1);
  for (int i=0; i < sites.size(); i++)
      sites[i] = lower + i;

  BuildTensorProductBlock(sites);

}

SpinBlock::SpinBlock(int start, int finish, const std::vector<int>& nonactive_orbs_, bool is_complement) :  
  name (rand()), 
  hasMemoryAllocated (false), 
  integralIndex(0),
  direct(false), leftBlock(0), rightBlock(0),
  nonactive_orbs(nonactive_orbs_)
{
  complementary = is_complement;
  normal = !is_complement;

  //this is used to make dot block and we make the 
  //additional operators by default because they are cheap
  default_op_components(is_complement, false);

  std::vector<int> sites; 
  //TODO
  //twoInt = boost::shared_ptr<TwoElectronArray>( &v_2[integralIndex], boostutils::null_deleter());
  //TODO
  //Set the integral connecting active space and nonactive space.
  twoInt = boost::shared_ptr<TwoElectronArray>( &v_2[integralIndex], boostutils::null_deleter());

  int lower = min(start, finish);
  int higher = max(start, finish);
  sites.resize(higher - lower + 1);
  for (int i=0; i < sites.size(); i++)
      sites[i] = lower + i;

  BuildTensorProductBlock(sites);
}

SpinBlock::SpinBlock (const SpinBlock& b) { *this = b; }

SpinBlock::SpinBlock(const StateInfo& s, int pintegralIndex)
{
  braStateInfo = s;
  ketStateInfo = s;
  sites.resize(0);
  integralIndex = pintegralIndex;
}



void SpinBlock::BuildTensorProductBlock(std::vector<int>& new_sites)
{

  if (twoInt.get() == 0 && dmrginp.use_partial_two_integrals()) { //this is when dummy block is being added for non zero spin
    std::vector<int> o;
    for (int i=dmrginp.spatial_to_spin()[new_sites[0]]; i<dmrginp.spatial_to_spin()[new_sites[new_sites.size()-1]+1]; i+=2)
      o.push_back(i/2);
    twoInt = boost::shared_ptr<PartialTwoElectronArray> (new PartialTwoElectronArray(o));
    twoInt->Load(dmrginp.load_prefix(), integralIndex);
#ifndef SERIAL
    mpi::communicator world;
    PartialTwoElectronArray& ar = dynamic_cast<PartialTwoElectronArray&>(*twoInt.get());
    mpi::broadcast(world, ar, 0);
    //world.broadcast(twoInt);
#endif
  }
  else if (twoInt.get() == 0)
    twoInt = boost::shared_ptr<TwoElectronArray>( &v_2[integralIndex], boostutils::null_deleter());

  name = get_name();


  std::vector< std::vector<Csf> > ladders;
  std::vector< Csf > dets; 

  if (dmrginp.spinAdapted()) {
    sites = new_sites;
    dets = CSFUTIL::spinfockstrings(new_sites, ladders);
  }
  else {
    for (int i=0; i<new_sites.size(); i++) {
      sites.push_back( dmrginp.spatial_to_spin()[new_sites[i]]   );
      sites.push_back( dmrginp.spatial_to_spin()[new_sites[i]]+1 );
    }
    // dets is   {0, alpha, beta, alphabeta;}
    // in StateInfo.C, they are sorted, they become {0,beta,alpha,alphabeta}
    dets = CSFUTIL::spinfockstrings(new_sites);
    for (int j=0; j<dets.size(); j++)
      ladders.push_back(std::vector<Csf>(1,dets[j]));
  }

  braStateInfo = StateInfo(dets);
  ketStateInfo = StateInfo(dets);

  setstoragetype(LOCAL_STORAGE);
  complementary_sites = make_complement(sites);

  //this is where we are building blocks from sites
  //currently only used for building single site blocks
  //temporarily disable screening for single site blocks
  double twoindex_ScreenTol = dmrginp.twoindex_screen_tol();
  double oneindex_ScreenTol = dmrginp.oneindex_screen_tol();
  if (new_sites.size() == 1 ){
    dmrginp.twoindex_screen_tol() = 0.0;
    dmrginp.oneindex_screen_tol() = 0.0;
  }

  build_iterators();

  if (new_sites.size() == 1 ) {
    dmrginp.twoindex_screen_tol() = twoindex_ScreenTol;
    dmrginp.oneindex_screen_tol() = oneindex_ScreenTol;
  }

  build_operators(dets, ladders);

}

std::vector<int> SpinBlock::make_complement(const std::vector<int>& sites)
{
  std::vector<int> complementary_sites;
  for (int i=0; i<dmrginp.last_site(); ++i)
    if (find(sites.begin(), sites.end(), i) == sites.end())
      complementary_sites.push_back(i);
  
  return complementary_sites;
  
}

void SpinBlock::build_iterators()
{
  dmrginp.builditeratorsT->start();
  for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
  {
    it->second->build_iterators(*this);
  }
  dmrginp.builditeratorsT->stop();
}


void SpinBlock::build_operators(std::vector< Csf >& dets, std::vector< std::vector<Csf> >& ladders)
{
  dmrginp.buildcsfops->start();
  for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
    {
      if(it->second->is_core()) {
        it->second->build_csf_operators(dets, ladders, *this);      
      }
    }
  dmrginp.buildcsfops->stop();
}
  


void SpinBlock::build_operators()
{
  for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
    {
      if(it->second->is_core()) {
	    it->second->build_operators(*this);
      }
    }
}


void SpinBlock::renormalise_transform(const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo)
{
  for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it) {
    if(it->second->is_core()) {
      it->second->renormalise_transform(rotateMatrix, stateinfo);
    }
  }
}


void SpinBlock::renormalise_transform(const std::vector<Matrix>& leftMat, const StateInfo *bra, const std::vector<Matrix>& rightMat, const StateInfo *ket)
{
  for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it) {
    if(it->second->is_core()) {
      it->second->renormalise_transform(leftMat, bra, rightMat, ket);
    }
  }
}


void SpinBlock::build_and_renormalise_operators(const std::vector<Matrix>& rotateMatrix, const StateInfo *newStateInfo)
{
  for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it) {
    opTypes ot = it->first;
    if(! it->second->is_core()) {
      it->second->build_and_renormalise_operators(*this, ot, rotateMatrix, newStateInfo);
    }
  }
}


void SpinBlock::build_and_renormalise_operators(const std::vector<Matrix>& leftMat, const StateInfo *bra, const std::vector<Matrix>& rightMat, const StateInfo *ket)
{
  for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it) {
    opTypes ot = it->first;
    if(! it->second->is_core()) {
      it->second->build_and_renormalise_operators(*this, ot, leftMat, bra, rightMat, ket);
    }
  }
}


void SpinBlock::clear()
{
  for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
    it->second->clear();
}


void SpinBlock::BuildSumBlockSkeleton(int condition, SpinBlock& lBlock, SpinBlock& rBlock, StateInfo* compState)
{

  name = get_name();
  p1out << "\t\t\t Building Sum Block " << name << endl;
  leftBlock = &lBlock;
  rightBlock = &rBlock;

  if(lBlock.nonactive_orbs.size()!=0 && rBlock.nonactive_orbs.size()==0 )
    nonactive_orbs= lBlock.nonactive_orbs;
  else if(lBlock.nonactive_orbs.size()!=0 && rBlock.nonactive_orbs.size()==0 )
    nonactive_orbs= rBlock.nonactive_orbs;
  else if(lBlock.nonactive_orbs.size()!=0 && rBlock.nonactive_orbs.size()!=0 ){
    if(lBlock.nonactive_orbs.size() != rBlock.nonactive_orbs.size()){
      pout << "Nonactive_orbs in left block and right block are different.";
      abort();
    }
    else{
      nonactive_orbs= rBlock.nonactive_orbs;
      for(int i=0; i<nonactive_orbs.size(); i++)
        if(lBlock.nonactive_orbs[i]!= rBlock.nonactive_orbs[i]){
          pout << "Nonactive_orbs in left block and right block are different.";
          abort();
        }
    }
  }
  sites.reserve (lBlock.sites.size () + rBlock.sites.size ());

  dmrginp.blockintegrals -> start();
  
  if (dmrginp.use_partial_two_integrals()) {
    if (rBlock.sites.size() == 1) {
      std::vector<int> o;
      for (int i=dmrginp.spatial_to_spin().at(rBlock.sites[0]); i<dmrginp.spatial_to_spin().at(rBlock.sites[0]+1); i+=2)
	o.push_back(i/2);
      twoInt = boost::shared_ptr<PartialTwoElectronArray> (new PartialTwoElectronArray(o));
      twoInt->Load(dmrginp.load_prefix(), integralIndex);
#ifndef SERIAL
      mpi::communicator world;
      PartialTwoElectronArray& ar = dynamic_cast<PartialTwoElectronArray&>(*twoInt.get());
      mpi::broadcast(world, ar, 0);
#endif

    }
    //pout << "Cannot use partial two electron integrals, when the dot block has more than one orbital"<<endl;
    //abort();

  }
  else 
    twoInt = boost::shared_ptr<TwoElectronArray>( &v_2[integralIndex],  boostutils::null_deleter());

  dmrginp.blockintegrals -> stop();

  dmrginp.blocksites -> start();

  sites = lBlock.sites;
  copy (rBlock.sites.begin(), rBlock.sites.end (), back_inserter (sites));
  sort(sites.begin(), sites.end());
  complementary_sites = make_complement(sites);
  p2out << "\t\t\t ";
  for (int i = 0; i < sites.size(); ++i) p2out << sites[i] << " ";
  p2out << endl;
  dmrginp.blocksites -> stop();

  dmrginp.statetensorproduct -> start();
  if(dmrginp.transition_diff_irrep()){
    if( condition== PARTICLE_SPIN_NUMBER_CONSTRAINT)
    // When bra and ket wavefuntion have different spatial or spin irrep,
    // Bra Stateinfo for the big block should not be used with quantum number of effective_molecule_quantum
      TensorProduct (lBlock.braStateInfo, rBlock.braStateInfo, dmrginp.bra_quantum(), EqualQ, braStateInfo);
    else if (condition== NO_PARTICLE_SPIN_NUMBER_CONSTRAINT)
      TensorProduct (lBlock.braStateInfo, rBlock.braStateInfo, dmrginp.bra_quantum(), LessThanQ, braStateInfo,compState);
    // When bra and ket wavefuntion have different spatial or spin irrep,
  }
  else
    TensorProduct (lBlock.braStateInfo, rBlock.braStateInfo, braStateInfo, condition, compState);

  TensorProduct (lBlock.ketStateInfo, rBlock.ketStateInfo, ketStateInfo, condition, compState);
  dmrginp.statetensorproduct -> stop();

  dmrginp.statecollectquanta -> start();
  if (!( (dmrginp.hamiltonian() == BCS && condition == SPIN_NUMBER_CONSTRAINT)  ||
	 (dmrginp.hamiltonian() != BCS && condition == PARTICLE_SPIN_NUMBER_CONSTRAINT))) {
    braStateInfo.CollectQuanta();
    ketStateInfo.CollectQuanta();
  }
  dmrginp.statecollectquanta -> stop();

}

//Build Sum Block with different quanta num in bra and ket.
void SpinBlock::BuildSumBlockSkeleton(int condition, SpinBlock& lBlock, SpinBlock& rBlock, const std::vector<SpinQuantum>& braquantum, const std::vector<SpinQuantum>& ketquantum)
{

  name = get_name();
  if (dmrginp.outputlevel() > 0) 
    pout << "\t\t\t Building Sum Block " << name << endl;
  leftBlock = &lBlock;
  rightBlock = &rBlock;

  if(lBlock.nonactive_orbs.size()!=0 && rBlock.nonactive_orbs.size()==0 )
    nonactive_orbs= lBlock.nonactive_orbs;
  else if(lBlock.nonactive_orbs.size()!=0 && rBlock.nonactive_orbs.size()==0 )
    nonactive_orbs= rBlock.nonactive_orbs;
  else if(lBlock.nonactive_orbs.size()!=0 && rBlock.nonactive_orbs.size()!=0 ){
    if(lBlock.nonactive_orbs.size() != rBlock.nonactive_orbs.size()){
      pout << "Nonactive_orbs in left block and right block are different.";
      abort();
    }
    else{
      nonactive_orbs= rBlock.nonactive_orbs;
      for(int i=0; i<nonactive_orbs.size(); i++)
        if(lBlock.nonactive_orbs[i]!= rBlock.nonactive_orbs[i]){
          pout << "Nonactive_orbs in left block and right block are different.";
          abort();
        }
    }
  }


  dmrginp.blockintegrals -> start();
  
  if (dmrginp.use_partial_two_integrals()) {
    if (rBlock.sites.size() == 1) {
      std::vector<int> o;
      for (int i=dmrginp.spatial_to_spin().at(rBlock.sites[0]); i<dmrginp.spatial_to_spin().at(rBlock.sites[0]+1); i+=2)
	o.push_back(i/2);
      twoInt = boost::shared_ptr<PartialTwoElectronArray> (new PartialTwoElectronArray(o));
      twoInt->Load(dmrginp.load_prefix(), integralIndex);
#ifndef SERIAL
      mpi::communicator world;
      PartialTwoElectronArray& ar = dynamic_cast<PartialTwoElectronArray&>(*twoInt.get());
      mpi::broadcast(world, ar, 0);
#endif

    }
    //pout << "Cannot use partial two electron integrals, when the dot block has more than one orbital"<<endl;
    //abort();

  }
  else 
    twoInt = boost::shared_ptr<TwoElectronArray>( &v_2[integralIndex],  boostutils::null_deleter());

  dmrginp.blockintegrals -> stop();

  dmrginp.blocksites -> start();

  sites.reserve (lBlock.sites.size () + rBlock.sites.size ());
  sites = lBlock.sites;
  copy (rBlock.sites.begin(), rBlock.sites.end (), back_inserter (sites));
  sort(sites.begin(), sites.end());
  complementary_sites = make_complement(sites);
  if (dmrginp.outputlevel() > 0) {
    pout << "\t\t\t ";
    for (int i = 0; i < sites.size(); ++i) pout << sites[i] << " ";
    pout << endl;
  }
  dmrginp.blocksites -> stop();

  dmrginp.statetensorproduct -> start();
  //TODO
  //It is very possible that several quantum numbers for bra.
  //For example, Perturber can have different Spatial Symmetry.
  //TODO
  //What does compState do? 
  //It is not used at all in the whole block code.
  TensorProduct (lBlock.braStateInfo, rBlock.braStateInfo, braquantum[0],condition, braStateInfo);
  TensorProduct (lBlock.ketStateInfo, rBlock.ketStateInfo, ketquantum[0],condition, ketStateInfo);
  dmrginp.statetensorproduct -> stop();

  dmrginp.statecollectquanta -> start();
  if (!( (dmrginp.hamiltonian() == BCS && condition == SPIN_NUMBER_CONSTRAINT)  ||
	 (dmrginp.hamiltonian() != BCS && condition == PARTICLE_SPIN_NUMBER_CONSTRAINT))) {
    braStateInfo.CollectQuanta();
    ketStateInfo.CollectQuanta();
  }
  dmrginp.statecollectquanta -> stop();

}

void SpinBlock::BuildSumBlock(int condition, SpinBlock& lBlock, SpinBlock& rBlock, StateInfo* compState)
{
  if (!(lBlock.integralIndex == rBlock.integralIndex && lBlock.integralIndex == integralIndex))  {
    pout << "The left, right and dot block should use the same integral indices"<<endl;
    pout << "ABORTING!!"<<endl;
    exit(0);
  }
  dmrginp.buildsumblock -> start();
  BuildSumBlockSkeleton(condition, lBlock, rBlock, compState);

  build_iterators();

  dmrginp.buildblockops -> start();
  build_operators();
  dmrginp.buildblockops -> stop();
  dmrginp.buildsumblock -> stop();
}

//Build Sum Block with different quanta num in bra and ket.
void SpinBlock::BuildSumBlock(int condition, SpinBlock& lBlock, SpinBlock& rBlock, const std::vector<SpinQuantum>& braquantum, const std::vector<SpinQuantum>& ketquantum)
{
  if (!(lBlock.integralIndex == rBlock.integralIndex && lBlock.integralIndex == integralIndex))  {
    pout << "The left, right and dot block should use the same integral indices"<<endl;
    pout << "ABORTING!!"<<endl;
    exit(0);
  }
  dmrginp.buildsumblock -> start();
  BuildSumBlockSkeleton(condition, lBlock, rBlock, braquantum,ketquantum);

  build_iterators();

  dmrginp.buildblockops -> start();
  build_operators();
  dmrginp.buildblockops -> stop();
  dmrginp.buildsumblock -> stop();
}


void SpinBlock::operator= (const SpinBlock& b)
{
  localstorage = b.localstorage;
  name = b.name;
  complementary = b.is_complementary();
  normal = b.is_normal();
  loopblock = b.is_loopblock();

  hasMemoryAllocated = b.hasMemoryAllocated;
  sites = b.sites;
  complementary_sites = b.complementary_sites;
  integralIndex = b.integralIndex;

  direct = b.is_direct();

  braStateInfo = b.braStateInfo;
  ketStateInfo = b.ketStateInfo;
  leftBlock = b.leftBlock;
  rightBlock = b.rightBlock;
  twoInt = b.twoInt;
  ops = b.ops;
}

void SpinBlock::initialise_op_array(opTypes optype, bool is_core)
{
  ops[optype] = make_new_op(optype, is_core);
  return;
}

void SpinBlock::multiplyOverlap(Wavefunction& c, Wavefunction* v, int num_threads) const
{
  if (mpigetrank() == 0) {
    boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_array(OVERLAP).get_local_element(0)[0]->getworkingrepresentation(leftBlock);
    
    boost::shared_ptr<SparseMatrix> overlap = rightBlock->get_op_array(OVERLAP).get_local_element(0)[0]->getworkingrepresentation(rightBlock);
    TensorMultiply(leftBlock, *op, *overlap, this, c, *v, op->get_deltaQuantum(0) ,1.0);  // dmrginp.ef
  }

}

void SpinBlock::multiplyCDD_sum(Wavefunction& c, Wavefunction* v, int num_threads) const
{

  SpinBlock* loopBlock=(leftBlock->is_loopblock()) ? leftBlock : rightBlock;
  SpinBlock* otherBlock = loopBlock == leftBlock ? rightBlock : leftBlock;

  Wavefunction *v_array=0, *v_distributed=0, *v_add=0;

  int maxt = 1;
  initiateMultiThread(v, v_array, v_distributed, MAX_THRD);
  dmrginp.oneelecT -> start();
  dmrginp.s0time -> start();

  //if (mpigetrank() == 0) {
    boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_array(CDD_SUM).get_local_element(0)[0]->getworkingrepresentation(leftBlock);
    boost::shared_ptr<SparseMatrix> overlap = rightBlock->get_op_array(OVERLAP).get_local_element(0)[0]->getworkingrepresentation(rightBlock);
    TensorMultiply(leftBlock, *op, *overlap, this, c, *v, op->get_deltaQuantum(0) ,1.0);  // dmrginp.effective_molecule_quantum() is never used in TensorMultiply
    
    overlap = leftBlock->get_op_array(OVERLAP).get_local_element(0)[0]->getworkingrepresentation(leftBlock);
    op = rightBlock->get_op_array(CDD_SUM).get_local_element(0)[0]->getworkingrepresentation(rightBlock);
    TensorMultiply(rightBlock, *op, *overlap, this, c, *v, op->get_deltaQuantum(0), 1.0);  
    //}

  dmrginp.s0time -> stop();
#ifndef SERIAL
  boost::mpi::communicator world;
  int size = world.size();
#endif

  SpinQuantum q= -getSpinQuantum( nonactive_orb(0));  dmrginp.s1time -> start();
  v_add =  leftBlock->get_op_array(CDD_CRE_DESCOMP).is_local() ? v_array : v_distributed;
  Functor f = boost::bind(&opxop::cdd_dxcdcomp, leftBlock, _1, this, boost::ref(c), v_add,q); 
  for_all_multithread(rightBlock->get_op_array(DES), f);

  v_add =  rightBlock->get_op_array(CDD_CRE_DESCOMP).is_local() ? v_array : v_distributed;
  f = boost::bind(&opxop::cdd_dxcdcomp, rightBlock, _1, this, boost::ref(c), v_add, q ); 
  for_all_multithread(leftBlock->get_op_array(DES), f);  
    
  v_add =  leftBlock->get_op_array(CDD_DES_DESCOMP).is_local() ? v_array : v_distributed;
  f = boost::bind(&opxop::cdd_cxddcomp, leftBlock, _1, this, boost::ref(c), v_add, q);
  for_all_multithread(rightBlock->get_op_array(CRE), f);
  
  v_add =  rightBlock->get_op_array(CDD_DES_DESCOMP).is_local() ? v_array : v_distributed;
  f = boost::bind(&opxop::cdd_cxddcomp, rightBlock, _1, this, boost::ref(c), v_add, q);
  for_all_multithread(leftBlock->get_op_array(CRE), f);
  

  accumulateMultiThread(v, v_array, v_distributed, MAX_THRD);
  dmrginp.oneelecT -> stop();

}

void SpinBlock::multiplyCCD_sum(Wavefunction& c, Wavefunction* v, int num_threads) const
{

  SpinBlock* loopBlock=(leftBlock->is_loopblock()) ? leftBlock : rightBlock;
  SpinBlock* otherBlock = loopBlock == leftBlock ? rightBlock : leftBlock;

  Wavefunction *v_array=0, *v_distributed=0, *v_add=0;

  int maxt = 1;
  initiateMultiThread(v, v_array, v_distributed, MAX_THRD);
  dmrginp.oneelecT -> start();
  dmrginp.s0time -> start();

  //if (mpigetrank() == 0) {
    boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_array(CCD_SUM).get_local_element(0)[0]->getworkingrepresentation(leftBlock);
    boost::shared_ptr<SparseMatrix> overlap = rightBlock->get_op_array(OVERLAP).get_local_element(0)[0]->getworkingrepresentation(rightBlock);
    TensorMultiply(leftBlock, *op, *overlap, this, c, *v, op->get_deltaQuantum(0) ,1.0);  // dmrginp.effective_molecule_quantum() is never used in TensorMultiply
    
    overlap = leftBlock->get_op_array(OVERLAP).get_local_element(0)[0]->getworkingrepresentation(leftBlock);
    op = rightBlock->get_op_array(CCD_SUM).get_local_element(0)[0]->getworkingrepresentation(rightBlock);
    TensorMultiply(rightBlock, *op, *overlap, this, c, *v, op->get_deltaQuantum(0), 1.0);  
    //}

  dmrginp.s0time -> stop();
#ifndef SERIAL
  boost::mpi::communicator world;
  int size = world.size();
#endif

  SpinQuantum q= getSpinQuantum( nonactive_orb(0));  dmrginp.s1time -> start();
  v_add =  leftBlock->get_op_array(CCD_CRE_DESCOMP).is_local() ? v_array : v_distributed;
  Functor f = boost::bind(&opxop::ccd_cxcdcomp, leftBlock, _1, this, boost::ref(c), v_add,q); 
  for_all_multithread(rightBlock->get_op_array(CRE), f);

  v_add =  rightBlock->get_op_array(CCD_CRE_DESCOMP).is_local() ? v_array : v_distributed;
  f = boost::bind(&opxop::ccd_cxcdcomp, rightBlock, _1, this, boost::ref(c), v_add, q ); 
  for_all_multithread(leftBlock->get_op_array(CRE), f);  
    
  v_add =  leftBlock->get_op_array(CCD_CRE_CRECOMP).is_local() ? v_array : v_distributed;
  f = boost::bind(&opxop::ccd_dxcccomp, leftBlock, _1, this, boost::ref(c), v_add, q);
  for_all_multithread(rightBlock->get_op_array(DES), f);
  
  v_add =  rightBlock->get_op_array(CCD_CRE_CRECOMP).is_local() ? v_array : v_distributed;
  f = boost::bind(&opxop::ccd_dxcccomp, rightBlock, _1, this, boost::ref(c), v_add, q);
  for_all_multithread(leftBlock->get_op_array(DES), f);
  

  accumulateMultiThread(v, v_array, v_distributed, MAX_THRD);

  dmrginp.oneelecT -> stop();

}


void SpinBlock::multiplyH(Wavefunction& c, Wavefunction* v, int num_threads) const
{

  SpinBlock* loopBlock=(leftBlock->is_loopblock()) ? leftBlock : rightBlock;
  SpinBlock* otherBlock = loopBlock == leftBlock ? rightBlock : leftBlock;

  Wavefunction *v_array=0, *v_distributed=0, *v_add=0;

  int maxt = 1;
  initiateMultiThread(v, v_array, v_distributed, MAX_THRD);
  dmrginp.oneelecT -> start();
  dmrginp.s0time -> start();
  //coreEnergy
  if (fabs(coreEnergy[integralIndex]) > TINY && mpigetrank() == 0) {
    Wavefunction vtemp = *v;
    vtemp.Clear();
    multiplyOverlap(c, &vtemp, num_threads);
    ScaleAdd(coreEnergy[integralIndex], vtemp, *v);
  }
  //if (mpigetrank() == 0) {
    boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_array(HAM).get_local_element(0)[0]->getworkingrepresentation(leftBlock);
    boost::shared_ptr<SparseMatrix> overlap = rightBlock->get_op_array(OVERLAP).get_local_element(0)[0]->getworkingrepresentation(rightBlock);
    TensorMultiply(leftBlock, *op, *overlap, this, c, *v, op->get_deltaQuantum(0) ,1.0);  // dmrginp.effective_molecule_quantum() is never used in TensorMultiply
    overlap = leftBlock->get_op_array(OVERLAP).get_local_element(0)[0]->getworkingrepresentation(leftBlock);
    op = rightBlock->get_op_array(HAM).get_local_element(0)[0]->getworkingrepresentation(rightBlock);
    TensorMultiply(rightBlock, *op, *overlap, this, c, *v, op->get_deltaQuantum(0), 1.0);  
   // }

  dmrginp.s0time -> stop();
#ifndef SERIAL
  boost::mpi::communicator world;
  int size = world.size();
#endif

  dmrginp.s1time -> start();
  v_add =  leftBlock->get_op_array(CRE_CRE_DESCOMP).is_local() && rightBlock->get_op_array(CRE).is_local() ? v_array : v_distributed;
  Functor f = boost::bind(&opxop::cxcddcomp, leftBlock, _1, this, boost::ref(c), v_add, dmrginp.effective_molecule_quantum() ); 
  for_all_multithread(rightBlock->get_op_array(CRE), f);

  v_add =  rightBlock->get_op_array(CRE_CRE_DESCOMP).is_local() && leftBlock->get_op_array(CRE).is_local() ? v_array : v_distributed;
  f = boost::bind(&opxop::cxcddcomp, rightBlock, _1, this, boost::ref(c), v_add, dmrginp.effective_molecule_quantum() ); 
  for_all_multithread(leftBlock->get_op_array(CRE), f);  

  dmrginp.s1time -> stop();

  dmrginp.oneelecT -> stop();

  dmrginp.twoelecT -> start();

  if (dmrginp.hamiltonian() != HUBBARD) {
    
    dmrginp.s0time -> start();
    v_add =  otherBlock->get_op_array(CRE_DESCOMP).is_local() ? v_array : v_distributed;
    f = boost::bind(&opxop::cdxcdcomp, otherBlock, _1, this, boost::ref(c), v_add, dmrginp.effective_molecule_quantum() );
    for_all_multithread(loopBlock->get_op_array(CRE_DES), f);
    
    v_add =  otherBlock->get_op_array(DES_DESCOMP).is_local() ? v_array : v_distributed;
    f = boost::bind(&opxop::ddxcccomp, otherBlock, _1, this, boost::ref(c), v_add, dmrginp.effective_molecule_quantum() );
    for_all_multithread(loopBlock->get_op_array(CRE_CRE), f);
    dmrginp.s0time -> stop();
  }
  
  dmrginp.twoelecT -> stop();


  accumulateMultiThread(v, v_array, v_distributed, MAX_THRD);

}

void SpinBlock::multiplyH_Q(Wavefunction& c, Wavefunction* v, int num_threads, SpinQuantum &Q) const
{

  SpinBlock* loopBlock=(leftBlock->is_loopblock()) ? leftBlock : rightBlock;
  SpinBlock* otherBlock = loopBlock == leftBlock ? rightBlock : leftBlock;

  Wavefunction *v_array=0, *v_distributed=0, *v_add=0;

  int maxt = 1;
  initiateMultiThread(v, v_array, v_distributed, MAX_THRD);
  dmrginp.oneelecT -> start();
  dmrginp.s0time -> start();

  //if (mpigetrank() == 0) {
    boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_array(HAM).get_local_element(0)[0]->getworkingrepresentation(leftBlock);
    boost::shared_ptr<SparseMatrix> overlap = rightBlock->get_op_array(OVERLAP).get_local_element(0)[0]->getworkingrepresentation(rightBlock);
    TensorMultiply(leftBlock, *op, *overlap, this, c, *v, op->get_deltaQuantum(0) ,1.0);  // dmrginp.effective_molecule_quantum() is never used in TensorMultiply
    
    overlap = leftBlock->get_op_array(OVERLAP).get_local_element(0)[0]->getworkingrepresentation(leftBlock);
    op = rightBlock->get_op_array(HAM).get_local_element(0)[0]->getworkingrepresentation(rightBlock);
    TensorMultiply(rightBlock, *op, *overlap, this, c, *v, op->get_deltaQuantum(0), 1.0);  
    //}

  dmrginp.s0time -> stop();
#ifndef SERIAL
  boost::mpi::communicator world;
  int size = world.size();
#endif

  dmrginp.s1time -> start();
  v_add =  leftBlock->get_op_array(CRE_CRE_DESCOMP).is_local() ? v_array : v_distributed;
  Functor f = boost::bind(&opxop::cxcddcomp, leftBlock, _1, this, boost::ref(c), v_add, Q ); 
  for_all_multithread(rightBlock->get_op_array(CRE), f);

  v_add =  rightBlock->get_op_array(CRE_CRE_DESCOMP).is_local() ? v_array : v_distributed;
  f = boost::bind(&opxop::cxcddcomp, rightBlock, _1, this, boost::ref(c), v_add, Q ); 
  for_all_multithread(leftBlock->get_op_array(CRE), f);  

  dmrginp.s1time -> stop();

  dmrginp.oneelecT -> stop();

  dmrginp.twoelecT -> start();

  if (dmrginp.hamiltonian() != HUBBARD) {
    
    dmrginp.s0time -> start();
    v_add =  otherBlock->get_op_array(CRE_DESCOMP).is_local() ? v_array : v_distributed;
    f = boost::bind(&opxop::cdxcdcomp, otherBlock, _1, this, boost::ref(c), v_add, Q );
    for_all_multithread(loopBlock->get_op_array(CRE_DES), f);
    
    v_add =  otherBlock->get_op_array(DES_DESCOMP).is_local() ? v_array : v_distributed;
    f = boost::bind(&opxop::ddxcccomp, otherBlock, _1, this, boost::ref(c), v_add, Q );
    for_all_multithread(loopBlock->get_op_array(CRE_CRE), f);
    dmrginp.s0time -> stop();
  }
  
  dmrginp.twoelecT -> stop();

  accumulateMultiThread(v, v_array, v_distributed, MAX_THRD);

}


void SpinBlock::diagonalH(DiagonalMatrix& e) const
{
  SpinBlock* loopBlock=(leftBlock->is_loopblock()) ? leftBlock : rightBlock;
  SpinBlock* otherBlock = loopBlock == leftBlock ? rightBlock : leftBlock;
  
  DiagonalMatrix *e_array=0, *e_distributed=0, *e_add=0;

  initiateMultiThread(&e, e_array, e_distributed, MAX_THRD);

  boost::shared_ptr<SparseMatrix> op =leftBlock->get_op_array(HAM).get_local_element(0)[0]->getworkingrepresentation(leftBlock);
  TensorTrace(leftBlock, *op, this, &(get_stateInfo()), e, 1.0);
  
  op = rightBlock->get_op_array(HAM).get_local_element(0)[0]->getworkingrepresentation(rightBlock);
  TensorTrace(rightBlock, *op, this, &(get_stateInfo()), e, 1.0);  
  for (int i=0; i<e.Nrows(); i++)
    e(i+1) += coreEnergy[integralIndex];


#ifndef SERIAL
  boost::mpi::communicator world;
  int size = world.size();
#endif

  e_add =  leftBlock->get_op_array(CRE_CRE_DESCOMP).is_local() ? e_array : e_distributed;
  Functor f = boost::bind(&opxop::cxcddcomp_d, leftBlock, _1, this, e_add); 
  for_all_multithread(rightBlock->get_op_array(CRE), f); //not needed in diagonal
  
  e_add =  rightBlock->get_op_array(CRE_CRE_DESCOMP).is_local() ? e_array : e_distributed;
  f = boost::bind(&opxop::cxcddcomp_d, rightBlock, _1, this, e_add); 
  for_all_multithread(leftBlock->get_op_array(CRE), f);  //not needed in diagonal

  if (dmrginp.hamiltonian() != HUBBARD) {
    
    e_add =  otherBlock->get_op_array(CRE_DESCOMP).is_local() ? e_array : e_distributed;
    f = boost::bind(&opxop::cdxcdcomp_d, otherBlock, _1, this, e_add);
    for_all_multithread(loopBlock->get_op_array(CRE_DES), f);  
    
    e_add =  otherBlock->get_op_array(DES_DESCOMP).is_local() ? e_array : e_distributed;
    f = boost::bind(&opxop::ddxcccomp_d, otherBlock, _1, this, e_add);
    for_all_multithread(loopBlock->get_op_array(CRE_CRE), f);  //not needed in diagonal
  }

  accumulateMultiThread(&e, e_array, e_distributed, MAX_THRD);

}

void SpinBlock::BuildSlaterBlock (std::vector<int> sts, std::vector<SpinQuantum> qnumbers, std::vector<int> distribution, bool random, const bool haveNormops)
{
  name = get_name();

  if (dmrginp.spinAdapted()) {
    sites = sts;
  }
  else {
    for (int i=0; i<sts.size(); i++) {
      sites.push_back( dmrginp.spatial_to_spin()[sts[i]]   );
      sites.push_back( dmrginp.spatial_to_spin()[sts[i]]+1 );
    }
  }

  complementary_sites = make_complement(sites);

  assert (sites.size () > 0);
  sort (sites.begin (), sites.end ());

  //always have implicit transpose in this case
  default_op_components(!haveNormops, true);
  
  setstoragetype(DISTRIBUTED_STORAGE);


  std::vector< Csf > dets;
  std::vector< Csf > det_ex;
  Timer slatertimer;

  for (int i = 0; i < qnumbers.size (); ++i)
    {
      if (distribution [i] == 0 || (qnumbers [i].get_n() > 2*sites.size ())) continue;

      if(dmrginp.spinAdapted()) 
	det_ex = Csf::distribute (qnumbers [i].get_n(), qnumbers [i].get_s().getirrep(), IrrepVector(qnumbers [i].get_symm().getirrep(), 0) , sites [0],
				  sites [0] + sites.size (), dmrginp.last_site(), integralIndex);
      else
	det_ex = Csf::distributeNonSpinAdapted (qnumbers [i].get_n(), qnumbers [i].get_s().getirrep(), IrrepVector(qnumbers [i].get_symm().getirrep(), 0) , sites [0],
						sites [0] + sites.size (), dmrginp.last_site(), integralIndex);

      multimap <double, Csf > slater_emap;

      for (int j = 0; j < det_ex.size(); ++j) {
	slater_emap.insert (pair <double, Csf > (csf_energy (det_ex[j], integralIndex), det_ex[j]));
      }

      multimap <double, Csf >::iterator m = slater_emap.begin();
      int sz = det_ex.size();
      det_ex.resize (min (distribution [i], sz));
      for (int j = 0; j < det_ex.size(); ++j)
        {
          det_ex[j] = m->second;
          ++m;
        }

      copy (det_ex.begin(), det_ex.end(), back_inserter (dets));
    }

  std::multimap<SpinQuantum, Csf > tmp;
  for (int i = 0; i < dets.size (); ++i) {
    tmp.insert(pair<SpinQuantum, Csf > (SpinQuantum ( (dets [i]).n, (dets [i]).S, (dets[i]).sym_is()), dets[i]));
  }
  std::multimap<SpinQuantum, Csf >::iterator tmpiter = tmp.begin();
  dets.clear();
  for (; tmpiter != tmp.end();)
  {
    dets.push_back(tmpiter->second);
    tmpiter++;
  }

  braStateInfo = StateInfo (dets);
  ketStateInfo = StateInfo (dets);

  twoInt = boost::shared_ptr<TwoElectronArray>( &v_2[integralIndex], boostutils::null_deleter());
  build_iterators();

  tcpu = slatertimer.elapsedcputime(); twall= slatertimer.elapsedwalltime();
  p3out << "\t\t\t time in slater distribution " << twall << " " << tcpu << endl;

  std::vector< std::vector<Csf> > ladders; ladders.resize(dets.size());
  for (int i=0; i< dets.size(); i++)
    ladders[i] = dets[i].spinLadder(min(2, dets[i].S.getirrep()));


  build_operators(dets, ladders);
  tcpu = slatertimer.elapsedcputime();twall= slatertimer.elapsedwalltime();
  p3out << "\t\t\t time in slater operator build " <<  twall << " " << tcpu << endl;


}

void SpinBlock::BuildSingleSlaterBlock(std::vector<int> sts) {
  name = get_name();
  if (dmrginp.spinAdapted()) {
    sites = sts;
  }
  else {
    for (int i=0; i<sts.size(); i++) {
      sites.push_back( dmrginp.spatial_to_spin()[sts[i]]   );
      sites.push_back( dmrginp.spatial_to_spin()[sts[i]]+1 );
    }
  }

  complementary_sites = make_complement(sites);
  assert (sites.size () > 0);
  sort (sites.begin (), sites.end ());

  int left = sites[0], right = sites[0] + sites.size(), edge = dmrginp.last_site();
  int n = 0, sp = 0;
  std::vector<bool> tmp(0);
  IrrepSpace irrep(0);

  if (dmrginp.spinAdapted()) {
    for (int orbI = left; orbI < right; ++orbI) {
      n += dmrginp.hf_occupancy()[dmrginp.spatial_to_spin()[orbI]] + dmrginp.hf_occupancy()[dmrginp.spatial_to_spin()[orbI]+1];
      sp += dmrginp.hf_occupancy()[dmrginp.spatial_to_spin()[orbI]] - dmrginp.hf_occupancy()[dmrginp.spatial_to_spin()[orbI]+1];

      // FIXME: NN wrote, follows don't work correctly for non-abelian symmetry
      if (dmrginp.hf_occupancy()[dmrginp.spatial_to_spin()[orbI]] == 1) {
        irrep = IrrepSpace(Symmetry::add(irrep.getirrep(),SymmetryOfSpatialOrb(orbI).getirrep())[0]);
      }
      if (dmrginp.hf_occupancy()[dmrginp.spatial_to_spin()[orbI]+1] == 1) {
        irrep = IrrepSpace(Symmetry::add(irrep.getirrep(),SymmetryOfSpatialOrb(orbI).getirrep())[0]);
      }
    }

    for (int i = 0; i < dmrginp.spatial_to_spin()[left]; ++i) {
      tmp.push_back(0);
    }
    for (int orbI = left; orbI < right; ++orbI) {
      tmp.push_back(dmrginp.hf_occupancy()[dmrginp.spatial_to_spin()[orbI]]);
      tmp.push_back(dmrginp.hf_occupancy()[dmrginp.spatial_to_spin()[orbI]+1]);
    }
    for (int i = 0; i < dmrginp.spatial_to_spin()[edge]-dmrginp.spatial_to_spin()[right]; ++i) {
      tmp.push_back(0);
    }
  } else {
    for (int i = 0; i < left; ++i) {
      tmp.push_back(0);
    }
    for (int orbI = left; orbI < right; ++orbI) {
      if (dmrginp.hf_occupancy()[orbI] == 1) {
        n += 1;
        sp += SpinOf(orbI);
      }
      tmp.push_back(dmrginp.hf_occupancy()[orbI]);
    }
    for (int i = 0; i < edge-right; ++i) {
      tmp.push_back(0);
    }
  }

  Slater new_det = Slater(Orbstring(tmp));
  map<Slater, double> m;
  m[new_det] = 1.0;

  Csf origin(m, n, SpinSpace(sp), sp, IrrepVector(irrep.getirrep(), 0));
  std::vector<Csf> dets(1, origin);
  braStateInfo = StateInfo (dets);
  ketStateInfo = StateInfo (dets);
  twoInt = boost::shared_ptr<TwoElectronArray>( &v_2[integralIndex], boostutils::null_deleter());
  build_iterators();

  std::vector< std::vector<Csf> > ladders; ladders.resize(dets.size());
  for (int i=0; i< dets.size(); i++)
    ladders[i] = dets[i].spinLadder(min(2, dets[i].S.getirrep()));

  build_operators(dets, ladders);

}
}
