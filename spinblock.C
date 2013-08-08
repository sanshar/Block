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
         cout << it->second->size()<<" :  "<<it->second->get_op_string()<<"  Core Operators  ";      
      else
         cout << it->second->size()<<" :  "<<it->second->get_op_string()<<"  Virtual Operators  ";      
      
      vector<int> numops(world.size(), 0);
      for (int proc = 0; proc <world.size(); proc++) {
         if (proc != 0) 
            receiveobject(numops[proc],proc);
         else 
            numops[proc] = it->second->get_size();
         cout <<numops[proc]<<"  ";
      }
      cout << endl;
    }
  }
#else
  for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::const_iterator it = ops.begin(); it != ops.end(); ++it)
  {
    if(it->second->is_core()) 
      cout << it->second->size()<<" :  "<<it->second->get_op_string()<<"  Core Operators  ";      
    else
      cout << it->second->size()<<" :  "<<it->second->get_op_string()<<"  Virtual Operators  ";      
    cout << endl;
  }
#endif
  
}
ostream& operator<< (ostream& os, const SpinBlock& b)
{
  os << "\t\t\t Sites ::  ";
  for (int i = 0; i < b.sites.size(); ++i) { os << b.sites[i] << " "; } 
  
  if (dmrginp.outputlevel() > 0) {
    os << endl;
    os << b.stateInfo;
  }
  else {
    os <<"    # states: "<<b.stateInfo.totalStates<<endl;
  }
  return os;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Constructors

SpinBlock::SpinBlock () : 
  localstorage(false),
  name (rand()), 
  hasMemoryAllocated (false),
  direct(false), complementary(false), normal(true), leftBlock(0), rightBlock(0) { }

SpinBlock::SpinBlock(int start, int finish, bool is_complement) :  
  name (rand()), 
  hasMemoryAllocated (false), 
  direct(false), leftBlock(0), rightBlock(0)
{
  complementary = is_complement;
  normal = !is_complement;
//MAW
pout << "calling default op_components...\n";
  default_op_components(is_complement);
  std::vector<int> sites; 
  if (dmrginp.use_partial_two_integrals()) {
    if (start != finish) {
      cout << "Cannot use partial two electron integrals, when making spin block with more than two orbitals"<<endl;
      abort();
    }
    std::vector<int> o;
    for (int i=dmrginp.spatial_to_spin()[start]; i<dmrginp.spatial_to_spin()[start+1]; i+=2)
      o.push_back(i/2);
    twoInt = boost::shared_ptr<PartialTwoElectronArray> (new PartialTwoElectronArray(o));
    twoInt->Load(dmrginp.load_prefix());
#ifndef SERIAL
    mpi::communicator world;
    PartialTwoElectronArray& ar = dynamic_cast<PartialTwoElectronArray&>(*twoInt.get());
    mpi::broadcast(world, ar, 0);
#endif
  }
  else
    twoInt = boost::shared_ptr<TwoElectronArray>( &v_2, boostutils::null_deleter());

  int lower = min(start, finish);
  int higher = max(start, finish);
  sites.resize(higher - lower + 1);
  for (int i=0; i < sites.size(); i++)
      sites[i] = lower + i;

//MAW
pout << "calling BuildTensorProductBlock(sites)...\n";   
  BuildTensorProductBlock(sites);   
pout << "done BuildTensorProductBlock(sites)!\n";   
}

SpinBlock::SpinBlock (const SpinBlock& b) { *this = b; }

SpinBlock::SpinBlock(const StateInfo& s)
{
  stateInfo = s;
  sites.resize(0);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinBlock::BuildTensorProductBlock(std::vector<int>& new_sites)
{
pout << "SpinBlock::BuildTensorProductBlock(std::vector<int>& new_sites)\n";

  if (twoInt.get() == 0 && dmrginp.use_partial_two_integrals()) { //this is when dummy block is being added for non zero spin
    std::vector<int> o;
    for (int i=dmrginp.spatial_to_spin()[new_sites[0]]; i<dmrginp.spatial_to_spin()[new_sites[new_sites.size()-1]+1]; i+=2)
      o.push_back(i/2);
    twoInt = boost::shared_ptr<PartialTwoElectronArray> (new PartialTwoElectronArray(o));
    twoInt->Load(dmrginp.load_prefix());
#ifndef SERIAL
    mpi::communicator world;
    PartialTwoElectronArray& ar = dynamic_cast<PartialTwoElectronArray&>(*twoInt.get());
    mpi::broadcast(world, ar, 0);
    //world.broadcast(twoInt);
#endif
  }
  else if (twoInt.get() == 0)
    twoInt = boost::shared_ptr<TwoElectronArray>( &v_2, boostutils::null_deleter());

  name = get_name();

  sites = new_sites;
  std::vector< Csf > dets = CSFUTIL::spinfockstrings(new_sites);

  stateInfo = StateInfo(dets);
//MAW note LOCAL_STORAGE for this (used by Dot block?)
  setstoragetype(LOCAL_STORAGE);
  complementary_sites = make_complement(sites);
//MAW
pout << "build_iterators....\n";
  build_iterators();
  std::vector< std::vector<Csf> > ladders; ladders.resize(dets.size());
  for (int i=0; i< dets.size(); i++)
    ladders[i] = dets[i].spinLadder(min(2,dets[i].S));
//MAW
pout << "build_operators(dets,ladders)....\n";
  build_operators(dets, ladders);
pout << "done SpinBlock::BuildTensorProductBlock(std::vector<int>& new_sites)\n";
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

std::vector<int> SpinBlock::make_complement(const std::vector<int>& sites)
{
  std::vector<int> complementary_sites;
  for (int i=0; i<dmrginp.last_site(); ++i)
    if (find(sites.begin(), sites.end(), i) == sites.end())
      complementary_sites.push_back(i);

  return complementary_sites;
  
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinBlock::build_iterators()
{
  for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
  {
    it->second->build_iterators(*this);
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinBlock::build_operators(std::vector< Csf >& dets, std::vector< std::vector<Csf> >& ladders)
{
  for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
    {
      if(it->second->is_core()) {
        // Output file for operators written to disk
        std::string ofile = it->second->get_filename();
        // Build operators from CSFs
        it->second->build_csf_operators(*this, ofile, dets, ladders);
      }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinBlock::build_virtual_operators()
{
  for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
    {
      opTypes ot = it->first;
//MAW note ! here
      if(! it->second->is_core()) {
        // Output file for operators written to disk
        std::string ofile = it->second->get_filename();
        // Input file for operators written to disk on sysblock
        std::string sysfile = get_leftBlock()->ops[ot]->get_filename();
        // Input file for operators written to disk on dotblock
        std::string dotfile = get_rightBlock()->ops[ot]->get_filename();
        // Build operators
        it->second->build_operators(*this, ot, ofile, sysfile, dotfile);
      }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinBlock::build_operators()
{
  for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
    {
      opTypes ot = it->first;
      if(it->second->is_core()) {
        // Output file for operators written to disk
        std::string ofile = it->second->get_filename();
        // Input file for operators written to disk on sysblock
        std::string sysfile = get_leftBlock()->ops[ot]->get_filename();
        // Input file for operators written to disk on dotblock
        std::string dotfile = get_rightBlock()->ops[ot]->get_filename();
        // Build operators
        it->second->build_operators(*this, ot, ofile, sysfile, dotfile);
      }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinBlock::renormalise_transform(const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo) 
{
  for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
    {
      if(it->second->is_core()) {
        it->second->renormalise_transform(rotateMatrix, stateinfo);
      }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinBlock::clear()
{
  for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
    it->second->clear();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinBlock::BuildSumBlockSkeleton(int condition, SpinBlock& lBlock, SpinBlock& rBlock, StateInfo* compState)
{
  name = get_name();
  if (dmrginp.outputlevel() > 0) 
    pout << "\t\t\t Building Sum Block " << name << endl;
  leftBlock = &lBlock;
  rightBlock = &rBlock;

  sites.reserve (lBlock.sites.size () + rBlock.sites.size ());

  if (dmrginp.use_partial_two_integrals()) {
    if (rBlock.sites.size() == 1) {
      std::vector<int> o;
      for (int i=dmrginp.spatial_to_spin().at(rBlock.sites[0]); i<dmrginp.spatial_to_spin().at(rBlock.sites[0]+1); i+=2)
	o.push_back(i/2);
      twoInt = boost::shared_ptr<PartialTwoElectronArray> (new PartialTwoElectronArray(o));
      twoInt->Load(dmrginp.load_prefix());
#ifndef SERIAL
      mpi::communicator world;
      PartialTwoElectronArray& ar = dynamic_cast<PartialTwoElectronArray&>(*twoInt.get());
      mpi::broadcast(world, ar, 0);
#endif

    }
    //cout << "Cannot use partial two electron integrals, when the dot block has more than one orbital"<<endl;
    //abort();

  }
  else 
    twoInt = boost::shared_ptr<TwoElectronArray>( &v_2,  boostutils::null_deleter());


  sites = lBlock.sites;
  copy (rBlock.sites.begin(), rBlock.sites.end (), back_inserter (sites));
  sort(sites.begin(), sites.end());
  complementary_sites = make_complement(sites);
  if (dmrginp.outputlevel() > 0) {
    pout << "\t\t\t ";
    for (int i = 0; i < sites.size(); ++i) pout << sites[i] << " ";
    pout << endl;
  }

  TensorProduct (lBlock.stateInfo, rBlock.stateInfo, stateInfo, condition, compState);
  if (condition != PARTICLE_SPIN_NUMBER_CONSTRAINT)
    stateInfo.CollectQuanta();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinBlock::BuildSumBlock(int condition, SpinBlock& lBlock, SpinBlock& rBlock, StateInfo* compState)
{
  dmrginp.buildsumblock -> start();
  BuildSumBlockSkeleton(condition, lBlock, rBlock, compState);

//pout << "maw SpinBlock::BuildSumBlock build_iterators()\n";
  build_iterators();

  dmrginp.buildblockops -> start();
//pout << "maw SpinBlock::BuildSumBlock build_operators()\n";
  // Operators builts from previous operators (not CSF!) (c.f.  build_operators(dets, ladders); )
  build_operators();
  dmrginp.buildblockops -> stop();
  dmrginp.buildsumblock -> stop();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

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

  direct = b.is_direct();

  stateInfo = b.stateInfo;
  leftBlock = b.leftBlock;
  rightBlock = b.rightBlock;
  twoInt = b.twoInt;
  ops = b.ops;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinBlock::multiplyH(Wavefunction& c, Wavefunction* v, int num_threads) const
{

  SpinBlock* loopBlock=(leftBlock->is_loopblock()) ? leftBlock : rightBlock;
  SpinBlock* otherBlock = loopBlock == leftBlock ? rightBlock : leftBlock;

  Wavefunction *v_array=0, *v_distributed=0, *v_add=0;
  
  int maxt = 1;
  initiateMultiThread(v, v_array, v_distributed, MAX_THRD);
  dmrginp.oneelecT -> start();
  dmrginp.s0time -> start();
  boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_array(HAM).get_local_element(0)[0];
  TensorMultiply(leftBlock, *op, this, c, *v, dmrginp.effective_molecule_quantum() ,1.0, MAX_THRD);

  op = rightBlock->get_op_array(HAM).get_local_element(0)[0];
  TensorMultiply(rightBlock, *op, this, c, *v, dmrginp.effective_molecule_quantum(), 1.0, MAX_THRD);  

  dmrginp.s0time -> stop();
#ifndef SERIAL
  boost::mpi::communicator world;
  int size = world.size();
#endif

  dmrginp.s1time -> start();
  v_add =  leftBlock->get_op_array(CRE_CRE_DESCOMP).is_local() ? v_array : v_distributed;
  Functor f = boost::bind(&opxop::cxcddcomp, leftBlock, _1, this, boost::ref(c), v_add, dmrginp.effective_molecule_quantum() ); 
  for_all_multithread(rightBlock->get_op_array(CRE), f);

  v_add =  rightBlock->get_op_array(CRE_CRE_DESCOMP).is_local() ? v_array : v_distributed;
  f = boost::bind(&opxop::cxcddcomp, rightBlock, _1, this, boost::ref(c), v_add, dmrginp.effective_molecule_quantum() ); 
  for_all_multithread(leftBlock->get_op_array(CRE), f);  

  dmrginp.s1time -> stop();

  dmrginp.oneelecT -> stop();

  dmrginp.twoelecT -> start();

  if (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY) {
    
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

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinBlock::diagonalH(DiagonalMatrix& e) const
{
  SpinBlock* loopBlock=(leftBlock->is_loopblock()) ? leftBlock : rightBlock;
  SpinBlock* otherBlock = loopBlock == leftBlock ? rightBlock : leftBlock;

  DiagonalMatrix *e_array=0, *e_distributed=0, *e_add=0;

  initiateMultiThread(&e, e_array, e_distributed, MAX_THRD);

  boost::shared_ptr<SparseMatrix> op =leftBlock->get_op_array(HAM).get_local_element(0)[0]->getworkingrepresentation(this);
  TensorTrace(leftBlock, *op, this, &(get_stateInfo()), e, 1.0);

  op = rightBlock->get_op_array(HAM).get_local_element(0)[0]->getworkingrepresentation(this);
  TensorTrace(rightBlock, *op, this, &(get_stateInfo()), e, 1.0);  


#ifndef SERIAL
  boost::mpi::communicator world;
  int size = world.size();
#endif

  std::vector< std::vector<int> > indices;
  e_add =  leftBlock->get_op_array(CRE_CRE_DESCOMP).is_local() ? e_array : e_distributed;
  indices = rightBlock->get_op_array(CRE).get_array();
  Functor f = boost::bind(&opxop::cxcddcomp_d, leftBlock, _1, this, e_add); 
  //for_all_multithread(rightBlock->get_op_array(CRE), f); //not needed in diagonal


  e_add =  rightBlock->get_op_array(CRE_CRE_DESCOMP).is_local() ? e_array : e_distributed;
  indices = leftBlock->get_op_array(CRE).get_array();
  f = boost::bind(&opxop::cxcddcomp_d, rightBlock, _1, this, e_add); 
  //for_all_multithread(leftBlock->get_op_array(CRE), f);  //not needed in diagonal



  if (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY) {
    
    e_add =  otherBlock->get_op_array(CRE_DESCOMP).is_local() ? e_array : e_distributed;
    indices = loopBlock->get_op_array(CRE_DES).get_array();
    f = boost::bind(&opxop::cdxcdcomp_d, otherBlock, _1, this, e_add);
    for_all_multithread(loopBlock->get_op_array(CRE_DES), f);  

    
    e_add =  otherBlock->get_op_array(DES_DESCOMP).is_local() ? e_array : e_distributed;
    indices = loopBlock->get_op_array(CRE_CRE).get_array();
    f = boost::bind(&opxop::ddxcccomp_d, otherBlock, _1, this, e_add);
    //for_all_multithread(loopBlock->get_op_array(CRE_CRE), f);  //not needed in diagonal 

  }

  accumulateMultiThread(&e, e_array, e_distributed, MAX_THRD);

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinBlock::BuildSlaterBlock (std::vector<int> sts, std::vector<SpinQuantum> qnumbers, std::vector<int> distribution, 
                                  bool random, const bool haveNormops)
{
  name = get_name();

  sites = sts;
  complementary_sites = make_complement(sites);

  assert (sites.size () > 0);
  sort (sites.begin (), sites.end ());

  default_op_components(!haveNormops);
  
  setstoragetype(DISTRIBUTED_STORAGE);


  std::vector< Csf > dets;
  std::vector< Csf > det_ex;
  Timer slatertimer;

  for (int i = 0; i < qnumbers.size (); ++i) {
      if (distribution [i] == 0 || (qnumbers [i].get_n() > 2*sites.size ())) continue;
      det_ex = Csf::distribute (qnumbers [i].get_n(), qnumbers [i].get_s(), IrrepVector(qnumbers [i].get_symm().getirrep(), 0) , sites [0],
                                   sites [0] + sites.size (), dmrginp.last_site());

      multimap <double, Csf > slater_emap;

      for (int j = 0; j < det_ex.size(); ++j) {
        slater_emap.insert (pair <double, Csf > (csf_energy (det_ex[j]), det_ex[j]));
      }

      multimap <double, Csf >::iterator m = slater_emap.begin();
      int sz = det_ex.size();
      det_ex.resize (min (distribution [i], sz));
      for (int j = 0; j < det_ex.size(); ++j) {
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

  stateInfo = StateInfo (dets);
  twoInt = boost::shared_ptr<TwoElectronArray>( &v_2, boostutils::null_deleter());
//MAW
  build_iterators();

  if (dmrginp.outputlevel() > 0) 
    pout << "\t\t\t time in slater distribution " << slatertimer.elapsedwalltime() << " " << slatertimer.elapsedcputime() << endl;

  std::vector< std::vector<Csf> > ladders; ladders.resize(dets.size());
  for (int i=0; i< dets.size(); i++)
    ladders[i] = dets[i].spinLadder(min(2, dets[i].S));

//MAW
  build_operators(dets, ladders);
  if (dmrginp.outputlevel() > 0) 
    pout << "\t\t\t time in slater operator build " << slatertimer.elapsedwalltime() << " " << slatertimer.elapsedcputime() << endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

}
