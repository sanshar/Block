/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "spinblock.h"
#include "wavefunction.h"
#include <boost/format.hpp>
#include <fstream>
#include <stdio.h>

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
#include "pario.h"
namespace SpinAdapted{


std::string SpinBlock::restore (bool forward, const vector<int>& sites, SpinBlock& b)
{
  Timer disktimer;
  std::string file;
  if (forward)
    file = str(boost::format("%s%s%d%s%d%s%d%s") % dmrginp.load_prefix() % "/SpinBlock-forward-" % sites[0] % "-" % sites[sites.size()-1] % "." % mpigetrank() % ".tmp" );
  else
    file = str(boost::format("%s%s%d%s%d%s%d%s") % dmrginp.load_prefix() % "/SpinBlock-backward-"% sites[0] % "-" % sites[sites.size()-1] % "." % mpigetrank() % ".tmp" );

  if (dmrginp.outputlevel() > 0) 
    pout << "\t\t\t Restoring block file :: " << file << endl;

  std::ifstream ifs(file.c_str(), std::ios::binary);
  //coutbuf = &ifs;
  b.Load (ifs);
  ifs.close();
  //coutbuf = 0;
  return file;
}


void SpinBlock::store (bool forward, const vector<int>& sites, SpinBlock& b)
{
  Timer disktimer;
  std::string file;
  if (forward)
    file = str(boost::format("%s%s%d%s%d%s%d%s") % dmrginp.save_prefix() % "/SpinBlock-forward-"% sites[0] % "-" % sites[sites.size()-1] % "." % mpigetrank() % ".tmp" );
  else
    file = str(boost::format("%s%s%d%s%d%s%d%s") % dmrginp.save_prefix() % "/SpinBlock-backward-"% sites[0] % "-" % sites[sites.size()-1] % "." % mpigetrank() % ".tmp" );
  if (dmrginp.outputlevel() > 0) 
    pout << "\t\t\t Saving block file :: " << file << endl;


  std::ofstream ofs(file.c_str(), std::ios::binary);
  b.Save (ofs);
  ofs.close();
  //pout << "\t\t\t block save disk time " << disktimer.elapsedwalltime() << " " << disktimer.elapsedcputime() << endl;
}

void SpinBlock::Save (std::ofstream &ofs)
{
  boost::archive::binary_oarchive save_block(ofs);
  save_block << *this;
}

//helper function
void SpinBlock::Load (std::ifstream & ifs)
{
  boost::archive::binary_iarchive load_block(ifs);
  load_block >> *this;
}


void SpinBlock::remove_normal_ops()
{
  ops.erase(CRE_CRE);
  ops.erase(CRE_DES);
}

void SpinBlock::sendcompOps(Op_component_base& opcomp, int I, int J, opTypes optype, int compsite)
{
#ifndef SERIAL
  boost::mpi::communicator world;
  std::vector<boost::shared_ptr<SparseMatrix> > oparray = opcomp.get_element(I,J);
  for(int i=0; i<oparray.size(); i++) {
    world.send(processorindex(compsite), (optype & OP_TYPE_MASK)+i+100*J+10000*I, *oparray[i]);
  }  
#endif
}

void SpinBlock::recvcompOps(Op_component_base& opcomp, int I, int J, opTypes optype)
{
#ifndef SERIAL
  boost::mpi::communicator world;
  std::vector<boost::shared_ptr<SparseMatrix> > oparray = opcomp.get_element(I,J);
  for(int i=0; i<oparray.size(); i++) {
    world.recv(processorindex(trimap(I, J, dmrginp.last_site())), (optype & OP_TYPE_MASK)+i+100*J+10000*I, *oparray[i]);
  }
#endif
}

void SpinBlock::addAdditionalCompOps(int iState, int jState)
{
#ifndef SERIAL
  boost::mpi::communicator world;
  if (world.size() == 1)
    return; //there is no need to have additional compops

  int length = dmrginp.last_site();
  int dotopindex = (sites[0] == 0) ? complementary_sites[0] : complementary_sites[complementary_sites.size()-1];

  // make op-type index for < iState | O | jState >
  opTypes state_index = make_state_index(iState, jState);

  if (!ops[CRE+state_index]->is_local()) {
    for(int i=0; i<get_sites().size(); i++) {
      if (ops[CRE+state_index]->has(sites[i])) {
	if (processorindex(sites[i]) != mpigetrank()) 
	  ops[CRE+state_index]->add_local_indices(sites[i]);
	mpi::broadcast(world, *(ops[CRE+state_index]->get_element(sites[i])[0]), processorindex(sites[i]));
      }
    }
  }

  for (int i=0; i<complementary_sites.size(); i++) {
    int compsite = complementary_sites[i];
    if (compsite == dotopindex) continue;
    int I = (compsite > dotopindex) ? compsite : dotopindex;
    int J = (compsite > dotopindex) ? dotopindex : compsite;
    if (processorindex(compsite) == processorindex(trimap(I, J, length)))
      continue;
    if (processorindex(compsite) == mpigetrank())
    {
      bool other_proc_has_ops = true;
      world.recv(processorindex(trimap(I, J, length)), 0, other_proc_has_ops);
      //this will potentially receive some ops
      if (other_proc_has_ops) {
	ops[CRE_DESCOMP+state_index]->add_local_indices(I, J);
	recvcompOps(*ops[CRE_DESCOMP+state_index], I, J, CRE_DESCOMP);
	ops[DES_DESCOMP+state_index]->add_local_indices(I, J);
	recvcompOps(*ops[DES_DESCOMP+state_index], I, J, DES_DESCOMP);
      }
    }
    else
    {
      
      //this will potentially send some ops
      if (processorindex(trimap(I, J, length)) == mpigetrank()) {
	bool this_proc_has_ops = ops[CRE_DESCOMP+state_index]->has_local_index(I, J);
	world.send(processorindex(compsite), 0, this_proc_has_ops);
	if (this_proc_has_ops) {
	  sendcompOps(*ops[CRE_DESCOMP+state_index], I, J, CRE_DESCOMP, compsite);
	  sendcompOps(*ops[DES_DESCOMP+state_index], I, J, DES_DESCOMP, compsite);
	}
      }
      else 
	continue;
    }
    //dmrginp.datatransfer.stop();
    //dmrginp.datatransfer -> stop(); //ROA
      
  }
#endif
}

void SpinBlock::transform_operators(std::vector<Matrix>& rotateMatrix)
{
  StateInfo oldStateInfo = stateInfo;
  std::vector<SpinQuantum> newQuanta;
  std::vector<int> newQuantaStates;
  std::vector<int> newQuantaMap;
  for (int Q = 0; Q < rotateMatrix.size (); ++Q)
    {
      if (rotateMatrix [Q].Ncols () != 0)
        {
          newQuanta.push_back (stateInfo.quanta [Q]);
          newQuantaStates.push_back (rotateMatrix [Q].Ncols ());
          newQuantaMap.push_back (Q);
        }
    }
  StateInfo newStateInfo = StateInfo (newQuanta, newQuantaStates, newQuantaMap);

  for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
    if (! it->second->is_core())
      for_all_operators_multithread(*it->second, bind(&SparseMatrix::build_and_renormalise_transform, _1, this, it->first, 
						       ref(rotateMatrix) , &newStateInfo));
  stateInfo = newStateInfo;
  stateInfo.AllocatePreviousStateInfo ();
  *stateInfo.previousStateInfo = oldStateInfo;

  for (int i = 0; i < newQuantaMap.size (); ++i)
    assert (stateInfo.quanta [i] == oldStateInfo.quanta [newQuantaMap [i]]);

  if (dmrginp.outputlevel() > 0) {
    pout << "\t\t\t total elapsed time " << globaltimer.totalwalltime() << " " << globaltimer.totalcputime() << " ... " 
	 << globaltimer.elapsedwalltime() << " " << globaltimer.elapsedcputime() << endl;
    pout << "\t\t\t Transforming to new basis " << endl;
  }
  Timer transformtimer;

  for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
    if ( it->second->is_core())
      for_all_operators_multithread(*it->second, bind(&SparseMatrix::renormalise_transform, _1, ref(rotateMatrix), (&this->stateInfo)));


  for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
    if (! it->second->is_core())
      ops[it->first]->set_core(true);

  this->direct = false;
  if (dmrginp.outputlevel() > 0)
    pout << "\t\t\t transform time " << transformtimer.elapsedwalltime() << " " << transformtimer.elapsedcputime() << endl;

  if (leftBlock)
    leftBlock->clear();
  if (rightBlock)
    rightBlock->clear();


}
}
