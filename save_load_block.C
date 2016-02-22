/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "sweep.h"
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


std::string SpinBlock::restore (bool forward, const vector<int>& sites, SpinBlock& b, int left, int right, char* name)
{
  Timer disktimer;
  std::string file;

  if (forward)
    file = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%s") % dmrginp.save_prefix() % "/SpinBlock-forward-"% sites[0] % "-" % sites[sites.size()-1] % "." % left % "." % right % "." %b.integralIndex % "." % mpigetrank() % ".tmp" );
  else
    file = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%s") % dmrginp.save_prefix() % "/SpinBlock-backward-"% sites[0] % "-" % sites[sites.size()-1] % "." % left % "." % right % "." %b.integralIndex % "." % mpigetrank() % ".tmp" );
  
  p1out << "\t\t\t Restoring block file :: " << file << endl;

  std::ifstream ifs(file.c_str(), std::ios::binary);

  int lstate =  left;
  int rstate =  right;

  if (mpigetrank() == 0) {
    StateInfo::restore(forward, sites, b.braStateInfo, lstate);
    StateInfo::restore(forward, sites, b.ketStateInfo, rstate);
  }
  
#ifndef SERIAL
  mpi::communicator world;
  mpi::broadcast(world, b.braStateInfo, 0);
  mpi::broadcast(world, b.ketStateInfo, 0);
#endif

  b.Load (ifs);
  ifs.close();


  return file;
}
  
void SpinBlock::store (bool forward, const vector<int>& sites, SpinBlock& b, int left, int right, char *name)
{
  Timer disktimer;
  std::string file;
  if(dmrginp.spinAdapted()) {
    if (forward)
      file = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%s") % dmrginp.save_prefix() % "/SpinBlock-forward-"% sites[0] % "-" % sites[sites.size()-1] % "." % left % "." % right % "." %b.integralIndex % "." % mpigetrank() % ".tmp" );
    else
      file = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%s") % dmrginp.save_prefix() % "/SpinBlock-backward-"% sites[0] % "-" % sites[sites.size()-1] % "." % left % "." % right % "." %b.integralIndex % "." % mpigetrank() % ".tmp" );
  }
  else {
    if (forward)
      file = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%s") % dmrginp.save_prefix() % "/SpinBlock-forward-"% (sites[0]/2) % "-" % (sites[sites.size()-1]/2) % "." % left % "." % right % "." %b.integralIndex % "." % mpigetrank() % ".tmp" );
    else
      file = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%s") % dmrginp.save_prefix() % "/SpinBlock-backward-"% (sites[0]/2) % "-" % (sites[sites.size()-1]/2) % "." % left % "." % right % "." %b.integralIndex % "." % mpigetrank() % ".tmp" );
  }
  
  p1out << "\t\t\t Saving block file :: " << file << endl;
  
  
  std::ofstream ofs(file.c_str(), std::ios::binary);
  
  int lstate =  left;
  int rstate =  right;
  
  if (mpigetrank()==0) {
    StateInfo::store(forward, sites, b.braStateInfo, lstate);
    StateInfo::store(forward, sites, b.ketStateInfo, rstate);
  }

  b.Save (ofs);
  ofs.close(); 
  //p1out << "\t\t\t block save disk time " << disktimer.elapsedwalltime() << " " << disktimer.elapsedcputime() << endl;
}

void SpinBlock::Save (std::ofstream &ofs)
{
  dmrginp.diskio->start();
  boost::archive::binary_oarchive save_block(ofs);
  save_block << *this;
  dmrginp.diskio->stop();
}

//helper function
void SpinBlock::Load (std::ifstream & ifs)
{
  dmrginp.diskio->start();
  boost::archive::binary_iarchive load_block(ifs);
  load_block >> *this;
  dmrginp.diskio->stop();
}


void SpinBlock::remove_normal_ops()
{
  ops.erase(CRE_CRE);
  ops.erase(CRE_DES);
}

void SpinBlock::sendcompOps(Op_component_base& opcomp, int I, int J, int optype, int compsite)
{
#ifndef SERIAL
  boost::mpi::communicator world;
  std::vector<boost::shared_ptr<SparseMatrix> > oparray = opcomp.get_element(I,J);
  for(int i=0; i<oparray.size(); i++) {
    world.send(processorindex(compsite), optype+i*10+1000*J+100000*I, *oparray[i]);
  }  
#endif
}

void SpinBlock::recvcompOps(Op_component_base& opcomp, int I, int J, int optype)
{
#ifndef SERIAL
  boost::mpi::communicator world;
  std::vector<boost::shared_ptr<SparseMatrix> > oparray = opcomp.get_element(I,J);
  for(int i=0; i<oparray.size(); i++) {
    world.recv(processorindex(trimap_2d(I, J, dmrginp.last_site())), optype+i*10+1000*J+100000*I, *oparray[i]);
  }
#endif
}

void SpinBlock::addAdditionalCompOps()
{
#ifndef SERIAL
  boost::mpi::communicator world;
  if (world.size() == 1)
    return; //there is no need to have additional compops

  int length = dmrginp.last_site();

  if (!ops[CRE]->is_local()) {
    for(int i=0; i<get_sites().size(); i++) {
      if (ops[CRE]->has(sites[i])) {
        if (processorindex(sites[i]) != mpigetrank()) ops[CRE]->add_local_indices(sites[i]);
        ops[CRE]->set_local() = true;
        mpi::broadcast(world, *(ops[CRE]->get_element(sites[i])[0]), processorindex(sites[i]));
      }
    }
  }

  if (has(DES)) {
    if (!ops[DES]->is_local()) {
      for(int i=0; i<get_sites().size(); i++) {
        if (ops[DES]->has(sites[i])) {
        if (processorindex(sites[i]) != mpigetrank()) ops[DES]->add_local_indices(sites[i]);
        ops[DES]->set_local() = true;
        mpi::broadcast(world, *(ops[DES]->get_element(sites[i])[0]), processorindex(sites[i]));
	      }
      }
    }
  }
  if(dmrginp.calc_type() == MPS_NEVPT)
  {
  //  if ( has(CDD_CRE_DESCOMP))
  //    if (!ops[CDD_CRE_DESCOMP]->is_local) {
  //      for(int i=0; i<get_sites().size(); i++) {
  //        if (ops[CDD_CRE_DESCOMP]->has(sites[i])) {
  //        if (processorindex(sites[i]) != mpigetrank()) ops[CDD_CRE_DESCOMP]->add_local_indices(sites[i]);
  //          ops[CDD_CRE_DESCOMP]->set_local() = true;
  //          mpi::broadcast(world, *(ops[CDD_CRE_DESCOMP]->get_element(sites[i])[0]), processorindex(sites[i]));
  //        }
  //      }
  //    }
    return ;
  }



  vector<int> dotindice;
  dotindice.push_back((sites[0] == 0) ? complementary_sites[0] : complementary_sites[complementary_sites.size()-1]);
  if (!dmrginp.spinAdapted()) { // when non-spinadapted, sites are spin orbitals
    dotindice.push_back((sites[0] == 0) ? complementary_sites[1] : complementary_sites[complementary_sites.size()-2]);    
  }

  for (int idx = 0; idx < dotindice.size(); ++idx) {
    int dotopindex = dotindice[idx];

    for (int i=0; i<complementary_sites.size(); i++) {
      int compsite = complementary_sites[i];
      if (std::find(dotindice.begin(), dotindice.end(), compsite) != dotindice.end())
        continue;
      int I = (compsite > dotopindex) ? compsite : dotopindex;
      int J = (compsite > dotopindex) ? dotopindex : compsite;
      if (processorindex(compsite) == processorindex(trimap_2d(I, J, length)))
        continue;
      if (processorindex(compsite) == mpigetrank()) {
        //this will potentially receive some ops        
        bool other_proc_has_ops = true;
        world.recv(processorindex(trimap_2d(I, J, length)), 0, other_proc_has_ops);
        if (other_proc_has_ops) {
	      ops[CRE_DESCOMP]->add_local_indices(I, J);
	      recvcompOps(*ops[CRE_DESCOMP], I, J, CRE_DESCOMP);
        }
        other_proc_has_ops = true;
        world.recv(processorindex(trimap_2d(I, J, length)), 0, other_proc_has_ops);
        if (other_proc_has_ops) {
	      ops[DES_DESCOMP]->add_local_indices(I, J);
	      recvcompOps(*ops[DES_DESCOMP], I, J, DES_DESCOMP);
        }
	if (has(DES)) {
	  other_proc_has_ops = true;
	  world.recv(processorindex(trimap_2d(I, J, length)), 0, other_proc_has_ops);
	  if (other_proc_has_ops) {
	    ops[CRE_CRECOMP]->add_local_indices(I, J);
	    recvcompOps(*ops[CRE_CRECOMP], I, J, CRE_CRECOMP);
	  }
	  other_proc_has_ops = true;
	  world.recv(processorindex(trimap_2d(I, J, length)), 0, other_proc_has_ops);
	  if (other_proc_has_ops) {
	    ops[DES_CRECOMP]->add_local_indices(I, J);
	    recvcompOps(*ops[DES_CRECOMP], I, J, DES_CRECOMP);
	  }
	}
      } 
      else {
        //this will potentially send some ops
        if (processorindex(trimap_2d(I, J, length)) == mpigetrank()) {
	  bool this_proc_has_ops = ops[CRE_DESCOMP]->has_local_index(I, J);
	  world.send(processorindex(compsite), 0, this_proc_has_ops);
	  if (this_proc_has_ops) {
	    sendcompOps(*ops[CRE_DESCOMP], I, J, CRE_DESCOMP, compsite);
	  }
          this_proc_has_ops = ops[DES_DESCOMP]->has_local_index(I, J);
	  world.send(processorindex(compsite), 0, this_proc_has_ops);
	  if (this_proc_has_ops) {
	    sendcompOps(*ops[DES_DESCOMP], I, J, DES_DESCOMP, compsite);     
	  }
	  if (has(DES)) {
	    this_proc_has_ops = ops[CRE_CRECOMP]->has_local_index(I, J);
	    world.send(processorindex(compsite), 0, this_proc_has_ops);
	    if (this_proc_has_ops) {
	      sendcompOps(*ops[CRE_CRECOMP], I, J, CRE_CRECOMP, compsite);     
	    }
	    this_proc_has_ops = ops[DES_CRECOMP]->has_local_index(I, J);
	    world.send(processorindex(compsite), 0, this_proc_has_ops);
	    if (this_proc_has_ops) {
	      sendcompOps(*ops[DES_CRECOMP], I, J, DES_CRECOMP, compsite);     
	    }
	  }
	  
        } 
	else 
	  continue;
      }
    }
  }
#endif
}

void SpinBlock::transform_operators(std::vector<Matrix>& rotateMatrix)
{

  StateInfo oldStateInfo = braStateInfo;
  std::vector<SpinQuantum> newQuanta;
  std::vector<int> newQuantaStates;
  std::vector<int> newQuantaMap;
  for (int Q = 0; Q < rotateMatrix.size (); ++Q)
  {
    if (rotateMatrix [Q].Ncols () != 0)
      {
	newQuanta.push_back (braStateInfo.quanta [Q]);
	newQuantaStates.push_back (rotateMatrix [Q].Ncols ());
	newQuantaMap.push_back (Q);
      }
  }
  StateInfo newStateInfo = StateInfo (newQuanta, newQuantaStates, newQuantaMap);

  build_and_renormalise_operators( rotateMatrix, &newStateInfo );

  braStateInfo = newStateInfo;
  braStateInfo.AllocatePreviousStateInfo ();
  *braStateInfo.previousStateInfo = oldStateInfo;
  ketStateInfo = braStateInfo;

  for (int i = 0; i < newQuantaMap.size (); ++i) {
    assert (braStateInfo.quanta [i] == oldStateInfo.quanta [newQuantaMap [i]]);
    assert (ketStateInfo.quanta [i] == oldStateInfo.quanta [newQuantaMap [i]]);
  }

  tcpu = globaltimer.totalcputime(); twall= globaltimer.totalwalltime();
  ecpu = globaltimer.elapsedcputime(); ewall= globaltimer.elapsedwalltime();
  p3out << "\t\t\t total elapsed time " <<  twall << " " << tcpu << " ... " << ewall << " " << ecpu << endl;
  p1out << "\t\t\t Transforming to new basis " << endl;
  Timer transformtimer;

  renormalise_transform( rotateMatrix, &this->braStateInfo );

  for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
    if (! it->second->is_core())
      ops[it->first]->set_core(true);

  this->direct = false;
  ecpu = transformtimer.elapsedcputime();ewall= transformtimer.elapsedwalltime();
  p3out << "\t\t\t transform time " << ewall << " " << ecpu << endl;

  if (leftBlock)
    leftBlock->clear();
  if (rightBlock)
    rightBlock->clear();


}

  void SpinBlock::transform_operators(std::vector<Matrix>& leftrotateMatrix, std::vector<Matrix>& rightrotateMatrix, bool clearRightBlock, bool clearLeftBlock)
{

  StateInfo oldbraStateInfo=braStateInfo, oldketStateInfo=ketStateInfo;
  StateInfo newbraStateInfo, newketStateInfo;
  StateInfo::transform_state(leftrotateMatrix, braStateInfo, newbraStateInfo);
  StateInfo::transform_state(rightrotateMatrix, ketStateInfo, newketStateInfo);


  build_and_renormalise_operators( leftrotateMatrix, &newbraStateInfo, rightrotateMatrix, &newketStateInfo );

  braStateInfo = newbraStateInfo;
  braStateInfo.AllocatePreviousStateInfo ();
  *braStateInfo.previousStateInfo = oldbraStateInfo;

  ketStateInfo = newketStateInfo;
  ketStateInfo.AllocatePreviousStateInfo ();
  *ketStateInfo.previousStateInfo = oldketStateInfo;


  tcpu = globaltimer.totalcputime(); twall= globaltimer.totalwalltime();
  ecpu = globaltimer.elapsedcputime(); ewall= globaltimer.elapsedwalltime();
  p3out << "\t\t\t total elapsed time " << twall << " " << tcpu << " ... " << ewall << " " << ecpu << endl;
  p1out << "\t\t\t Transforming to new basis " << endl;

  Timer transformtimer;

  renormalise_transform( leftrotateMatrix, &this->braStateInfo, rightrotateMatrix, &this->ketStateInfo );

  for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
    if (! it->second->is_core())
      ops[it->first]->set_core(true);

  this->direct = false;
  ecpu = transformtimer.elapsedcputime(); ewall= transformtimer.elapsedwalltime();
  p3out << "\t\t\t transform time " << ewall << " " << ecpu << endl;

  if (leftBlock && clearLeftBlock)
    leftBlock->clear();
  if (rightBlock && clearRightBlock)
    rightBlock->clear();


}

}
