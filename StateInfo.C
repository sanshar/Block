/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "StateInfo.h"
#include <boost/format.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>




SpinAdapted::StateInfo::StateInfo () :  leftStateInfo (0), rightStateInfo (0), totalStates (0), initialised (false), hasAllocatedMemory (false), hasCollectedQuanta (false), hasPreviousStateInfo (false) {}

SpinAdapted::StateInfo::StateInfo (const int n, const SpinQuantum q [], const int qS []) : leftStateInfo (0), rightStateInfo (0), totalStates (0), initialised (true), 
  hasAllocatedMemory (false), hasCollectedQuanta (false), hasPreviousStateInfo (false)
{
  quanta.resize (n);
  quantaStates.resize (n);
  for (int i = 0; i < n; ++i)
    {
      quanta [i] = q [i];
      quantaStates [i] = qS [i];
      totalStates += qS [i];
    }
  UnBlockIndex ();
}

SpinAdapted::StateInfo::StateInfo (const std::vector<SpinQuantum>& q, const std::vector<int>& qS, const std::vector<int>& nMap) : leftStateInfo (0), rightStateInfo (0), totalStates (0), initialised (true), hasAllocatedMemory (false), hasCollectedQuanta (false), hasPreviousStateInfo (false)
{
  quanta = q;
  quantaStates = qS;
  newQuantaMap = nMap;
  totalStates = accumulate (quantaStates.begin (), quantaStates.end (), 0);
  UnBlockIndex ();
}


SpinAdapted::StateInfo::StateInfo (const std::vector< Csf >& dets, bool addWavefunctionQuanta) : totalStates (0), initialised (true), hasAllocatedMemory (false), hasCollectedQuanta (false), hasPreviousStateInfo (false)
{
  std::vector<SpinQuantum> tmp;
  for (int i = 0; i < dets.size (); ++i) {
    SpinQuantum s2 = SpinQuantum ( (dets [i]).n, (dets [i]).S, (dets[i]).sym_is());
    tmp.push_back(s2);
  }
  sort (tmp.begin (), tmp.end ());
  unique_copy (tmp.begin (), tmp.end (), back_inserter (quanta));
  sort (quanta.begin (), quanta.end ());

  for (int i = 0; i < quanta.size (); ++i)
  {
    int cnt = 0;
    cnt = count (tmp.begin (), tmp.end (), quanta [i]);
    quantaStates.push_back (cnt);
  }

  totalStates = accumulate (quantaStates.begin (), quantaStates.end (), 0);
  UnBlockIndex ();
}




void SpinAdapted::StateInfo::UnBlockIndex ()
{
  unBlockedIndex.resize (quanta.size ());
  for (int i = 0; i < unBlockedIndex.size (); ++i)
    unBlockedIndex [i] = (i > 0) ? unBlockedIndex [i - 1] + quantaStates [i - 1] : 0; 
}

namespace SpinAdapted {
ostream& operator<< (ostream& os, const StateInfo& s) {
  os <<"\t\t\t Total number of Quanta and States : "<< s.quanta.size()<<" "<<s.totalStates<<endl;
    os <<  "\t\t\t States per Quantum :: \t";
    for (int i = 0; i < s.quanta.size(); ++i) os << s.quanta[i] << " =  "<<s.quantaStates[i]<<", ";
    os << endl;
  return os;
}


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



void TensorProduct (StateInfo& a, StateInfo& b, const SpinQuantum q, const int constraint, StateInfo& c, StateInfo* compState) {
  c.leftStateInfo = &a;
  c.rightStateInfo = &b;
  c.allowedQuanta.ReSize (a.quanta.size (), b.quanta.size ());
  c.quantaMap.ReSize (a.quanta.size (), b.quanta.size ());
  c.totalStates = 0;
  c.initialised = true;

  //if non-conventional nevpt2 is invoked, states with up to 2 particles more 
  //are required to be generated.
  bool DoNonConvNevpt2 = ((dmrginp.nevpt2())&&(!dmrginp.read_higherpdm()));
  
  if (DoNonConvNevpt2){
    for (int i = 0; i < a.quanta.size (); ++i)
      for (int j = 0; j < b.quanta.size (); ++j)
        if (  constraint == LessThanQ && 
              (  a.quanta [i].get_n() + b.quanta [j].get_n() > q.get_n() ))

          {
          if ((  a.quanta [i].get_n() + b.quanta [j].get_n() > q.get_n() + 2 ))continue;
          else{
              vector<SpinQuantum> v = a.quanta[i] + b.quanta[j];
              for (int vq=0; vq< v.size(); vq++) {
                bool include = false;
                if (compState != 0) {
                  for (int k=0; k<compState->quanta.size(); k++)
                  {
                    if (q.allow(v[vq], compState->quanta[k])){
                      include = true;
                      break;
                    }
                  }
                }
                else
                  include = true;
                if (!include) continue;
                c.quanta.push_back(v[vq]);
                c.quantaStates.push_back(a.quantaStates[i] * b.quantaStates[j]);
                c.totalStates += a.quantaStates[i]*b.quantaStates[j];
                c.allowedQuanta(i,j) = true;
                c.quantaMap(i,j).push_back(c.quanta.size() - 1);
                c.leftUnMapQuanta.push_back(i);
                c.rightUnMapQuanta.push_back(j);
              }//vk
            //continue;
            }
          }
        else if (constraint == EqualQ && q.allow(a.quanta [i] , b.quanta [j]) )
        {
          c.quantaMap(i,j).resize(1);

          c.quanta.push_back(q);
          c.quantaStates.push_back(a.quantaStates[i] * b.quantaStates[j]);
          c.totalStates += a.quantaStates[i]*b.quantaStates[j];
          c.allowedQuanta(i,j) = true;
          c.quantaMap(i,j)[0] = c.quanta.size() - 1;
          c.leftUnMapQuanta.push_back(i);
          c.rightUnMapQuanta.push_back(j);
        }
        else if (constraint == LessThanQ)
        {
          vector<SpinQuantum> v = a.quanta[i] + b.quanta[j];
          for (int vq=0; vq< v.size(); vq++) {
            if ( (v[vq].get_n() > q.get_n()+2)
                 || ( abs(v[vq].get_s().getirrep()-q.get_s().getirrep()) > (q.get_n()+2-v[vq].get_n())  ))
                 continue;
            bool include = false;
            if (compState != 0) {
              for (int k=0; k<compState->quanta.size(); k++)
              {
                if (q.allow(v[vq], compState->quanta[k])){
                  include = true;
                  break;
                }
              }
            }
            else
              include = true;
            if (!include) continue;
            c.quanta.push_back(v[vq]);
            c.quantaStates.push_back(a.quantaStates[i] * b.quantaStates[j]);
            c.totalStates += a.quantaStates[i]*b.quantaStates[j];
            c.allowedQuanta(i,j) = true;
            c.quantaMap(i,j).push_back(c.quanta.size() - 1);
            c.leftUnMapQuanta.push_back(i);
            c.rightUnMapQuanta.push_back(j);
            }
        }
  }//non-conventional nevpt2 invoked
  else{
    for (int i = 0; i < a.quanta.size (); ++i)
      for (int j = 0; j < b.quanta.size (); ++j)
        if (constraint == LessThanQ && 
                (a.quanta [i].get_n() + b.quanta [j].get_n() > q.get_n() )) {	  
              continue;
            } else if (constraint == EqualQ && q.allow(a.quanta [i] , b.quanta [j]) ) {
              c.quantaMap(i,j).resize(1);
              c.quanta.push_back(q);
              c.quantaStates.push_back(a.quantaStates[i] * b.quantaStates[j]);
              c.totalStates += a.quantaStates[i]*b.quantaStates[j];
              c.allowedQuanta(i,j) = true;
              c.quantaMap(i,j)[0] = c.quanta.size() - 1;
              c.leftUnMapQuanta.push_back(i);
              c.rightUnMapQuanta.push_back(j);
        } 
        else if (constraint == EqualS) {
          vector<SpinQuantum> v = a.quanta[i] + b.quanta[j];
          for (int vq=0; vq<v.size(); vq++) {
            if (v[vq].get_n() > q.get_n() || v[vq].get_s() != q.get_s() || v[vq].get_s().getirrep() != q.get_s().getirrep() || v[vq].get_symm() != q.get_symm())
              continue;

            c.quanta.push_back(v[vq]);
                c.quantaStates.push_back(a.quantaStates[i] * b.quantaStates[j]);
                c.totalStates += a.quantaStates[i]*b.quantaStates[j];
                c.allowedQuanta(i,j) = true;
                c.quantaMap(i,j).push_back(c.quanta.size() - 1);
                c.leftUnMapQuanta.push_back(i);
                c.rightUnMapQuanta.push_back(j);
          }
        } 
        else if (constraint == LessThanN) {
          vector<SpinQuantum> v = a.quanta[i] + b.quanta[j];
          for (int vq=0; vq<v.size(); vq++) {
            if (v[vq].get_n() > q.get_n())
              continue;

            c.quanta.push_back(v[vq]);
            c.quantaStates.push_back(a.quantaStates[i] * b.quantaStates[j]);
            c.totalStates += a.quantaStates[i]*b.quantaStates[j];
            c.allowedQuanta(i,j) = true;
            c.quantaMap(i,j).push_back(c.quanta.size() - 1);
            c.leftUnMapQuanta.push_back(i);
            c.rightUnMapQuanta.push_back(j);
          }
        } 
        else if (constraint == LessThanQ) {
              vector<SpinQuantum> v = a.quanta[i] + b.quanta[j];
              for (int vq=0; vq< v.size(); vq++) {
                if ( (v[vq].get_n() > q.get_n()) || (v[vq].get_n()==q.get_n() && v[vq].get_s() != q.get_s())
                  || ( abs(v[vq].get_s().getirrep()-q.get_s().getirrep()) > (q.get_n()-v[vq].get_n())  )
                  || (v[vq].get_n() == q.get_n() && v[vq].get_symm() != q.get_symm()) )
                  continue;
                bool include = false;
                if (compState != 0) {
                  for (int k=0; k<compState->quanta.size(); k++) {
                    if (q.allow(v[vq], compState->quanta[k])){
                          include = true;
                          break;
                    }
                  }
                } else
                  include = true;
                if (!include)
              continue;
                c.quanta.push_back(v[vq]);
                c.quantaStates.push_back(a.quantaStates[i] * b.quantaStates[j]);
                c.totalStates += a.quantaStates[i]*b.quantaStates[j];
                c.allowedQuanta(i,j) = true;
                c.quantaMap(i,j).push_back(c.quanta.size() - 1);
                c.leftUnMapQuanta.push_back(i);
                c.rightUnMapQuanta.push_back(j);
              }
        }
  }//regular DMRG
  c.UnBlockIndex ();
}


void TensorProduct (StateInfo& a, StateInfo& b, StateInfo& c, const int constraint, StateInfo* compState)
{
  ObjectMatrix<char> dummy;
  assert (constraint != WITH_LIST);


  if (constraint == NO_PARTICLE_SPIN_NUMBER_CONSTRAINT) {
    TensorProduct (a, b, dmrginp.effective_molecule_quantum(), LessThanQ, c, compState);
  } else if (constraint == PARTICLE_SPIN_NUMBER_CONSTRAINT) {
    TensorProduct (a, b, dmrginp.effective_molecule_quantum(), EqualQ, c);
  } else if (constraint == PARTICLE_NUMBER_CONSTRAINT) {
    TensorProduct (a, b, dmrginp.effective_molecule_quantum(), LessThanN, c);
  } else if (constraint == SPIN_NUMBER_CONSTRAINT) {
    TensorProduct (a, b, dmrginp.effective_molecule_quantum(), EqualS, c);
  }
}

}

void SpinAdapted::StateInfo::CollectQuanta ()
{
  StateInfo uniqueStateInfo;
  uniqueStateInfo.initialised = true;
  uniqueStateInfo.leftStateInfo = leftStateInfo;
  uniqueStateInfo.rightStateInfo = rightStateInfo;
  uniqueStateInfo.leftUnMapQuanta = leftUnMapQuanta;
  uniqueStateInfo.rightUnMapQuanta = rightUnMapQuanta;

  std::vector<SpinQuantum> duplicateQuanta = quanta;
  sort (duplicateQuanta.begin (), duplicateQuanta.end ());
  unique_copy (duplicateQuanta.begin (), duplicateQuanta.end (), back_inserter (uniqueStateInfo.quanta));
  sort (uniqueStateInfo.quanta.begin (), uniqueStateInfo.quanta.end ());

  uniqueStateInfo.quantaStates.resize (uniqueStateInfo.quanta.size ());
  uniqueStateInfo.oldToNewState.resize (uniqueStateInfo.quanta.size ());
  
  for (int i = 0; i < uniqueStateInfo.quanta.size (); ++i)
    for (int j = 0; j < quanta.size (); ++j)
      if (quanta [j] == uniqueStateInfo.quanta [i])
	{
	  uniqueStateInfo.quantaStates [i] += quantaStates [j];
	  uniqueStateInfo.oldToNewState [i].push_back (j);
	}
  
  uniqueStateInfo.totalStates = accumulate (uniqueStateInfo.quantaStates.begin (), uniqueStateInfo.quantaStates.end (), 0);
  uniqueStateInfo.AllocateUnCollectedStateInfo ();
  *uniqueStateInfo.unCollectedStateInfo = *this;
  *this = uniqueStateInfo;
  UnBlockIndex ();
}

void SpinAdapted::StateInfo::AllocateUnCollectedStateInfo ()
{
  unCollectedStateInfo = boost::shared_ptr<StateInfo> (new StateInfo);
  hasCollectedQuanta = true;
}

void SpinAdapted::StateInfo::quanta_distribution (std::vector<SpinQuantum>& qnumbers, std::vector<int>& distribution, const bool complement)
{
  // first extract the unique quantum numbers

  std::vector<SpinQuantum> nosymquanta;
  std::vector<int> nosymquantastates;

  nosymquanta = quanta;
  nosymquantastates = quantaStates;
  

  qnumbers.resize (0);
  distribution.resize (0);
  
  for (int i = 0; i < nosymquanta.size (); ++i) {
    if (complement) {
      if (SpinQuantum::can_complement (nosymquanta [i])) {
        std::vector<SpinQuantum> complements = nosymquanta[i].get_complement();
        for (int j=0; j<complements.size(); j++) {
          qnumbers.push_back (complements[j]);
          distribution.push_back (nosymquantastates [i]);
        }
        //qnumbers.rbegin ()->complementize ();
      }
    } else {
      qnumbers.push_back (nosymquanta [i]);
      distribution.push_back (nosymquantastates [i]);
    }
  }
}

void SpinAdapted::StateInfo::Allocate ()
{
  hasAllocatedMemory = true;
  leftStateInfo = new StateInfo;
  leftStateInfo->leftStateInfo = new StateInfo;
  leftStateInfo->rightStateInfo = new StateInfo;
  rightStateInfo = new StateInfo;
  rightStateInfo->leftStateInfo = new StateInfo;
  rightStateInfo->rightStateInfo = new StateInfo;
}

void SpinAdapted::StateInfo::Free ()
{
  hasAllocatedMemory = false;
  delete (leftStateInfo->leftStateInfo);
  delete (leftStateInfo->rightStateInfo);
  delete (leftStateInfo);
  delete (rightStateInfo->leftStateInfo);
  delete (rightStateInfo->rightStateInfo);
  delete (rightStateInfo);
}

void SpinAdapted::StateInfo::AllocatePreviousStateInfo ()
{
  previousStateInfo = boost::shared_ptr<StateInfo> (new StateInfo);
  hasPreviousStateInfo = true;
}

void SpinAdapted::StateInfo::UnMapQuantumState (const int QS, const int secondQSTotal, int& firstQS, int& secondQS) const
{
  firstQS = QS / secondQSTotal;
  secondQS = QS % secondQSTotal;
}

//takes the old state and a rotaition matrix and gives back a rotated renormalized state
void SpinAdapted::StateInfo::transform_state(const std::vector<Matrix>& rotateMatrix, const StateInfo& stateInfo, StateInfo& newStateInfo)
{
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
  newStateInfo = StateInfo (newQuanta, newQuantaStates, newQuantaMap);
}

void SpinAdapted::StateInfo::restore(bool forward, const vector<int>& sites, StateInfo& stateInfo, int left)
{
  std::string file;
  int first = min(sites[0], *sites.rbegin()), last = max(sites[0], *sites.rbegin());
  if (forward)
    file = str(boost::format("%s%s%d%s%d%s%d%s%d%s") % dmrginp.load_prefix() % "/StateInfo-forward-" % first % "-" % last % "." % left % "." % mpigetrank() % ".tmp" );
  else
    file = str(boost::format("%s%s%d%s%d%s%d%s%d%s") % dmrginp.load_prefix() % "/StateInfo-backward-" % first % "-" % last % "." % left % "." % mpigetrank() % ".tmp" );

  p1out << "\t\t\t Loading state file :: " << file << endl;

  std::ifstream ifs(file.c_str(), std::ios::binary);
  boost::archive::binary_iarchive load_state(ifs);

  load_state >> stateInfo;
  
  ifs.close();

}

void SpinAdapted::StateInfo::store(bool forward, const vector<int>& sites, StateInfo& stateInfo, int left)
{
  
  std::string file;
  int first = min(sites[0], *sites.rbegin()), last = max(sites[0], *sites.rbegin());
  if (dmrginp.spinAdapted()) {
    if (forward)
      file = str(boost::format("%s%s%d%s%d%s%d%s%d%s") % dmrginp.load_prefix() % "/StateInfo-forward-" % first % "-" % last % "." % left % "." % mpigetrank() % ".tmp" );
    else
      file = str(boost::format("%s%s%d%s%d%s%d%s%d%s") % dmrginp.load_prefix() % "/StateInfo-backward-" % first % "-" % last % "." % left % "." % mpigetrank() % ".tmp" );
  }
  else {
    if (forward)
      file = str(boost::format("%s%s%d%s%d%s%d%s%d%s") % dmrginp.load_prefix() % "/StateInfo-forward-" % (first/2) % "-" % (last/2) % "." % left % "." % mpigetrank() % ".tmp" );
    else
      file = str(boost::format("%s%s%d%s%d%s%d%s%d%s") % dmrginp.load_prefix() % "/StateInfo-backward-" % (first/2) % "-" % (last/2) % "." % left % "." % mpigetrank() % ".tmp" );
  }
  
  p1out << "\t\t\t Saving state file :: " << file << endl;

  std::ofstream ofs(file.c_str(), std::ios::binary);
  boost::archive::binary_oarchive save_state(ofs);
  
  save_state << stateInfo;
  
  ofs.close();

}

