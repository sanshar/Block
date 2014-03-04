/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_STATE_INFO_HEADER
#define SPIN_STATE_INFO_HEADER
#include "ObjectMatrix.h"
#include "SpinQuantum.h"
#include "csf.h"
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>

namespace SpinAdapted{
enum {
  NO_PARTICLE_SPIN_NUMBER_CONSTRAINT,
  PARTICLE_SPIN_NUMBER_CONSTRAINT,
  SPIN_NUMBER_CONSTRAINT
};

enum {
  AnyQ,  
  LessThanQ, 
  EqualQ,
  EqualS,
  WITH_LIST /*< states are added together if the are allowed by the quantaList */
};


class StateInfo
{
 private:
  friend class boost::serialization::access;
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
    ar & totalStates & initialised & unBlockedIndex & quanta & quantaStates
       & leftUnMapQuanta & rightUnMapQuanta & allowedQuanta
      & quantaMap & oldToNewState & hasCollectedQuanta & newQuantaMap;
    if (hasCollectedQuanta)
      ar & unCollectedStateInfo;
    ar & hasPreviousStateInfo;
    if (hasPreviousStateInfo)
      ar & previousStateInfo;
    ar & hasAllocatedMemory;
    }

 public:
  bool hasAllocatedMemory; 
  bool hasCollectedQuanta; 
  bool hasPreviousStateInfo;

  StateInfo* leftStateInfo;
  StateInfo* rightStateInfo;
  
  boost::shared_ptr<StateInfo> unCollectedStateInfo;
  boost::shared_ptr<StateInfo> previousStateInfo;

  std::vector<SpinQuantum> quanta; //the quantas present
  std::vector<int> quantaStates;  //the number of each quanta
  ObjectMatrix<char> allowedQuanta; //tensor product of left and right stateinfo states

  //takes i and j quanta of the left and the right state and gives a vector of all resulting quantas. (for abelian symmetry the vector would be of length 1)
  ObjectMatrix< std::vector<int> > quantaMap; 

  //for the uncollectedstateinfo each quanta is made up of a left and a right quanta, these two vectors unmap an uncollected quata back to the left and the right quanta.
  std::vector<int> leftUnMapQuanta;  
  std::vector<int> rightUnMapQuanta;

  int totalStates;

  //the next three are concerned with blocking and unblocking the states
  std::vector<int> unBlockedIndex;
  std::vector< std::vector<int> > oldToNewState; //quanta[I] is the same as quanta[Ij] where Ij are the indices inside the Ith vector.
  std::vector<int> newQuantaMap;
  bool initialised;
 public:
  StateInfo();
  StateInfo (const int n, const SpinQuantum q [], const int qS []); 
  StateInfo(const std::vector< Csf >& dets, bool addWavefunctionQuanta = false);
  StateInfo (const std::vector<SpinQuantum>& q, const std::vector<int>& qS, const std::vector<int>& nMap);
  void UnBlockIndex();

  friend ostream& operator<<(ostream& os, const StateInfo& s);

  friend void TensorProduct (StateInfo& a, StateInfo& b, const SpinQuantum q, const int constraint, StateInfo& c, StateInfo* compState=0);
  // interface to the above function
  friend void TensorProduct (StateInfo& a, StateInfo& b, StateInfo& c, const int constraint, StateInfo* compState=0);
  void quanta_distribution (std::vector<SpinQuantum>& qnumbers, std::vector<int>& distribution, const bool complement);

  void Allocate();
  void Free();
  void AllocatePreviousStateInfo ();
  void CollectQuanta();
  void AllocateUnCollectedStateInfo ();
  int getquantastates(int i) {return quantaStates.at(i);}
  int getquantastates(int i) const {return quantaStates.at(i);}
  void UnMapQuantumState (const int QS, const int secondQSTotal, int& firstQS, int& secondQS) const;
  static void restore(bool forward, const vector<int>& sites, StateInfo& states, int left, int right);
  static void store(bool forward, const vector<int>& sites, StateInfo& states, int left, int right);
  static void transform_state(std::vector<Matrix>& rotateMatrix, StateInfo& stateInfo, StateInfo& newStateInfo);
};
}

#endif
