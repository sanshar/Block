/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is free software: you can redistribute it and/or modify         
it under the terms of the GNU General Public License as published by         
the Free Software Foundation, either version 3 of the License, or            
(at your option) any later version.                                          
                                                                             
This program is distributed in the hope that it will be useful,              
but WITHOUT ANY WARRANTY; without even the implied warranty of               
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
GNU General Public License for more details.                                 
                                                                             
You should have received a copy of the GNU General Public License            
along with this program.  If not, see <http://www.gnu.org/licenses/>.        

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
  PARTICLE_SPIN_NUMBER_CONSTRAINT
};

enum {
  AnyQ,  
  LessThanQ, 
  EqualQ, 
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
  ObjectMatrix<char> allowedQuanta; //
  ObjectMatrix< std::vector<int> > quantaMap;
  std::vector<int> leftUnMapQuanta;
  std::vector<int> rightUnMapQuanta;

  int totalStates;
  std::vector<int> unBlockedIndex;
  std::vector< std::vector<int> > oldToNewState;
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
  static void store(bool forward, const vector<int>& sites, const vector<StateInfo>& states);
  static void restore(bool forward, const vector<int>& sites, vector<StateInfo>& states);

};
}

#endif
