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
  SPIN_NUMBER_CONSTRAINT,
  PARTICLE_NUMBER_CONSTRAINT
};

enum {
  LessThanQ, 
  EqualQ,
  EqualS,
  LessThanN,
  AnyQ,  
  WITH_LIST /*< states are added together if the are allowed by the quantaList */
};


/// @brief StateInfo carries the groups of quantum numbers (and their dimensions) formed
///        for a block of sites. This is the "information" on the states for the block.
/// 
///
/// A StateInfo object can arise in multiple ways. Different indexing information is needed
/// in these different situations and this is the reason for the large number
/// member variables. Roughly speaking, there are StateInfo arising from
/// two situations: when there is a product Hilbert space, and with no product Hilbert space.
/// When there is a product structure, StateInfo carries information on how to go from 
/// quantum numbers ql, qr, to the compound quantum number (ql+qr), and backwards. There
/// can be various constraints on ql+qr (e.g. a total particle number constraint). The
/// different kinds of StateInfo are:
///
/// 1. StateInfo from a set of sites, where the states have undergone no renormalization.
///    This might be the StateInfo for a single site. In this case, we always use a state
///    of the site, i.e. its full Hilbert Space. We do not use the product structure
///    of the Hilbert space.
/// 2. StateInfo from combining two blocks (referred to as left and right)  
///    This StateInfo defines information for the LR product structure
///    of the Hilbert space:  `leftStateInfo` and `rightStateInfo` (pointers), as well as
///    `leftUnMapQuanta`, `rightUnMapQuanta`, `allowedQuanta`. This type of StateInfo
///    is often constructed with "constraints", i.e. only certain combinations
///    of ql+qr are allowed (e.g. constraints on total particle number).
///    Constructed using `TensorProduct` function.
/// 3. StateInfo from combining two blocks, and with collection.
///    This defines `unCollectedStateInfo`, which points to the StateInfo before collection,
///    and `oldToNewState` which maps the collected and uncollected StateInfo.
///    Convert from 2. to 3. using `CollectQuanta()`.
/// 4. StateInfo after combined blocks have been renormalized.
///    This has no LR product structure. `newQuantaMap` maps indices between the
///    states before and after renormalization (when some states are thrown out).
///    Made using `transform_state`.

 class StateInfo;
 
  void TensorProduct (StateInfo& a, StateInfo& b, const SpinQuantum q, const int constraint, StateInfo& c, StateInfo* compState=0);

  /// Interface to other TensorProduct function.
  void TensorProduct (StateInfo& a, StateInfo& b, StateInfo& c, const int constraint, StateInfo* compState=0);

 
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

  /// For "product" StateInfo's (type 2., 3. above) the StateInfo of the left Block.
  StateInfo* leftStateInfo;

  /// For "product" StateInfo's (type 2., 3. above) the StateInfo of the right Block.
  StateInfo* rightStateInfo;
  
  /// For "product StateInfo", after `CollectQuanta` is called (type 3. above), points to StateInfo before call.
  boost::shared_ptr<StateInfo> unCollectedStateInfo;

  boost::shared_ptr<StateInfo> previousStateInfo;

  /// Quantum numbers for a set of sites.
  std::vector<SpinQuantum> quanta; 

  /// Number of states per quantum number.
  std::vector<int> quantaStates;  

  /// @brief Truth matrix of left and right stateinfo state (for type 2. StateInfo).
  /// `allowedQuanta(i,j)==true` if qi+qj satisfies constraints at time of type 2. StateInfo construction.
  ObjectMatrix<char> allowedQuanta;

  /// @brief Maps compound left, right quanta (ql,qr), to index of ql+qr. 
  /// (For non-Abelian
  /// quantum numbers, ql+qr gives a vector of quantum numbers, thus a vector of indices).
  ObjectMatrix< std::vector<int> > quantaMap; 

  /// @brief Quanta on left block, used during product StateInfo construction.
  /// For a "product StateInfo" (type 2.) each quantum is a sum of left and right quanta (ql, qr). 
  /// However, not all ql and qr contribute, because certain combinations (ql+qr) violate the
  /// constraints during StateInfo construction (e.g. total particle or spin constraints). This
  /// vector contains ql, where ql+qr did not violate the constraints.
  std::vector<int> leftUnMapQuanta;  

  /// @brief Quanta on right block, used during product StateInfo construction.
  /// qr, where ql+qr did not violate constraint. See documentation for `leftUnMapQuanta`.
  std::vector<int> rightUnMapQuanta;

  int totalStates;

  //the next three are concerned with blocking and unblocking the states
  std::vector<int> unBlockedIndex;
  std::vector< std::vector<int> > oldToNewState; //quanta[I] is the same as quanta[Ij] where Ij are the indices inside the Ith vector.
  std::vector<int> newQuantaMap;

  /// flag for StateInfo in valid state
  bool initialised;
 public:
  StateInfo();
  StateInfo (const int n, const SpinQuantum q [], const int qS []); 

  /// Type 1. Constructor 
  /// \param[in] CSF's on a set of sites
  StateInfo(const std::vector< Csf >& dets, bool addWavefunctionQuanta = false);


  /// Type 1. Constructor 
  /// \param[in] explicit quantum numbers, states, newQuantaMap. 
  /// 
  /// Used only rarely to 
  /// construct a valid StateInfo for a block after renormalization (which provides the
  /// newQuantaMap).
  StateInfo (const std::vector<SpinQuantum>& q, const std::vector<int>& qS, const std::vector<int>& nMap);
  void UnBlockIndex();

  friend ostream& operator<<(ostream& os, const StateInfo& s);

  /// Make type 2. StateInfo, by combining a and b.
  /// \param[in] a Left StateInfo.
  /// \param[in] b Right StateInfo.
  /// \param[in[ q SpinQuantum (type 2.) 
  /// \param[in] constraint One of `AnyQ` (no constraint), `LessThanQ`, `EqualS`.
  /// \param[out] c Output (type 2.) Product StateInfo.
  /// \param[in] compState Pointer to a StateInfo for a set of complementary quantum numbers. 
  ///        Implemented only
  ///        with `LessThanQ`, only constructs a+b if a+b+compStateQ in q. Usually 0, 
  ///        used only in warmup.
  friend void TensorProduct (StateInfo& a, StateInfo& b, const SpinQuantum q, const int constraint, StateInfo& c, StateInfo* compState);

  /// Interface to other TensorProduct function.
  friend void TensorProduct (StateInfo& a, StateInfo& b, StateInfo& c, const int constraint, StateInfo* compState);

  friend void makeStateInfo(StateInfo& s, int site);

  void quanta_distribution (std::vector<SpinQuantum>& qnumbers, std::vector<int>& distribution, const bool complement);

  void Allocate();
  void Free();
  void AllocatePreviousStateInfo ();
  void CollectQuanta();
  void AllocateUnCollectedStateInfo ();
  int getquantastates(int i) {return quantaStates.at(i);}
  int getquantastates(int i) const {return quantaStates.at(i);}
  void UnMapQuantumState (const int QS, const int secondQSTotal, int& firstQS, int& secondQS) const;
  static void restore(bool forward, const vector<int>& sites, StateInfo& states, int state);
  static void store(bool forward, const vector<int>& sites, StateInfo& states, int state);

  /// Make type 4. StateInfo
  /// \param[in] rotateMatrix Rotation matrix from DMRG truncation.
  /// \param[in] stateInfo Starting StateInfo.
  /// \param[out] newStateInfo New StateInfo after DMRG truncation.
  static void transform_state(const std::vector<Matrix>& rotateMatrix, const StateInfo& stateInfo, StateInfo& newStateInfo);
};
}

#endif
