/* 
 * File:   ripdm_opconstruct.h
 * Author: roemelt
 *
 * Created on October 31, 2013, 2:18 PM
 */

#ifndef NEVPT2_OPCONSTRUCT_H
#define	NEVPT2_OPCONSTRUCT_H

#define _INITIAL_ITER_ 0
#define _REGULAR_ITER_ 1
#define _FINAL_ITER_   2

#include "nevpt2_info.h"

namespace SpinAdapted{

    void CheckOperator(SpinBlock &big, Wavefunction &WF, int t, int u, int v, SparseMatrix &Op, bool left);
    void CheckOperator(SpinBlock &big, Wavefunction &WF, int t,SparseMatrix &Op, bool left);
    void CheckOperator_(SpinBlock &big, Wavefunction &WF, int t, int u, int v, SparseMatrix &Op, bool left);

  //============================================================================
  // determine whether a given set of orbital indices is valid at this point in 
  // the sweep. This is to avoid double counting of certain index combinations
  //============================================================================
   bool CheckAllowedOperator(int t, int u, int v, int DotIndex, int iterCase, 
                             bool two_indexes_on_right=false);
   
  //============================================================================
  // Construct a+a+a where two indices are on the dot and the other on the
  // left block
  //============================================================================
  void ConstructCreCreDes_1_2_0(SpinBlock &big, SpinBlock *lBlock, SpinBlock *dotBlock,
                                ThreeIndOpArray &CCD, int *OrbWin);
  
  //============================================================================
  // Construct a+a+a where two indices are on the left and the other on the
  // dot block
  //============================================================================
  void ConstructCreCreDes_2_1_0(SpinBlock &big, SpinBlock *lBlock, SpinBlock *dotBlock,
                                ThreeIndOpArray &CCD, int *OrbWin);

  
  //============================================================================
  // Construct a+a+a on a single block
  //============================================================================
  void ConstructCreCreDesSingleBlock(SpinBlock &big, SpinBlock *block, ThreeIndOpArray &CCD,
                                     IntegralContainer &IKJL, int iterCase);

  //============================================================================
  // Construct a+a+a|psi> on two blocks block
  //============================================================================
  void ConstructCreCreDesTwoBlocks(SpinBlock &big, Wavefunction &WF, vector<WavefunctionArray> &T, 
                                 IntegralContainer &IKJL,int iterCase, int *OrbWin);
  
  //============================================================================
  //Construct the contributions to the perturber functions V(i) where all indices
  //are on one block
  //============================================================================
    void ConstructViSingleBlock(SpinBlock &big, Wavefunction &WF, ThreeIndOpArray &CCD,
                              ThreeIndOpArray CCD_, vector<WavefunctionArray> &Ti, 
                              IntegralContainer &IKJL, NEVPT2Info &Info);
  
  //============================================================================
  //Construct the one-electron part of the V(i) operator
  //============================================================================
  void ConstructCre(SpinBlock &big, Matrix &heff, Wavefunction &WF, vector<WavefunctionArray> &T);
  
  //============================================================================
  //the driver for the construction of the a+a+a|psi> operators
  //============================================================================
  void ConstructCreCreDes(SpinBlock &big, vector<Wavefunction> &WF, ThreeIndOpArray &CCD,
                          IntegralContainer &IKJL, SweepParams &sweepParams,NEVPT2Info &Info);
  
  
  
  //============================================================================
  // Construct a+aa where two indices are on the dot and the other on the
  // left block
  //============================================================================
  void ConstructCreDesDes_1_2_0(SpinBlock &big, SpinBlock *lBlock, SpinBlock *dotBlock,
                                ThreeIndOpArray &CDD, IntegralContainer &IKJA, int *OrbWin);
  
  //============================================================================
  // Construct a+aa where two indices are on the leftBlock and the other on the
  // dot Block
  //============================================================================
  void ConstructCreDesDes_2_1_0(SpinBlock &big, SpinBlock *lBlock, SpinBlock *dotBlock,
                                ThreeIndOpArray &CDD, IntegralContainer &IKJA, int *OrbWin);
  
  
  //============================================================================
  // Construct a+aa on a single block
  //============================================================================
  void ConstructCreDesDesSingleBlock(SpinBlock &big, SpinBlock *block, 
                                     ThreeIndOpArray &CDD, int iterCase);
  
  //============================================================================
  //Construct the one-electron part of the V(a) operator
  //============================================================================
  void ConstructDes(SpinBlock &big, Matrix &heff_, Wavefunction &WF, vector<WavefunctionArray> &T,
                    IntegralContainer &IKJA, int *OrbWin);
  
  //============================================================================
  // Construct a+aa|psi> on two blocks block
  //============================================================================
  void ConstructCreDesDesTwoBlocks(SpinBlock &big, Wavefunction &WF, vector<WavefunctionArray> &T, 
                                 IntegralContainer &IKJA,int iterCase, int *OrbWin);
  
  //============================================================================
  //Construct the contributions to the perturber functions V(a) where all indices
  //are on one block
  //============================================================================
  void ConstructVaSingleBlock(SpinBlock &big, Wavefunction &WF, ThreeIndOpArray &CDD,
                              ThreeIndOpArray &CDD_,vector<WavefunctionArray> &Ta, 
                              IntegralContainer &IKJA,NEVPT2Info &Info);
  
  //============================================================================
  //the driver for the construction of the a+aa|psi> operators
  //============================================================================
  void ConstructCreDesDes(SpinBlock &big, vector<Wavefunction> &WF,ThreeIndOpArray &CDD, 
                          IntegralContainer &IKJA, SweepParams &sweepParams, NEVPT2Info &Info);

  
}


#endif	/* RIPDM_OPCONSTRUCT_H */

