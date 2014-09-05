/* 
 * File:   nevpt2_renormalize.h
 * Author: roemelt
 *
 * Created on February 4, 2014, 1:40 PM
 */

#ifndef NEVPT2_RENORMALIZE_H
#define	NEVPT2_RENORMALIZE_H

#include "spinblock.h"
#include "wavefunction.h"
#include "nevpt2_operators.h"
#include "nevpt2_info.h"

namespace SpinAdapted{
  
  //============================================================================
  // Evaluate the rotation matrices for the states that are out of scope of the 
  // regular DMRG algorithm because of their quanta. These states are important
  // subsequent NEVPT2 calculations, in particular for the V(i) class of perturber 
  // functions
  //============================================================================
  void NEVPT2_AddToRotationMat(SpinBlock &big, vector<Matrix> &RotMatrix, SpinBlock &system,
                               vector<Wavefunction> &WF, int sweepIter);
  
  
  //============================================================================
  // Generate a (normalized) density from a set of perturber functions and the
  // reference function
  //============================================================================
  void NEVPT2MakeSpecialDensity(SpinBlock &big,vector <vector<WavefunctionArray> > &T,NEVPT2Info &Info, 
                                DensityMatrix &TotalDensity, vector<Wavefunction> &references);
  
  //============================================================================
  // Add the first order interacting space density to the density of the ground 
  // state (with a weight of 10%)
  //============================================================================
  void AddFOISDensity(SpinBlock &big, DensityMatrix &D, vector<Wavefunction> &WF,
                      NEVPT2Info &Info, int SweepIter, int BlockIter);
  
  
  
}




#endif	/* NEVPT2_RENORMALIZE_H */

