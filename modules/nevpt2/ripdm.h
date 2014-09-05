/* 
 * File:   ripdm.h
 * Author: roemelt
 *
 * Created on April 3, 2013, 4:35 PM
 */

#ifndef RIPDM_H
#define	RIPDM_H

#include "nevpt2_operators.h"

//==============================================================================
// Module to generate reduced density matrices up to fourth order using  
// configurational RI
//==============================================================================
namespace SpinAdapted {
  
  // Main driver
  void generate_RI_density_matrices(std::vector<Wavefunction>& wavefunctions,
          SpinBlock& big);
  
  // generation and gathering routines for the operators needed
  void generate_terminal_replacement_operators(const SpinBlock &big, const Wavefunction &WF,
          WavefunctionArray &T);
  void gather_terminal_operators(SpinBlock &big, WavefunctionArray &T);
  
  //generation routines for the different density matrices
  void generate_twopdm_RI_a(const SpinBlock &big, const WavefunctionArray &T_rep,
                           array_4d<double>& twopdm, array_4d<double>& twopdm_C, Matrix &onepdm);
  void generate_twopdm_RI_b(const SpinBlock &big, const WavefunctionArray &T_rep,
                          array_4d<double>& twopdm, array_4d<double>& twopdm_spatial);
  void generate_threepdm_RI(SpinBlock &big, WavefunctionArray& T, array_6d& threepdm,
                            array_6d& threepdm_C, const array_4d<double>& twopdm, const Matrix& onepdm, Wavefunction& WF);
  void generate_fourpdm_RI(SpinBlock &big, WavefunctionArray& T,array_4d<double>& fourpdm, Wavefunction &WF);
  void generate_A_RI(SpinBlock &big, WavefunctionArray& T,array_6d& A, array_6d& A_, Wavefunction &WF,array_6d &threepdm);
      
  // the routines that write the density matrices to disk
  void save_onepdm_RI_text(Matrix &onepdm, int root);
  void save_twopdm_RI_text(array_4d<double>& twopdm_spatial, int root, bool predensity=false);
  void save_threepdm_RI_text(array_6d& threepdm, int root, bool predensity=false,int atype=_3PDM_NORMAL_);
  void save_fourpdm_RI_text(SpinBlock& big, array_4d<double>& fourpdm, int root);
  void save_onepdm_RI_binary(Matrix &onepdm, int root);
  void save_twopdm_RI_binary(array_4d<double> &twopdm, int root, bool predensity=false);
  void save_threepdm_RI_binary(array_6d threepdm, int root, bool predensity=false,int atype=_3PDM_NORMAL_);
  void save_fourpdm_RI_binary(array_4d<double> &fourpdm, int root);
  
  // the routines that synchronize or sum up the parallel pdm's
  void SynchronizePDM(array_4d<double> &pdm);
  void SynchronizePDM(array_6d &pdm);
  void SumPDM(array_4d<double> &pdm);
  void SumPDM(array_6d &pdm);
  
}



#endif	/* RIPDM_H */

