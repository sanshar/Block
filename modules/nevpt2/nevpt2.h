/* 
 * File:   ripdm_nevpt2.h
 * Author: roemelt
 *
 * Created on July 26, 2013, 4:07 PM
 */

#ifndef NEVPT2_H
#define	NEVPT2_H

#include "nevpt2_info.h"

namespace SpinAdapted{
  //============================================================================
  // the main driver for BLOCK-NEVPT2 calculations
  // INPUT      big     the spinblock of the whole lattice
  //            WF      the vector wavefunctions for which you want to evaluate
  //                    the NEVPT2 energy
  // also you need to supply an input file with name "dmrg.nevpt2.inp" that 
  // specifies (in that order) the name of your calculation, first and last
  // indices for internal, active and external orbital spaces, the static energy
  // of the nuclei and if the Overlap should be calculated conventionally or not. 
  //============================================================================
  void NEVPT2_Driver(SpinBlock &big, std::vector<Wavefunction> &WF, NEVPT2Info &Info);

  //============================================================================
  // Evaluate the energy contribution from the eight NEVPT2 classes
  //============================================================================
  double V_ijab(int *OrbWin, vector<double> &EOrb, SpinBlock &big, Wavefunction &WF, IntegralContainer &IAJB);
  double V_iab(int *OrbWin, vector<double> &EOrb, SpinBlock &big, Wavefunction &WF, 
               IntegralContainer &IAJB, IntegralContainer &IKJL, Matrix &D1, 
               array_4d<double> &D2, Matrix &Heff);
  double V_ija(int *OrbWin, vector<double> &EOrb, SpinBlock &big, Wavefunction &WF, 
               IntegralContainer &IKJA, IntegralContainer &IKJL, IntegralContainer &IJKL,
               Matrix &D1, Matrix &D1_, array_4d<double> &D2, Matrix &Heff);
  double V_ab(int *OrbWin, const vector<double> &EOrb, const SpinBlock &big,
              Wavefunction &WF, IntegralContainer &IAJB, IntegralContainer &IKJL,
              const array_4d<double> &D2, const array_6d &D3, const Matrix &heff, 
              bool Conventional, char *BaseName, double E0, bool ConventionalOverlap,
              int MaxCore);
  double V_ij(int *OrbWin, const vector<double> &EOrb, const SpinBlock &big,
              Wavefunction &WF, IntegralContainer &IKJL, const array_4d<double> &D2_, 
              const array_6d &D3_, const Matrix &heff_, bool Conventional, double E0, 
              bool ConventionalOverlap, int MaxCore);
    double V_ia(int *OrbWin, const vector<double> &EOrb, Wavefunction &WF, const SpinBlock &big,
              IntegralContainer &IKJA, IntegralContainer &IKJL, IntegralContainer &IJKA,
              boost::shared_ptr<IntegralContainer> KIAJ, const Matrix &D1, 
              const array_4d<double> &DC2, const array_6d &DC3, const Matrix &heff,
              const Matrix &heff_, bool Conventional,char *BaseName, double E0, 
              bool ConventionalOverlap,int MaxCore);
  double V_a(int *OrbWin, const vector<double> &EOrb, const SpinBlock &big,
             Wavefunction &WF, IntegralContainer &IKJL, IntegralContainer & IKJA,
             boost::shared_ptr<IntegralContainer> KIAJ, const array_6d &DC3, 
             const array_4d<double> &DC2, const Matrix &D1, const array_6d &AuxA, 
             const Matrix &heff_, bool Conventional, bool ConventionalOverlap);
  double V_a(vector<WavefunctionArray> &T, SpinBlock &big, Wavefunction &WF, NEVPT2Info &Info, int iroot);
  double V_i(int *OrbWin, const vector<double> &EOrb, const SpinBlock &big,
             Wavefunction &WF, IntegralContainer &IKJL, 
             const array_6d &DC3, const array_4d<double> &DC2, const Matrix &D1, 
             const Matrix &D1_, array_6d &AuxA_, const Matrix &heff_, 
             const Matrix &heff, bool Conventional, bool ConventionalOverlap);
  double V_i(vector<WavefunctionArray> &T, SpinBlock &big, Wavefunction &WF, NEVPT2Info &Info, int iroot);
}



#endif	/* NEVPT2_H */

