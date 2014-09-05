/* 
 * File:   ripdm_util.h
 * Author: roemelt
 *
 * Created on September 5, 2013, 1:03 PM
 */

#ifndef NEVPT2_UTIL_H
#define	NEVPT2_UTIL_H

#include "nevpt2_info.h"

namespace SpinAdapted{

  //====================================
  //get the Order of the active orbitals
  //====================================
  void GetOrder(vector<int> &reorder, int NActive);
  
  //========================================
  //read the orbital spaces and the BaseName
  //========================================
  void ReadInput(char *BaseName, int *OrbWin, double& ENuc, bool& ConventionalOverlap);
  
  //===============================================
  //read the one-body and two-body density matrices
  //===============================================
  void ReadDensityMatrices(Matrix &D,array_4d<double> &D2, int root, const vector<int> &reorder);
  
  //============================================================================
  // Read the orbital energies from disk
  //============================================================================
  void ReadOrbEnergies(vector<double> &EOrb,const char* BaseName, int *OrbWin, const vector<int> &ReOrder);
  
  //============================================================================
  // Read the previously stored MO integrals
  //    IJKL => (i(1)j(1)|k(2)l(2)) stored as matrices of the form Jij(k,l)
  //    IKJL => (i(1)k(1)|j(2)l(2)) stored as matrices of the form Kij(k,l)
  //    IAJB => (i(1)a(1)|j(2)b(2)) stored as matrices of the form Kij(a,b)
  //    IJKA => (i(1)j(1)|k(2)a(2)) stored as matrices of the form Jij(k,a)
  //    IKJA => (i(1)k(1)|j(2)a(2)) stored as matrices of the form Kij(k,a)
  //============================================================================
  void ReadIntegrals(IntegralContainer &IJKL, IntegralContainer &IKJL, 
                     IntegralContainer &IAJB, IntegralContainer &IJAB, 
                     IntegralContainer &IJKA, IntegralContainer &IKJA, 
                     Matrix &h, int *OrbWin, const char *BaseName,const vector<int> &ReOrder);
  
  //============================================================================
  // Read the prestored MO integrals
  //    IKJL => (i(1)k(1)|j(2)l(2)) stored as matrices of the form Kij(k,l)
  //    IKJA => (i(1)k(1)|j(2)a(2)) stored as matrices of the form Kij(k,a)
  //============================================================================
  void ReadK(IntegralContainer &IKJL,IntegralContainer &IKJA,int *OrbWin, const char *BaseName,
             const vector<int> &ReOrder);
  
  //============================================================================
  // Read the effective one-electron matrices from disk
  //    heff(p,q)  = h(p,q) + 2*(IJ|pq) - (Ip|Jq)
  //    heff_(p,q) = h(p,q) + 2*(IJ|pq) - (Ip|Jq) - (Tp|Tq)
  //============================================================================
  void ReadHeff(int *OrbWin, Matrix &Heff, Matrix &Heff_,const char *BaseName, const vector<int> &ReOrder);
  
  //============================================================================
  // Generate the two effective one-electron matrices
  //    heff(p,q)  = h(p,q) + 2*(IJ|pq) - (Ip|Jq)
  //    heff_(p,q) = h(p,q) + 2*(IJ|pq) - (Ip|Jq) - (Tp|Tq)
  //============================================================================
  void GenerateHeff(int *OrbWin,Matrix &h, Matrix &h_eff, Matrix &h_eff_, 
                    IntegralContainer &IJKL, IntegralContainer &IKJL,
                    IntegralContainer &IJAB, IntegralContainer &IAJB,
                    IntegralContainer &IJKA, IntegralContainer &IKJA);
  
  //============================================================================
  // Build hole densities according to 
  //
  //    D1_(p,q)     = 2*delta(a,b) - D1(q,p)
  //    D2_(p,q,r,s) = D2(r,s,p,q) + delta(p,s)*D1(r,q) + delta(q,r)*D1(s,p)
  //                   -2*delta(p,r)*D1(s,q) -2*delta(q,s)D1(r,p)
  //============================================================================
  void BuildHoleDensities(int *OrbWin,const Matrix &D1, Matrix &D1_, 
                          const array_4d<double> &D2, array_4d<double> &D2_);
  
  //============================================================================
  // Construct the leading term of the 3- and 2-PDM and the three hole density 
  //            DC3(i,j,k,l,m,n) = <psi|E(i,l)E(j,m)E(k,n)|psi>
  //            DC2(i,j,k,l)     = <psi|E(i,k)E(j,l)|psi>
  //============================================================================
  void ConstructAuxPDM(array_6d &DC3, array_6d &D3_, const array_6d &D3, 
                       const array_4d<double> &D2, array_4d<double> &DC2, 
                       const array_4d<double> &D2_, const Matrix &D1, const Matrix &D1_,
                       bool Conventional);
  
  //============================================================================
  // Calculate the energy that arises from taking into account the core orbitals
  // E(core) = h(I,I) + 2*(II|JJ) - (IJ|IJ) + 2(II|TU) D(T,U) - (IT|JU) * D(T,U)
  //============================================================================
  void CalcCoreEnergy(double &CoreEnergy, Matrix &h, Matrix &D, int i0,
                      int i1, int t0, int t1, IntegralContainer &IJKL,
                      IntegralContainer &IKJL, array_4d<double> &D2,double &Eval);

  //============================================================================
  // Calculate the orbital energies from MO integrals according to 
  //    e(p) = h(p,p) + 2(pp|II) - (pI|pI) + D(T,T) {(pp|TT) - 0.5(pT|pT)}
  //============================================================================
  void  CalcOrbEnergies(int*OrbWin, Matrix &h, IntegralContainer &IJKL, IntegralContainer &IKJL, 
                        IntegralContainer &IJAB, IntegralContainer &IAJB, vector<double> &EOrb,
                        Matrix &D); 
  
  //============================================================================
  // a function that evaluates factors for complementary operators in NEVPT2
  //============================================================================
  void EvalCompFactors(vector<double> &Fac, int Gamma, int S_, int S);
  
  //============================================================================
  //a little help function that facilitates the measurement of time
  //============================================================================
  double GetTime();
  
  //=========================
  // print a Matrix to a file
  //=========================
  void PrintMatrix(const Matrix &M, const char *Name);
  void PrintMatrix(const Matrix &M, FILE *f);
  
  //===========================================
  // print a wavefunction or operator to a file
  //===========================================
  void PrintWavefunction(const Baseoperator<Matrix> &WF,const char *FileName);
  void PrintWavefunction(const Baseoperator<Matrix> &WF,const char *FileName, const SpinBlock &big);
  void PrintWavefunction(const Baseoperator<Matrix> &WF,const char *FileName, const StateInfo &Info);
  void PrintWavefunctionProperties(const Baseoperator<Matrix> &WF,const StateInfo *lS,const StateInfo *rS);
  void PrintOperator(const Baseoperator<Matrix> &WF,const char *FileName, const StateInfo &Info);
  
  //=========================
  // print a vector to a file
  //=========================
  void PrintVector(const vector<double> &v, const char* name);
  void PrintVector(const vector<int> &v, const char* name);
  
  //=======================
  // Print a 2pdm to a file
  //=======================
  void Print2PDM(const array_4d<double> &D2, const char *name);
  
  //=======================
  // Print a 3pdm to a file
  //=======================
  void Print3PDM(const array_6d &D3, const char *name);
  
  //============================
  // Print the quanta of a block
  //============================
  void PrintStateInfo(const SpinBlock &block);
  void PrintStateInfo(const StateInfo &Info);
  
  //===============================================
  //a small function that transposes a given Matrix
  //===============================================
  void Transpose(Matrix &M);
  
  //===================================
  // Set all values in a matrix to zero
  //===================================
  void Initialize(Matrix& M);
  void Initialize(array_4d<double> &M);
  
  //========================
  // Read the 3PDM from disk
  //========================
  void Read3PDM(array_6d &D3, const char* FileName);
  
  //===============================================
  // evaluate the factor of complimentary operators
  //===============================================
  double CalcCompfactor(TensorOp& op1, TensorOp& op2, double &iljk);

  //============================================================================
  //give the output for an unconventional NEVPT2 calculation
  //============================================================================
  void GiveNEVPT2Output(NEVPT2Info &Info);
  
  //============================================================================
  //the Kronecker delta
  //============================================================================
  double delta(const int&i, const int &j);
  
  //============================================================================
  // Print A set of rotation matrices
  //============================================================================
  void PrintRotMat(const vector<Matrix> &RotMat, const char *Name);
  void PrintRotMatProperties(const vector<Matrix> &RotMat);
  
  //============================================================================
  // Add a Matrix to a another matrix (taking into account transposition)
  //============================================================================
  void ScaleAdd_(double d, const SparseMatrix& a, SparseMatrix& b);
  
  //============================================================================
  // Transpose a wavefunction
  //============================================================================
  void Transpose(const Wavefunction &oldWave, Wavefunction &NewWave, SpinBlock &big);
  
  
  //============================================================================
  // Evaluate the action of H on a function V(tu), meaning:
  //            sigma(t,u) = H * V(t,u)
  //============================================================================
  void GenerateActionOfH(Wavefunction &WF,vector<boost::shared_ptr<WavefunctionArray> > &Vtu, 
                         vector<boost::shared_ptr <WavefunctionArray> > &Sigmatu,
                         SpinQuantum &TripOpQ, SpinQuantum &SingOpQ, const SpinBlock &big,
                         int *OrbWin);
  
  //============================================================================
  //Save and load the rotation matrices
  //============================================================================
  void NEVPT2SaveRotationMatrix(const std::vector<int>& sites, const std::vector<Matrix>& m1, int state, const char *name);
  void NEVPT2LoadRotationMatrix (const std::vector<int>& sites, std::vector<Matrix>& m1, int state, const char *name);
  
  //============================================================================
  //Determine whether this process should skip generating operators based on 
  //E(dotsite,dotsite)
  //============================================================================
  bool SkipOperator(const vector<int> &sites, const int &dotsite);
  
  //============================================================================
  // Determine the size of a file
  //============================================================================
  int FileSize(const char *FileName);
  
  //============================================================================
  // Check the size of an array
  //============================================================================
  void CheckSize(ThreeIndOpArray &CCD);

  //============================================================================
  // Check if one of a given set of SpinQuanta equals another SpinQuantum
  // return the position of the first element in the vector that equals the
  // reference
  //============================================================================
  int CheckEquality(const SpinQuantum &Reference, const vector<SpinQuantum> &Test);

  //============================================================================
  //establish batching over second index b that runs from start to end-1
  //============================================================================
  void EstablishBatching(vector<pair<int,int> > &batches, int MaxCore, int M, int start,
                     int end, int NOperators);
  
}  




#endif	/* NEVPT2_UTIL_H */

