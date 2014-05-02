/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef TWOPDM_HEADER_H
#define TWOPDM_HEADER_H
#include "spinblock.h"
#include "wavefunction.h"
#include "BaseOperator.h"
#include <vector>
#include <multiarray.h>

namespace SpinAdapted{
enum Oporder {CD_CD, CC_DD, CDt_CD, C_CD_D, D_CD_C, CC_D_D, D_CC_D, CD_D_C};

void assign_antisymmetric(array_4d<double>& twopdm, const int i, const int j, const int k, const int l, const double val); // done

 void compute_twopdm_sweep(std::vector<Wavefunction>& solutions, const SpinBlock& system, const SpinBlock& systemDot, const SpinBlock& newSystem, const SpinBlock& newEnvironment, const SpinBlock& big, const int numprocs, int state); // done
 void compute_twopdm_initial(std::vector<Wavefunction>& solutions, const SpinBlock& system, const SpinBlock& systemDot, const SpinBlock& newSystem, const SpinBlock& newEnvironment, const SpinBlock& big, const int numprocs, int state); // done
 void compute_twopdm_final(std::vector<Wavefunction>& solutions, const SpinBlock& system, const SpinBlock& systemDot, const SpinBlock& newSystem, const SpinBlock& newEnvironment, const SpinBlock& big, const int numprocs, int state); // done

void compute_two_pdm_0_2_2(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm); // done
void compute_two_pdm_2_0_2(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm); // done
void compute_two_pdm_2_2_0(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm); // done

void compute_two_pdm_1_2_1(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm); // done
void compute_two_pdm_1_1_2(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm); // done
void compute_two_pdm_2_1_1(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm); // done

void compute_two_pdm_0_4_0(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm); // done
void compute_two_pdm_0_0_4(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm); // done
void compute_two_pdm_4_0_0(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm); // done

void compute_two_pdm_0_3_1_notranspose(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm);
void compute_two_pdm_0_3_1(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm); // done
void compute_two_pdm_1_3_0(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm); // done
void compute_two_pdm_3_1_0(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm); // done
void compute_two_pdm_3_0_1(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm); // done
void compute_two_pdm_1_3(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm);//0_1_3 and 1_0_3 // done

void spinExpectation(Wavefunction& wave1, Wavefunction& wave2, SparseMatrix &leftOp, SparseMatrix& dotOp, SparseMatrix& rightOp, const SpinBlock& big, vector<double>& expectations, bool doTranspose);  // Done
void FormLeftOp(const SpinBlock* leftBlock, const SparseMatrix& leftOp, const SparseMatrix& dotOp, SparseMatrix& Aop, int totalspin); // Done
void spin_to_nonspin(vector<int>& indices, vector<double>& coeffs, array_4d<double>& twopdm, Oporder order, bool dotranspose);  // Done

void save_twopdm_text(const array_4d<double>& twopdm, const int &i, const int &j);
void save_spatial_twopdm_text(const array_4d<double>& twopdm, const int &i, const int &j);
void save_spatial_twopdm_binary(const array_4d<double>& twopdm, const int &i, const int &j);
void save_twopdm_binary(const array_4d<double>& twopdm, const int &i, const int &j);
void load_twopdm_binary(array_4d<double>& twopdm, const int &i, const int &j);
void save_averaged_twopdm(const int &nroots);
void accumulate_twopdm(array_4d<double>& twopdm);
double DotProduct(const Wavefunction& w1, const Wavefunction& w2, double Sz, const SpinBlock& big);
std::vector<int> distribute_procs(const int numprocs, const int numjobs);
std::vector<int> distribute_procs(const int numprocs, const std::vector<int>& sites);
}
#endif
