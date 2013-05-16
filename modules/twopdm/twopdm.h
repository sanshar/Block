/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef TWOPDM_HEADER_H
#define TWOPDM_HEADER_H

#include <vector>
#include <multiarray.h>
#include "spinblock.h"
#include "wavefunction.h"
#include "BaseOperator.h"
#include "npdm_patterns.h"

namespace SpinAdapted{

//===========================================================================================================================================================

void compute_twopdm_sweep(std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int state, int sweepPos, int endPos);

void FormLeftOp(const SpinBlock* leftBlock, const SparseMatrix& leftOp, const SparseMatrix& dotOp, SparseMatrix& Aop, int totalspin);
double DotProduct(const Wavefunction& w1, const Wavefunction& w2, double Sz, const SpinBlock& big);
double spinExpectation(Wavefunction& wave1, Wavefunction& wave2, SparseMatrix &leftOp, SparseMatrix& dotOp, SparseMatrix& rightOp, const SpinBlock& big);

void loop_over_block_operators( Wavefunction & wavefunction, 
                                const SpinBlock & big, 
                                std::vector<Npdm::CD> & lhs_cd,
                                std::vector<Npdm::CD> & dot_cd,
                                std::vector<Npdm::CD> & rhs_cd,
                                array_4d<double> & twopdm );

void assign_antisymmetric(array_4d<double>& twopdm, const int i, const int j, const int k, const int l, const double val);
void save_twopdm_text(const array_4d<double>& twopdm, const int &i, const int &j);
void save_spatial_twopdm_text(const array_4d<double>& twopdm, const int &i, const int &j);
void save_spatial_twopdm_binary(const array_4d<double>& twopdm, const int &i, const int &j);
void save_twopdm_binary(const array_4d<double>& twopdm, const int &i, const int &j);
void load_twopdm_binary(array_4d<double>& twopdm, const int &i, const int &j);
void save_averaged_twopdm(const int &nroots);
void accumulate_twopdm(array_4d<double>& twopdm);

std::vector<int> distribute_procs(const int numprocs, const int numjobs);
std::vector<int> distribute_procs(const int numprocs, const std::vector<int>& sites);

//===========================================================================================================================================================

}

#endif

