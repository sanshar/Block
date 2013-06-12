/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef FOURPDM_HEADER_H
#define FOURPDM_HEADER_H

#include <vector>
#include <multiarray.h>
#include "spinblock.h"
#include "wavefunction.h"
#include "BaseOperator.h"
#include "npdm_patterns.h"

namespace SpinAdapted{

//===========================================================================================================================================================

void compute_fourpdm_sweep(std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int state, int sweepPos, int endPos);

void fourpdm_loop_over_block_operators( Wavefunction & wavefunction, 
                                         const SpinBlock & big, 
                                         std::vector<Npdm::CD> & lhs_cd,
                                         std::vector<Npdm::CD> & dot_cd,
                                         std::vector<Npdm::CD> & rhs_cd,
                                         array_8d<double> & fourpdm );

void assign_fourpdm_antisymmetric(array_8d<double>& fourpdm, const int i, const int j, const int k, const int l, const double val);
void save_fourpdm_text(const array_8d<double>& fourpdm, const int &i, const int &j);
void save_spatial_fourpdm_text(const array_8d<double>& fourpdm, const int &i, const int &j);
void save_spatial_fourpdm_binary(const array_8d<double>& fourpdm, const int &i, const int &j);
void save_fourpdm_binary(const array_8d<double>& fourpdm, const int &i, const int &j);
void load_fourpdm_binary(array_8d<double>& fourpdm, const int &i, const int &j);
void save_averaged_fourpdm(const int &nroots);
void accumulate_fourpdm(array_8d<double>& fourpdm);

std::vector<int> distribute_procs(const int numprocs, const int numjobs);
std::vector<int> distribute_procs(const int numprocs, const std::vector<int>& sites);

//===========================================================================================================================================================

}

#endif

