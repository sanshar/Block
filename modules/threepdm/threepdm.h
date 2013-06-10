/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef THREEPDM_HEADER_H
#define THREEPDM_HEADER_H

#include <vector>
#include <multiarray.h>
#include "spinblock.h"
#include "wavefunction.h"
#include "BaseOperator.h"
#include "npdm_patterns.h"

namespace SpinAdapted{

//===========================================================================================================================================================

void compute_threepdm_sweep(std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int state, int sweepPos, int endPos);

void threepdm_loop_over_block_operators( Wavefunction & wavefunction, 
                                         const SpinBlock & big, 
                                         std::vector<Npdm::CD> & lhs_cd,
                                         std::vector<Npdm::CD> & dot_cd,
                                         std::vector<Npdm::CD> & rhs_cd,
                                         array_6d<double> & threepdm );

void assign_threepdm_antisymmetric(array_6d<double>& threepdm, const int i, const int j, const int k, const int l, const double val);
void save_threepdm_text(const array_6d<double>& threepdm, const int &i, const int &j);
void save_spatial_threepdm_text(const array_6d<double>& threepdm, const int &i, const int &j);
void save_spatial_threepdm_binary(const array_6d<double>& threepdm, const int &i, const int &j);
void save_threepdm_binary(const array_6d<double>& threepdm, const int &i, const int &j);
void load_threepdm_binary(array_6d<double>& threepdm, const int &i, const int &j);
void save_averaged_threepdm(const int &nroots);
void accumulate_threepdm(array_6d<double>& threepdm);

std::vector<int> distribute_procs(const int numprocs, const int numjobs);
std::vector<int> distribute_procs(const int numprocs, const std::vector<int>& sites);

//===========================================================================================================================================================

}

#endif

