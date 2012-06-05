/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is free software: you can redistribute it and/or modify         
it under the terms of the GNU General Public License as published by         
the Free Software Foundation, either version 3 of the License, or            
(at your option) any later version.                                          
                                                                             
This program is distributed in the hope that it will be useful,              
but WITHOUT ANY WARRANTY; without even the implied warranty of               
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
GNU General Public License for more details.                                 
                                                                             
You should have received a copy of the GNU General Public License            
along with this program.  If not, see <http://www.gnu.org/licenses/>.        
*/


#ifndef ONEPDM_HEADER_H
#define ONEPDM_HEADER_H
#include "spinblock.h"
#include "wavefunction.h"
#include "BaseOperator.h"
#include <vector>
#include <newmat.h>

namespace SpinAdapted{
void compute_onepdm(std::vector<Wavefunction>& solutions, const SpinBlock& system, const SpinBlock& systemDot, const SpinBlock& newSystem, const SpinBlock& newEnvironment, const SpinBlock& big, const int numprocs);

void compute_one_pdm_0_2(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, Matrix& onepdm);
void compute_one_pdm_2_0(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, Matrix& onepdm);
void compute_one_pdm_1_1(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, Matrix& onepdm);

void save_onepdm_spatial_text(const Matrix& onepdm, const int &i, const int &j);
void save_onepdm_spatial_binary(const Matrix& onepdm, const int &i, const int &j);
void save_onepdm_binary(const Matrix& onepdm, const int &i, const int &j);
void load_onepdm_binary(Matrix& onepdm, const int &i, const int &j);
void save_onepdm_text(const Matrix& onepdm, const int &i, const int &j);
void accumulate_onepdm(Matrix& onepdm);
std::vector<int> distribute_procs(const int numprocs, const int numjobs);
}
#endif
