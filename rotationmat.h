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

This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_ROTATION_MAT_HEADER
#define SPIN_ROTATION_MAT_HEADER 
#include <vector>
#include "Operators.h"
#include "multiarray.h"

namespace SpinAdapted{

void SaveRotationMatrix (const std::vector<int>& sites, const std::vector<Matrix>& m1, int state);
void LoadRotationMatrix (const std::vector<int>& sites, std::vector<Matrix>& m1, int state);
void diagonalise_dm(SparseMatrix& tracedMatrix, SparseMatrix& transformMatrix, std::vector<DiagonalMatrix>& eigenMatrix);
void sort_weights(std::vector<DiagonalMatrix>& eigenMatrix, vector<pair<int, int> >& inorderwts, vector<vector<int> >& weightsbyquanta);
double assign_matrix_by_dm(std::vector<Matrix>& rotatematrix, std::vector<DiagonalMatrix>& eigenmatrix, SparseMatrix& transformmatrix, vector<pair<int, int> >& inorderwts, vector<vector<int> >& wtsbyquanta, int totalstatesbydm, int totalstatesbyquanta, int left_block_size, int right_block_size);
int get_total_states(const int &this_size, const int &other_size);
bool can_connect(int n, int spin, int right_block_size);
}
#endif
