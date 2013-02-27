/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_ROTATION_MAT_HEADER
#define SPIN_ROTATION_MAT_HEADER 
#include <vector>
#include "Operators.h"
#include "multiarray.h"
#include "modules/lrt/lrt_rotationmat.h"

namespace SpinAdapted {

namespace LRT {

double assign_matrix_by_dm
(std::vector<DiagonalMatrix>& eigenmatrix,
 std::vector<Matrix>& rotatematrix, std::vector< std::vector<double> >& selectedwts,
 std::vector<Matrix>& rejectedbasis, std::vector< std::vector<double> >& rejectedwts,
 SparseMatrix& transformmatrix, vector<pair<int, int> >& inorderwts, vector<vector<int> >& wtsbyquanta,
 int totalstatesbydm, int totalstatesbyquanta, int left_block_size, int right_block_size);

double SpinAdapted::LRT::assign_matrix_by_dm_deriv
(const std::vector<Matrix>& rotatematrix, const std::vector< std::vector<double> >& selectedwts,
 const std::vector<Matrix>& rejectedbasis, const std::vector< std::vector<double> >& rejectedwts,
 const SparseMatrix& density_deriv, std::vector<Matrix>& rotatematrix_deriv);

void project_onto_rejectedspace
(const Wavefunction& c, const std::vector<Matrix>& rejectedbasis, const bool& dot_with_sys, Wavefunction& c_projected);

};

};

#endif
