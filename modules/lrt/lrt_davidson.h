/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

// Written by N.N. for DMRG-LRT

#ifndef LRT_SPIN_DAVIDSON_HEADER
#define LRT_SPIN_DAVIDSON_HEADER

#include "davidson.h"
#include <boost/format.hpp>
#include <fstream>
#include <stdio.h>

namespace SpinAdapted {

namespace LRT {

class multiply_h_left : public Davidson_functor
{
private:
  const SpinBlock& block;
public:
  multiply_h_left(const SpinBlock& b, const bool& onedot_);
  void operator() (Wavefunction& c, Wavefunction& v, int state, bool conjugate = false);
  const SpinBlock& get_block() { return block; }
};

class multiply_h_total : public Davidson_functor
{
private:
  const SpinBlock& block;
public:
  multiply_h_total(const SpinBlock& b, const bool& onedot_);
  void operator() (Wavefunction& c, Wavefunction& v, int state, bool conjugate = false);
  const SpinBlock& get_block() { return block; }
};

void LoadDavidsonInfo
(Matrix& a_subspace, Matrix& b_subspace, Matrix& s_subspace, Matrix& d_subspace,
 std::vector<double>& eigenvalues, std::vector<double>& rnorm, std::vector<double>& ynorm,
 int& mroots, int& i_conv_root, bool& deflation_sweep);

void SaveDavidsonInfo
(Matrix& a_subspace, Matrix& b_subspace, Matrix& s_subspace, Matrix& d_subspace,
 std::vector<double>& eigenvalues, std::vector<double>& rnorm, std::vector<double>& ynorm,
 int& mroots, int& i_conv_root, bool& deflation_sweep);


};

};

#endif
