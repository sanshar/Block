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
  void operator() (Wavefunction& c, Wavefunction& v, int iState, int jState);
  const SpinBlock& get_block() { return block; }
};

class multiply_h_total : public Davidson_functor
{
private:
  const SpinBlock& block;
public:
  multiply_h_total(const SpinBlock& b, const bool& onedot_);
  void operator() (Wavefunction& c, Wavefunction& v, int iState, int jState);
  const SpinBlock& get_block() { return block; }
};

void LoadDavidsonInfo(Matrix& h_subspace, Matrix& s_subspace, int& mroots, int& i_conv_root, bool& deflation_sweep)
{
  std::string file;
  file = str(boost::format("%s%s") % dmrginp.load_prefix() % "/scratch_lrt_davidson.tmp" );

  if(!mpigetrank()) {
    std::ifstream ifs(file.c_str(), std::ios::binary);
    boost::archive::binary_iarchive load_scr(ifs);
    load_scr >> h_subspace >> s_subspace >> mroots >> i_conv_root >> deflation_sweep;
  }
}

void SaveDavidsonInfo(Matrix& h_subspace, Matrix& s_subspace, int& mroots, int& i_conv_root, bool& deflation_sweep)
{
  std::string file;
  file = str(boost::format("%s%s") % dmrginp.load_prefix() % "/scratch_lrt_davidson.tmp" );

  if(!mpigetrank()) {
    std::ifstream ofs(file.c_str(), std::ios::binary);
    boost::archive::binary_oarchive save_scr(ofs);
    save_scr << h_subspace << s_subspace << mroots << i_conv_root << deflation_sweep;
  }
}

};

};

#endif
