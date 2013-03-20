/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "modules/lrt/lrt_davidson.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "BaseOperator.h"
#include "MatrixBLAS.h"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

SpinAdapted::LRT::multiply_h_left::multiply_h_left(const SpinBlock& b, const bool &onedot_) : block(b){}
SpinAdapted::LRT::multiply_h_total::multiply_h_total(const SpinBlock& b, const bool &onedot_) : block(b){}

void SpinAdapted::LRT::multiply_h_left::operator()(Wavefunction& c, Wavefunction& v, int state, bool conjugate)
{
  block.multiplyH_lrt_left( c, &v, state, conjugate, MAX_THRD);
}

void SpinAdapted::LRT::multiply_h_total::operator()(Wavefunction& c, Wavefunction& v, int state, bool conjugate)
{
  block.multiplyH_lrt_total( c, &v, state, conjugate, MAX_THRD);
}

void SpinAdapted::LRT::LoadDavidsonInfo
(Matrix& a_subspace, Matrix& b_subspace, Matrix& s_subspace, Matrix& d_subspace,
 std::vector<double>& eigenvalues, std::vector<double>& rnorm, std::vector<double>& ynorm,
 int& mroots, int& i_conv_root, bool& deflation_sweep)
{
  std::string file;
  file = str(boost::format("%s%s") % dmrginp.load_prefix() % "/scratch_lrt_davidson.tmp" );
  if(!mpigetrank()) {
    std::ifstream ifs(file.c_str(), std::ios::binary);
    boost::archive::binary_iarchive load_scr(ifs);
    load_scr >> a_subspace
             >> b_subspace
             >> s_subspace
             >> d_subspace
             >> eigenvalues
             >> rnorm
             >> ynorm
             >> mroots
             >> i_conv_root
             >> deflation_sweep;
  }
}

void SpinAdapted::LRT::SaveDavidsonInfo
(Matrix& a_subspace, Matrix& b_subspace, Matrix& s_subspace, Matrix& d_subspace,
 std::vector<double>& eigenvalues, std::vector<double>& rnorm, std::vector<double>& ynorm,
 int& mroots, int& i_conv_root, bool& deflation_sweep)
{
  std::string file;
  file = str(boost::format("%s%s") % dmrginp.load_prefix() % "/scratch_lrt_davidson.tmp" );
  if(!mpigetrank()) {
    std::ofstream ofs(file.c_str(), std::ios::binary);
    boost::archive::binary_oarchive save_scr(ofs);
    save_scr << a_subspace
             << b_subspace
             << s_subspace
             << d_subspace
             << eigenvalues
             << rnorm
             << ynorm
             << mroots
             << i_conv_root
             << deflation_sweep;
  }
}

