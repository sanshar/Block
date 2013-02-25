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

namespace SpinAdapted{

struct Davidson_functor
{
  virtual void operator()(Wavefunction& c, Wavefunction& v) = 0;
  virtual const SpinBlock& get_block() = 0;
};

class multiply_h : public Davidson_functor
{
 private:
  const SpinBlock& block;
 public:
  multiply_h(const SpinBlock& b, const bool &onedot_);
  void operator()(Wavefunction& c, Wavefunction& v);
  const SpinBlock& get_block() {return block;}
};

};

#endif
