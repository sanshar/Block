#ifndef SPIN_DAVIDSON_HEADER
#define SPIN_DAVIDSON_HEADER
#include "wavefunction.h"
#include "spinblock.h"

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
}

#endif
