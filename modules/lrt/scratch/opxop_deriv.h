#ifndef SPIN_OPXOP_HEADER
#define SPIN_OPXOP_HEADER
 
#include "spinblock.h"
#include "module/dmrg_lrt/spinblock_derived.h"
#include "Operators.h"
#include "module/dmrg_lrt/Operators_derived.h"

namespace SpinAdapted{
namespace opxop
{
  void cxcddcomp_derived(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, SparseMatrix* o);
  void dxccdcomp_derived(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, SparseMatrix* o);
  void cdxcdcomp_derived(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, SparseMatrix* o);
  void dcxdccomp_derived(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, SparseMatrix* o);
  void ddxcccomp_derived(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, SparseMatrix* o);
  void ccxddcomp_derived(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, SparseMatrix* o);
  
  void cxcdcomp_derived(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, int I, SparseMatrix* o, double factor);
  void dxdccomp_derived(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, int I, SparseMatrix* o, double factor);
  void dxcccomp_derived(const SpinBlock* otherBlock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, int I, SparseMatrix* o, double scale);
  void cxddcomp_derived(const SpinBlock* otherBlock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, int I, SparseMatrix* o, double scale);

  void cxcddcomp_derived(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q);
  void dxccdcomp_derived(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q);
  void cdxcdcomp_derived(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q);
  void dcxdccomp_derived(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q);
  void ddxcccomp_derived(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q);
  void ccxddcomp_derived(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q);
}
}
#endif 
