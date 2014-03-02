#ifndef SPIN_OPXOP_HEADER
#define SPIN_OPXOP_HEADER
 
#include "spinblock.h"
#include "Operators.h"

namespace SpinAdapted{
namespace opxop
{
  void cxcddcomp(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, SparseMatrix* o);
  void cdxcdcomp(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, SparseMatrix* o);
  void ddxcccomp(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, SparseMatrix* o);
  
  void cxcddcomp_d(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, DiagonalMatrix* e);
  void cdxcdcomp_d(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, DiagonalMatrix* e);
  void ddxcccomp_d(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, DiagonalMatrix* e);
  

  void cxcdcomp(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, int I, SparseMatrix* o, double factor);
  void dxcccomp(const SpinBlock* otherBlock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, int I, SparseMatrix* o, double scale);

  void cxcddcomp(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q);
  void cdxcdcomp(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q);
  void ddxcccomp(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q); 
}
}
#endif 
