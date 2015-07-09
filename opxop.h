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

  void dxcdcomp(const SpinBlock* otherBlock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, int I, SparseMatrix* o, double scale);
  void cxddcomp(const SpinBlock* otherBlock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, int K, SparseMatrix* o, double scale);

  void cxcddcomp(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q);
  void cdxcdcomp(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q);
  void ddxcccomp(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q); 

  
  //these are only used when left and right states are different
  void dcxdccomp(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, SparseMatrix* o);
  void dxccdcomp(const SpinBlock* CDDblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, SparseMatrix* o);
  void dxccdcomp(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q);
  void dcxdccomp(const SpinBlock* otherblock, std::vector< boost::shared_ptr<SparseMatrix> >& op1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q);

void cdd_cxddcomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, SparseMatrix* o);

void cdd_dxcdcomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, SparseMatrix* o);

void cdd_cxddcomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q);

void cdd_dxcdcomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q);

void ccd_dxcccomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, SparseMatrix* o);

void ccd_cxcdcomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, SparseMatrix* o);

void ccd_dxcccomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q);

void ccd_cxcdcomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q);

}
}
#endif 
