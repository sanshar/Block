#ifndef SPIN_OPXOP_GENERIC_HEADER
#define SPIN_OPXOP_GENERIC_HEADER
 
#include "spinblock.h"
#include "Operators.h"

// generalized opxop, e.g., < LR(I) | A B | LR(J) > := < L(I) | A | L(J) > * < R(I) | B | R(J) >

namespace SpinAdapted {

namespace opxop {

namespace generic {

// for solving correction equation
void cxcddcomp(const SpinBlock* otherblock,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1_trans,
               const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q, opTypes state_index_ij);

void cdxcdcomp(const SpinBlock* otherblock,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1_trans,
               const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q, opTypes state_index_ij);

void ddxcccomp(const SpinBlock* otherblock,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1_trans,
               const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q, opTypes state_index_ij);

// for blocking
void cxcddcomp(const SpinBlock* otherblock,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1_trans,
               const SpinBlock* b, SparseMatrix* o, opTypes state_index_ij);

void cdxcdcomp(const SpinBlock* otherblock,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1_trans,
               const SpinBlock* b, SparseMatrix* o, opTypes state_index_ij);

void ddxcccomp(const SpinBlock* otherblock,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1_trans,
               const SpinBlock* b, SparseMatrix* o, opTypes state_index_ij);
  
void cxcdcomp (const SpinBlock* otherblock,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1,
               const SpinBlock* b, int I, SparseMatrix* o, double scale, opTypes state_index_ij);

void dxcccomp (const SpinBlock* otherBlock,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1_trans,
               const SpinBlock* b, int I, SparseMatrix* o, double scale, opTypes state_index_ij);

};

};

};

#endif 
