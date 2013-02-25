#ifndef SPIN_OPXOP_GENERIC_HEADER
#define SPIN_OPXOP_GENERIC_HEADER
 
#include "spinblock.h"
#include "Operators.h"

// generalized opxop, e.g., < LR(I) | A B | LR(J) > := < L(I) | A | L(J) > * < R(I) | B | R(J) >

namespace SpinAdapted {

namespace opxop {

// for solving correction equation
void cxcddcomp(const SpinBlock* otherblock,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1_trans,
               const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q, int iState = 0, int jState = 0);

void cdxcdcomp(const SpinBlock* otherblock,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1_trans,
               const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q, int iState = 0, int jState = 0);

void ddxcccomp(const SpinBlock* otherblock,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1_trans,
               const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q, int iState = 0, int jState = 0);

// for blocking
void cxcddcomp(const SpinBlock* otherblock,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1_trans,
               const SpinBlock* b, SparseMatrix* o, int iState = 0, int jState = 0);

void cdxcdcomp(const SpinBlock* otherblock,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1_trans,
               const SpinBlock* b, SparseMatrix* o, int iState = 0, int jState = 0);

void ddxcccomp(const SpinBlock* otherblock,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1_trans,
               const SpinBlock* b, SparseMatrix* o, int iState = 0, int jState = 0);
  
void cxcdcomp (const SpinBlock* otherblock,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1_trans,
               const SpinBlock* b, int I, SparseMatrix* o, double scale, int iState = 0, int jState = 0);

void dxcccomp (const SpinBlock* otherBlock,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1,
                     std::vector< boost::shared_ptr<SparseMatrix> >& op1_trans,
               const SpinBlock* b, int I, SparseMatrix* o, double scale, int iState = 0, int jState = 0);

};

};

#endif 
