#include "opxop.h"
#include "operatorfunctions.h"
#include "wavefunction.h"
#include "BaseOperator.h"
#include "tensor_operator.h"
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0 
#endif
#include "guess_wavefunction.h"
#include "pario.h"

//using namespace operatorfunctions;


/********************************************
Formulas for making hamiltonian matrix while blocking a block with a dot block
********************************************/


void SpinAdapted::opxop::cdxcdcomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, SparseMatrix* o) {
  int ilock = omp_get_thread_num();
  int numthrds = 1;//MAX_THRD;
  const SpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();    

  for (int opind=0; opind<opvec1.size(); opind++) { // this is CreDes_{ij}
    boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind)->getworkingrepresentation(loopblock);
    int i = op1->get_orbs(0);
    int j = op1->get_orbs(1);
    if (!otherblock->get_op_array(CRE_DESCOMP).has_local_index(i,j))
      return;
    boost::shared_ptr<SparseMatrix> op2 = otherblock->get_op_array(CRE_DESCOMP).get_element(i, j).at(opind)->getworkingrepresentation(otherblock);
    double factor = 1.0;
    SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], factor, numthrds);
    if (i != j) {
      //If we have all the operators we dont have to take transposes, useful for <bra|H|ket> evaluation
      if (otherblock->has(DES_CRECOMP)) {
	    op1 = loopblock->get_op_array(DES_CRE).get_element(i,j).at(opind)->getworkingrepresentation(loopblock);
	    double parity = 1.0;
	    if (dmrginp.spinAdapted() == true && dmrginp.hamiltonian() != BCS)
	      parity = getCommuteParity(-getSpinQuantum(i), getSpinQuantum(j), op1->get_deltaQuantum()[0]);
	    op2 = otherblock->get_op_array(DES_CRECOMP).get_element(i,j).at(opind)->getworkingrepresentation(otherblock);
	    SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], factor*parity, numthrds);
      } else {
	    SpinAdapted::operatorfunctions::TensorProduct(otherblock, Transposeview(*op2), Transposeview(*op1), b, &(b->get_stateInfo()), o[ilock], factor, numthrds);
      }
    }
  }
}

void SpinAdapted::opxop::ddxcccomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, SparseMatrix* o)
{
  int ilock = omp_get_thread_num();
  int numthrds = 1;//MAX_THRD;
  const SpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
  
  for (int opind=0; opind<opvec1.size(); opind++) {
    boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind)->getworkingrepresentation(loopblock); // CC_{ij}
    int i = op1->get_orbs(0);
    int j = op1->get_orbs(1);
    if (!otherblock->get_op_array(DES_DESCOMP).has_local_index(i,j))
      return;
    double factor = 2.0; if (i==j) factor = 1.0;
    boost::shared_ptr<SparseMatrix> op2 = otherblock->get_op_array(DES_DESCOMP).get_element(i, j).at(opind)->getworkingrepresentation(otherblock);

    double scale = 1.0;
    double parity = 1.0;
    if (otherblock == b->get_leftBlock())
      parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0));

    SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], parity*factor, numthrds);

    //If we have all the operators we dont have to take transposes, useful for <bra|H|ket> evaluation
    if (otherblock->has(CRE_CRECOMP)) {
      op2 = otherblock->get_op_array(CRE_CRECOMP).get_element(i, j).at(opind)->getworkingrepresentation(otherblock);
      op1 = loopblock->get_op_array(DES_DES).get_element(i, j).at(opind)->getworkingrepresentation(loopblock);
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], parity*factor, numthrds);
    } else {
      Transposeview top1 = Transposeview(*op1);
      Transposeview top2 = Transposeview(*op2);
      
      SpinQuantum sq1 = op1->get_deltaQuantum(0);
      SpinQuantum sq2 = op2->get_deltaQuantum(0);
      double parity2 =TensorOp::getTransposeFactorDD(i, j, sq1.get_s().getirrep(), sq1.get_symm().getirrep());
      parity2*=TensorOp::getTransposeFactorDD(i, j, sq2.get_s().getirrep(), sq2.get_symm().getirrep());
      
      parity *= parity2;
      
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, top2, top1, b, &(b->get_stateInfo()), o[ilock], parity*factor, numthrds);
    }
  }
}



void SpinAdapted::opxop::cxcddcomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, SparseMatrix* o)
{
  int ilock = omp_get_thread_num();
  int numthrds = 1;//MAX_THRD;
  const SpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  for (int opind=0; opind<opvec1.size(); opind++) {
    boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind)->getworkingrepresentation(loopblock); // CRE_i
    int i = op1->get_orbs(0);
    if (!otherblock->get_op_array(CRE_CRE_DESCOMP).has_local_index(i))
      return;

    //If we have all the operators we dont have to take transposes, useful for <bra|H|ket> evaluation
    if (loopblock->has(DES)) {
      boost::shared_ptr<SparseMatrix> op2 = otherblock->get_op_array(CRE_DES_DESCOMP).get_element(i).at(opind)->getworkingrepresentation(otherblock);
      double scale = 1.0;
      double parity = 1.0;
      if (otherblock == b->get_leftBlock()) parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0));
      else parity = 1.0;
      
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], scale*parity, numthrds);	    

      if (otherblock == b->get_rightBlock()) parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0));
      else parity = 1.0;

      op1 = loopblock->get_op_array(DES).get_element(i).at(opind)->getworkingrepresentation(loopblock);
      op2 = otherblock->get_op_array(CRE_CRE_DESCOMP).get_element(i).at(opind)->getworkingrepresentation(otherblock);      
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], scale*parity, numthrds);      
    } else {
      Transposeview top1 = Transposeview(op1);  // DES_i
      boost::shared_ptr<SparseMatrix> op2 = otherblock->get_op_array(CRE_CRE_DESCOMP).get_element(i).at(opind)->getworkingrepresentation(otherblock); // CCD_i
      
      double scale = 1.0;
      double parity = 1.0;
      
      if (otherblock == b->get_rightBlock()) parity = getCommuteParity(-op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0));
      
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, top1, b, &(b->get_stateInfo()), o[ilock], scale*parity, numthrds);
      
      // complex conjugate
      Transposeview top2 = Transposeview(op2); // CDD_i
      if (otherblock == b->get_leftBlock()) parity = getCommuteParity(-op2->get_deltaQuantum(0), op1->get_deltaQuantum(0), o->get_deltaQuantum(0));
      else parity = 1.0;
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, top2, *op1, b, &(b->get_stateInfo()), o[ilock], scale*parity, numthrds);	    
    }
  }
}

//***************************************************************************************





/********************************************
Formulas for multiplying hamiltonian with wavefunction without ever making the hamiltonian explicitly
********************************************/

void SpinAdapted::opxop::cdxcdcomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q)
{
  int ilock = omp_get_thread_num();
  int numthrds = 1;//MAX_THRD;
  SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));
  const SpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
    
  for (int opind=0; opind<opvec1.size(); opind++) {
    boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind)->getworkingrepresentation(loopblock);
    int i = op1->get_orbs(0);
    int j = op1->get_orbs(1);
    if (!otherblock->get_op_array(CRE_DESCOMP).has_local_index(i,j))
      return;
    boost::shared_ptr<SparseMatrix> op2 = otherblock->get_op_array(CRE_DESCOMP).get_element(i, j).at(opind)->getworkingrepresentation(otherblock);
    double factor = 1.0;
    SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v[ilock], hq, factor);
    if (i != j) {
      //If we have all the operators we dont have to take transposes, useful for <bra|H|ket> evaluation
      if (otherblock->has(DES_CRECOMP)) {
	    op1 = loopblock->get_op_array(DES_CRE).get_element(i,j).at(opind)->getworkingrepresentation(loopblock);
	    double parity = 1.0;
	    if (dmrginp.spinAdapted() == true && dmrginp.hamiltonian() != BCS)
	      parity = getCommuteParity(-getSpinQuantum(i), getSpinQuantum(j), op1->get_deltaQuantum()[0]);
	    op2 = otherblock->get_op_array(DES_CRECOMP).get_element(i,j).at(opind)->getworkingrepresentation(otherblock);
        SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v[ilock], hq, parity*factor);
      }
      else 
	    SpinAdapted::operatorfunctions::TensorMultiply(otherblock, Transposeview(*op2), Transposeview(*op1), b, c, v[ilock], hq, factor);
    }
  }
}

void SpinAdapted::opxop::ddxcccomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q)
{
  int ilock = omp_get_thread_num();
  int numthrds = 1;//MAX_THRD;
  SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));
  const SpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
  
  for (int opind=0; opind<opvec1.size(); opind++) {
    boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind)->getworkingrepresentation(loopblock);
    int i = op1->get_orbs(0);
    int j = op1->get_orbs(1);
    if (!otherblock->get_op_array(DES_DESCOMP).has_local_index(i,j))
      return;
    double factor = 2.0; if (i==j) factor = 1.0;
    boost::shared_ptr<SparseMatrix> op2 = otherblock->get_op_array(DES_DESCOMP).get_element(i, j).at(opind)->getworkingrepresentation(otherblock);

    double scale = 1.0;
    double parity = 1.0;
    if (otherblock == b->get_leftBlock())
      parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), hq);

    SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v[ilock], hq, factor*parity);

    //If we have all the operators we dont have to take transposes, useful for <bra|H|ket> evaluation
    if (otherblock->has(CRE_CRECOMP)) {
      op2 = otherblock->get_op_array(CRE_CRECOMP).get_element(i, j).at(opind)->getworkingrepresentation(otherblock);
      op1 = loopblock->get_op_array(DES_DES).get_element(i, j).at(opind)->getworkingrepresentation(loopblock);
      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v[ilock], hq, factor*parity);
    } else {
      Transposeview top1 = Transposeview(*op1);
      Transposeview top2 = Transposeview(*op2);
      
      SpinQuantum sq1 = op1->get_deltaQuantum(0);
      SpinQuantum sq2 = op2->get_deltaQuantum(0);
      double parity2 =TensorOp::getTransposeFactorDD(i, j, sq1.get_s().getirrep(), sq1.get_symm().getirrep());
      parity2*=TensorOp::getTransposeFactorDD(i, j, sq2.get_s().getirrep(), sq2.get_symm().getirrep());
      
      parity *= parity2;
      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, top2, top1, b, c, v[ilock], hq, factor*parity);
    }
  }
}



void SpinAdapted::opxop::cxcddcomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q)
{
  int ilock = omp_get_thread_num();
  int numthrds = 1;//MAX_THRD;
  SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));  // in get_parity, number part is not used
  const SpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  for (int opind=0; opind<opvec1.size(); opind++) {
    boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind);
    int i = op1->get_orbs(0);
    if (!otherblock->get_op_array(CRE_CRE_DESCOMP).has_local_index(i))
      return;

    //If we have all the operators we dont have to take transposes, useful for <bra|H|ket> evaluation
    if (loopblock->has(DES) ) {
      double scale = 1.0;
      double parity = 1.0;
      {
	boost::shared_ptr<SparseMatrix> op2 = otherblock->get_op_array(CRE_DES_DESCOMP).get_element(i).at(opind)->getworkingrepresentation(otherblock);
	boost::shared_ptr<SparseMatrix> op1rep = op1->getworkingrepresentation(loopblock);
	
	if (otherblock == b->get_leftBlock())
	  parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), hq);
	
	SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1rep, b, c, v[ilock], hq, scale*parity);	    
      }
      {
	boost::shared_ptr<SparseMatrix> op1 = loopblock->get_op_array(DES).get_element(i).at(opind)->getworkingrepresentation(loopblock);
	boost::shared_ptr<SparseMatrix> op2 = otherblock->get_op_array(CRE_CRE_DESCOMP).get_element(i).at(opind)->getworkingrepresentation(otherblock);

	if (otherblock == b->get_rightBlock()) parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), hq);
	else parity = 1.0;
	
	SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v[ilock], hq, scale*parity);
      }
    } else {
      boost::shared_ptr<SparseMatrix> op1rep = op1->getworkingrepresentation(loopblock);
      Transposeview top1 = Transposeview(*op1rep);
      boost::shared_ptr<SparseMatrix> op2 = otherblock->get_op_array(CRE_CRE_DESCOMP).get_element(i).at(opind)->getworkingrepresentation(otherblock);
      
      double scale = 1.0;
      double parity = 1.0;
      if (otherblock == b->get_rightBlock())
	parity = getCommuteParity(-op1rep->get_deltaQuantum(0), op2->get_deltaQuantum(0), hq);
      
      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, top1, b, c, v[ilock], hq, scale*parity);
      Transposeview top2 = Transposeview(*op2);    
      if (otherblock == b->get_leftBlock()) parity = getCommuteParity(op1rep->get_deltaQuantum(0), -op2->get_deltaQuantum(0), hq);
      else parity = 1.0;
      
      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, top2, *op1rep, b, c, v[ilock], hq, scale*parity);	    
    }
  }
}


//***************************************************************************************************

/********************************************
Formulas for making diagonal hamiltonian matrix while blocking system and environment blocks
********************************************/


void SpinAdapted::opxop::cdxcdcomp_d(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, DiagonalMatrix* e)
{
  int ilock = omp_get_thread_num();
  int numthrds = 1;//MAX_THRD;
  const SpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  for (int opind=0; opind<opvec1.size(); opind++) {
    boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind)->getworkingrepresentation(loopblock);
    int i = op1->get_orbs(0);
    int j = op1->get_orbs(1);
    if (!otherblock->get_op_array(CRE_DESCOMP).has_local_index(i,j))
      return;
    boost::shared_ptr<SparseMatrix> op2 = otherblock->get_op_array(CRE_DESCOMP).get_element(i, j).at(opind)->getworkingrepresentation(otherblock);
    double factor = 1.0;
    SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), e[ilock], factor);
    if (i != j)
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, Transposeview(*op2), Transposeview(*op1), b, &(b->get_stateInfo()), e[ilock], factor);
  }
}

void SpinAdapted::opxop::ddxcccomp_d(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, DiagonalMatrix* e)
{
  int ilock = omp_get_thread_num();
  int numthrds = 1;//MAX_THRD;
  SpinQuantum q(0,SpinSpace(0),IrrepSpace(0));
  const SpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
  
  for (int opind=0; opind<opvec1.size(); opind++) {
    boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind)->getworkingrepresentation(loopblock);
    int i = op1->get_orbs(0);
    int j = op1->get_orbs(1);
    if (!otherblock->get_op_array(DES_DESCOMP).has_local_index(i,j))
      return;
    double factor = 2.0; if (i==j) factor = 1.0;
    boost::shared_ptr<SparseMatrix> op2 = otherblock->get_op_array(DES_DESCOMP).get_element(i, j).at(opind)->getworkingrepresentation(otherblock);

    double scale = 1.0;
    
    SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), e[ilock], scale*factor);
    
    Transposeview top1 = Transposeview(*op1);
    Transposeview top2 = Transposeview(*op2);
    SpinAdapted::operatorfunctions::TensorProduct(otherblock, top2, top1, b, &(b->get_stateInfo()), e[ilock], scale*factor);
  }
}

void SpinAdapted::opxop::cxcddcomp_d(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, DiagonalMatrix* e)
{
  int ilock = omp_get_thread_num();
  int numthrds = 1;//MAX_THRD;
  const SpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  for (int opind=0; opind<opvec1.size(); opind++) {
    boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind)->getworkingrepresentation(loopblock);
    int i = op1->get_orbs(0);
    if (!otherblock->get_op_array(CRE_CRE_DESCOMP).has_local_index(i))
      return;
    Transposeview top1 = Transposeview(op1);
    boost::shared_ptr<SparseMatrix> op2 = otherblock->get_op_array(CRE_CRE_DESCOMP).get_element(i).at(opind)->getworkingrepresentation(otherblock);

    double scale = 1.0;
    double parity = 1.0;
    //if (otherblock == b->get_rightBlock())
      //parity = getCommuteParity(op1->get_deltaQuantum(), op2->get_deltaQuantum(), q);
    SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, top1, b, &(b->get_stateInfo()), e[ilock], scale*parity);

    Transposeview top2 = Transposeview(op2);

    SpinAdapted::operatorfunctions::TensorProduct(otherblock, top2, *op1, b, &(b->get_stateInfo()), e[ilock], scale*parity);	    
  }
}


//************************************************************************************



/********************************************
Formulas for making CCdcomp operators while blocking a block with a dot block
********************************************/

void SpinAdapted::opxop::cxcdcomp(const SpinBlock* otherBlock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, int I, SparseMatrix* o, double scale)
{
  int ilock = 0;//omp_get_thread_num();
  int numthrds = 1;
  const SpinBlock* loopblock = (otherBlock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  if (opvec1[0]->get_orbs(0) >= I) // opvec1 is CRE
  {
    for (int opind=0; opind<opvec1.size(); opind++) {    
      boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind)->getworkingrepresentation(loopblock);
      if (!otherBlock->get_op_array(CRE_DESCOMP).has_local_index(op1->get_orbs(0), I)) // we have c_J d_I *c_K d_L
	    return;

      const std::vector<boost::shared_ptr<SparseMatrix> >& opvec2 = otherBlock->get_op_array(CRE_DESCOMP).get_element(op1->get_orbs(0), I); // CD_comp(j,i) have multiple matrices because of spin adaption
      for (int opind2 = 0; opind2<opvec2.size(); opind2++) {
	    boost::shared_ptr<SparseMatrix> op2 = opvec2.at(opind2)->getworkingrepresentation(otherBlock); // CD
	    vector<SpinQuantum> op2q = op2->get_deltaQuantum(), op1q = op1->get_deltaQuantum(), oq = o->get_deltaQuantum(); // o is the resulted CCD
	    int j2 = op2q[0].get_s().getirrep(), j1 = op1q[0].get_s().getirrep(), j21 = oq[0].get_s().getirrep();
	    int l2 = op2q[0].get_symm().getirrep(), l1 = op1q[0].get_symm().getirrep(), l21 = oq[0].get_symm().getirrep(), l3 = (-SymmetryOfOrb(I)).getirrep();
	    double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
      if (NonabelianSym)
	    factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());

	    double parity = 1.0;
	    if (otherBlock == b->get_rightBlock())
	      parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0)); // doesn't depend on nelec
	    factor*= parity;

	    SpinAdapted::operatorfunctions::TensorProduct(otherBlock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], factor*scale, numthrds); // CD*C
      }
    }
  } else {
    for (int opind=0; opind<opvec1.size(); opind++) {    
      boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind)->getworkingrepresentation(loopblock);
      if (!otherBlock->get_op_array(CRE_DESCOMP).has_local_index(I, op1->get_orbs(0)))
	    return;
      const std::vector<boost::shared_ptr<SparseMatrix>>& opvec2 = otherBlock->has(DES_CRECOMP) ? \
          otherBlock->get_op_array(DES_CRECOMP).get_element(I, op1->get_orbs(0)) : \
          otherBlock->get_op_array(CRE_DESCOMP).get_element(I, op1->get_orbs(0));

      for (int opind2 = 0; opind2<opvec2.size(); opind2++) {
	boost::shared_ptr<SparseMatrix> op2 = opvec2.at(opind2)->getworkingrepresentation(otherBlock);
	vector<SpinQuantum> op2q = op2->get_deltaQuantum(), op1q = op1->get_deltaQuantum(), oq = o->get_deltaQuantum();
	int j2 = op2q[0].get_s().getirrep(), j1 = op1q[0].get_s().getirrep(), j21 = oq[0].get_s().getirrep();
	int l2 = (-op2q[0].get_symm()).getirrep(), l1 = op1q[0].get_symm().getirrep(), l21 = oq[0].get_symm().getirrep(), l3 = (-SymmetryOfOrb(I)).getirrep();
	double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((1+1+0+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
      if (NonabelianSym)
	factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());
	
	double parity = 1.0;
	if (otherBlock == b->get_rightBlock() && !otherBlock->has(DES_CRECOMP))
	  parity *= getCommuteParity(op1->get_deltaQuantum(0), -op2->get_deltaQuantum(0), o->get_deltaQuantum(0));
	else if (otherBlock == b->get_rightBlock())
	  parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0));
	
	parity *= TensorOp::getTransposeFactorCD(I, op1->get_orbs(0), j2, l2);
	
	//If we have all the operators we dont have to take transposes, useful for <bra|H|ket> evaluation
	if (!otherBlock->has(DES_CRECOMP)) {
	  SpinAdapted::operatorfunctions::TensorProduct(otherBlock, Transposeview(*op2), *op1, b, &(b->get_stateInfo()), o[ilock], factor*scale*parity, numthrds);	
        } else {
	  double parity = 1.0;
	  if (otherBlock == b->get_rightBlock())
	        parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0)); // doesn't depend on nelec
	      SpinAdapted::operatorfunctions::TensorProduct(otherBlock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], factor*scale*parity, numthrds);
	    }
      }
    }
  }
}

void SpinAdapted::opxop::dxcccomp(const SpinBlock* otherBlock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, int K, SparseMatrix* o, double scale)
{ 
  int ilock = 0;//omp_get_thread_num();
  int numthrds = 1;
  //int numthrds = dmrginp.thrds_per_node()[mpigetrank()];
  const SpinBlock* loopblock = (otherBlock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  for (int opind=0; opind<opvec1.size(); opind++) {

    //If we have all the operators we dont have to take transposes, useful for <bra|H|ket> evaluation
    if (loopblock->has(DES) ) {
      boost::shared_ptr<SparseMatrix> op1 = loopblock->get_op_array(DES).get_element(opvec1.at(opind)->get_orbs(0)).at(opind)->getworkingrepresentation(loopblock); // DES_j
      
      bool transpose = false;
      int k = K, i = op1->get_orbs(0); // P_{ij}=-P_{ji} so only one of them is stored --- P_{ij} where i>j
      if (k < i) { 
        k=i; i=K; transpose = true;
      }
      SpinQuantum iq = getSpinQuantum(i), kq = getSpinQuantum(k);
      
      if (!otherBlock->get_op_array(CRE_CRECOMP).has_local_index(k,i))
	    return;
      
      const std::vector<boost::shared_ptr<SparseMatrix> >& opvec2 = otherBlock->get_op_array(CRE_CRECOMP).get_element(k,i); // P_{ki}
      for (int opind2 = 0; opind2<opvec2.size(); opind2++) {
	    boost::shared_ptr<SparseMatrix> op2 = opvec2.at(opind2)->getworkingrepresentation(otherBlock);  // P_{ki}^\dagger
	    
	    SpinQuantum op2q = op2->get_deltaQuantum(0), op1q = op1->get_deltaQuantum(0), oq = o->get_deltaQuantum(0);
	    int j2 = op2q.get_s().getirrep(), j1 = op1q.get_s().getirrep(), j21 = oq.get_s().getirrep();
	    int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = (-SymmetryOfOrb(K)).getirrep();
	    double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
	    if (NonabelianSym)
	    factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());
	    
	    double parity = 1.0;

	    if (!transpose)  //this is because when you go from CC_{ij} to CC_{ji} there is a phase factor
	      parity *= getCommuteParity(iq, kq, op2->get_deltaQuantum(0)); 
	    if (loopblock == b->get_leftBlock()) //this is because you have CC_{ji} d_j 
	      parity*= getCommuteParity(-iq, op2->get_deltaQuantum(0), kq); 

	    SpinAdapted::operatorfunctions::TensorProduct(otherBlock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], parity*factor*scale, numthrds);
      }
    } else {
      boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind)->getworkingrepresentation(loopblock); // CRE_j
      
      bool transpose = false;
      int k = K, i = op1->get_orbs(0); // P_{ij}=-P_{ji} so only one of them is stored --- P_{ij} where i>j
      if (k < i) { 
        k=i; i=K; transpose = true;
      }
      SpinQuantum iq = getSpinQuantum(i), kq = getSpinQuantum(k);
      
      if (!otherBlock->get_op_array(DES_DESCOMP).has_local_index(k,i))
	    return;
      
      const std::vector<boost::shared_ptr<SparseMatrix> >& opvec2 = otherBlock->get_op_array(DES_DESCOMP).get_element(k,i); // P_{ki}
      for (int opind2 = 0; opind2<opvec2.size(); opind2++) {
	    Transposeview top = Transposeview(opvec2.at(opind2)->getworkingrepresentation(otherBlock));  // P_{ki}^\dagger
	    
	    SpinQuantum op2q = top.get_deltaQuantum(0), op1q = -op1->get_deltaQuantum(0), oq = o->get_deltaQuantum(0);
	    int j2 = op2q.get_s().getirrep(), j1 = op1q.get_s().getirrep(), j21 = oq.get_s().getirrep();
	    int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = (-SymmetryOfOrb(K)).getirrep();
	    double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
	    if (NonabelianSym)
	    factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());
	    
	    double parity = 1.0;
	    //this is because DD_ij^dag = - CC_ij for spin 0 and certain spatial irreps
	    parity*=TensorOp::getTransposeFactorDD(K, op1->get_orbs(0), j2, l2);
	    if (transpose)  //this is because when you go from CC_{ij} to CC_{ji} there is a phase factor
	      parity *= getCommuteParity(iq, kq, top.get_deltaQuantum(0)); 
	    if (loopblock == b->get_leftBlock()) //this is because you have CC_{ji} d_j 
	      parity*= getCommuteParity(-iq, top.get_deltaQuantum(0), kq); 
	    
	    //pout << k<<"  "<<i<<"  "<<factor<<"  "<<scale<<"  "<<parity<<endl;
	    SpinAdapted::operatorfunctions::TensorProduct(otherBlock, top, Transposeview(op1), b, &(b->get_stateInfo()), o[ilock], parity*factor*scale, numthrds);
      }
    }
  }
}


/********************************************
Formulas for making Cddcomp operators while blocking a block with a dot block
********************************************/

void SpinAdapted::opxop::dxcdcomp(const SpinBlock* otherBlock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, int I, SparseMatrix* o, double scale)
{
  int ilock = 0;//omp_get_thread_num();
  int numthrds = 1;
  const SpinBlock* loopblock = (otherBlock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  for (int opind=0; opind<opvec1.size(); opind++) {    
    boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind)->getworkingrepresentation(loopblock);
    if (!otherBlock->get_op_array(CRE_DESCOMP).has_local_index(max(I, op1->get_orbs(0)), min(I, op1->get_orbs(0)))) // we have c_J d_I *c_K d_L
      return;
    
    std::vector<boost::shared_ptr<SparseMatrix> > opvec2 ; 
    if (I > opvec1[0]->get_orbs(0) ) // opvec1 is DES
      opvec2 = otherBlock->get_op_array(CRE_DESCOMP).get_element(I, op1->get_orbs(0)); 
    else
      opvec2 = otherBlock->get_op_array(DES_CRECOMP).get_element(op1->get_orbs(0), I); 
    
    for (int opind2 = 0; opind2<opvec2.size(); opind2++) {
      boost::shared_ptr<SparseMatrix> op2 = opvec2.at(opind2)->getworkingrepresentation(otherBlock); // CD
      vector<SpinQuantum> op2q = op2->get_deltaQuantum(), op1q = op1->get_deltaQuantum(), oq = o->get_deltaQuantum(); // o is the resulted CCD
      int j2 = op2q[0].get_s().getirrep(), j1 = op1q[0].get_s().getirrep(), j21 = oq[0].get_s().getirrep();
      int l2 = op2q[0].get_symm().getirrep(), l1 = op1q[0].get_symm().getirrep(), l21 = oq[0].get_symm().getirrep(), l3 = (-SymmetryOfOrb(I)).getirrep();
      double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
      if (NonabelianSym)
      factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());
      
      double parity = 1.0;
      if (otherBlock == b->get_leftBlock())
	parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0)); // doesn't depend on nelec
      
      SpinAdapted::operatorfunctions::TensorProduct(otherBlock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], factor*scale*parity, numthrds); // CD*D
    }
  
  } 
}

void SpinAdapted::opxop::cxddcomp(const SpinBlock* otherBlock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, int K, SparseMatrix* o, double scale)
{ 
  int ilock = 0;//omp_get_thread_num();
  int numthrds = 1;
  //int numthrds = dmrginp.thrds_per_node()[mpigetrank()];
  const SpinBlock* loopblock = (otherBlock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  for (int opind=0; opind<opvec1.size(); opind++) {    
    
    boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind)->getworkingrepresentation(loopblock); // CRE_j
    
    bool transpose = false;
    int k = K, i = op1->get_orbs(0); // P_{ij}=-P_{ji} so only one of them is stored --- P_{ij} where i>j
    if (k < i) 
      { k=i; i=K; transpose = true;}
    SpinQuantum iq = getSpinQuantum(i), kq = getSpinQuantum(k);
    
    if (!otherBlock->get_op_array(DES_DESCOMP).has_local_index(k,i))
      return;
    
    const std::vector<boost::shared_ptr<SparseMatrix> >& opvec2 = otherBlock->get_op_array(DES_DESCOMP).get_element(k,i); // P_{ki}
    for (int opind2 = 0; opind2<opvec2.size(); opind2++) {
      boost::shared_ptr<SparseMatrix> op2 = opvec2.at(opind2)->getworkingrepresentation(otherBlock);  // P_{ki}^\dagger
      
      SpinQuantum op2q = op2->get_deltaQuantum(0), op1q = op1->get_deltaQuantum(0), oq = o->get_deltaQuantum(0);
      int j2 = op2q.get_s().getirrep(), j1 = op1q.get_s().getirrep(), j21 = oq.get_s().getirrep();
      int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = (-oq.get_symm()).getirrep(), l3 = (-SymmetryOfOrb(K)).getirrep();
      double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
      if (NonabelianSym)
	factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());
      
      double parity = 1.0;
      
      if (transpose)  //this is because when you go from CC_{ij} to CC_{ji} there is a phase factor
	parity *= getCommuteParity(iq, kq, -op2->get_deltaQuantum(0)); 
      if (loopblock == b->get_rightBlock()) //this is because you have CC_{ij} d_j 
	parity*= getCommuteParity(iq, op2->get_deltaQuantum(0), -kq); 

      SpinAdapted::operatorfunctions::TensorProduct(otherBlock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], parity*factor*scale, numthrds);

    }
  }
}


/********************************************
In nevpt2, to calculate V_a subspace, CDD_sum is calculatd by C*DD and D*CD.
CD and DD are complementary operator of C.
********************************************/


//**********************************************************************************************************

void SpinAdapted::opxop::cdd_cxddcomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, SparseMatrix* o)
{
  int ilock = omp_get_thread_num();
  int numthrds = 1;//MAX_THRD;
  const SpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  for (int opind1=0; opind1<opvec1.size(); opind1++) {
    boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind1)->getworkingrepresentation(loopblock); // CRE_i
    int i = op1->get_orbs(0);
    if (!otherblock->get_op_array(CDD_DES_DESCOMP).has_local_index(i))
      return;

    const std::vector<boost::shared_ptr<SparseMatrix> >& opvec2 = otherblock->get_op_array(CDD_DES_DESCOMP).get_element(i); // P_{ki}
    for (int opind2=0; opind2<opvec2.size(); opind2++) {
      boost::shared_ptr<SparseMatrix> op2 = opvec2.at(opind2)->getworkingrepresentation(otherblock);
      
	    SpinQuantum op1q = op1->get_deltaQuantum()[0], op2q = op2->get_deltaQuantum()[0], oq = -getSpinQuantum(b->nonactive_orb(0)); // o is the resulted CDD
	    int j2 = op2q.get_s().getirrep(), j1 = op1q.get_s().getirrep(), j21 = oq.get_s().getirrep();
	    int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = oq.get_symm().getirrep();
	    double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
	    factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());

      double parity = 1.0;
      if (otherblock == b->get_leftBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0));
      
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], factor*parity, numthrds);	    
    }

  }
}

void SpinAdapted::opxop::cdd_dxcdcomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, SparseMatrix* o)
{
  int ilock = omp_get_thread_num();
  int numthrds = 1;//MAX_THRD;
  const SpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  for (int opind1=0; opind1<opvec1.size(); opind1++) {
    boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind1)->getworkingrepresentation(loopblock); // CRE_i
    int i = op1->get_orbs(0);
    if (!otherblock->get_op_array(CDD_CRE_DESCOMP).has_local_index(i))
      return;

    const std::vector<boost::shared_ptr<SparseMatrix> >& opvec2 = otherblock->get_op_array(CDD_CRE_DESCOMP).get_element(i); // P_{ki}
    for (int opind2=0; opind2<opvec2.size(); opind2++) {
      boost::shared_ptr<SparseMatrix> op2 = opvec2.at(opind2)->getworkingrepresentation(otherblock);

	    SpinQuantum op1q = op1->get_deltaQuantum()[0], op2q = op2->get_deltaQuantum()[0], oq = -getSpinQuantum(b->nonactive_orb(0)); // o is the resulted CDD
	    int j2 = op2q.get_s().getirrep(), j1 = op1q.get_s().getirrep(), j21 = oq.get_s().getirrep();
	    int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = oq.get_symm().getirrep();
	    double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
	    factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());

    
      double parity = 1.0;
      if (otherblock == b->get_leftBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0));
      
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], factor*parity, numthrds);	    
    }
  }
}

void SpinAdapted::opxop::cdd_cxddcomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q)
{
  int ilock = omp_get_thread_num();
  int numthrds = 1;//MAX_THRD;
  const SpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
    
  for (int opind1=0; opind1<opvec1.size(); opind1++) {
    boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind1)->getworkingrepresentation(loopblock);
    int i = op1->get_orbs(0);
    if (!otherblock->get_op_array(CDD_DES_DESCOMP).has_local_index(i))
      return;
    const std::vector<boost::shared_ptr<SparseMatrix> >& opvec2 = otherblock->get_op_array(CDD_DES_DESCOMP).get_element(i); // P_{ki}
    for (int opind2=0; opind2<opvec2.size(); opind2++) {
      boost::shared_ptr<SparseMatrix> op2 = opvec2.at(opind2)->getworkingrepresentation(otherblock);

	    SpinQuantum op1q = op1->get_deltaQuantum()[0], op2q = op2->get_deltaQuantum()[0], oq = -getSpinQuantum(b->nonactive_orb(0)); // o is the resulted CDD
	    int j2 = op2q.get_s().getirrep(), j1 = op1q.get_s().getirrep(), j21 = oq.get_s().getirrep();
	    int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = oq.get_symm().getirrep();
	    double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
	    factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());

      double parity = 1.0;
      if (otherblock == b->get_leftBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), q);

      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v[ilock], q, factor*parity);
    }
  }
}

void SpinAdapted::opxop::cdd_dxcdcomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q)
{
  int ilock = omp_get_thread_num();
  int numthrds = 1;//MAX_THRD;
  const SpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
    
  for (int opind1=0; opind1<opvec1.size(); opind1++) {
    boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind1)->getworkingrepresentation(loopblock);
    int i = op1->get_orbs(0);
    if (!otherblock->get_op_array(CDD_CRE_DESCOMP).has_local_index(i))
      return;
    const std::vector<boost::shared_ptr<SparseMatrix> >& opvec2 = otherblock->get_op_array(CDD_CRE_DESCOMP).get_element(i); // P_{ki}
    for (int opind2=0; opind2<opvec2.size(); opind2++) {
      boost::shared_ptr<SparseMatrix> op2 = opvec2.at(opind2)->getworkingrepresentation(otherblock);

	    SpinQuantum op1q = op1->get_deltaQuantum()[0], op2q = op2->get_deltaQuantum()[0], oq = -getSpinQuantum(b->nonactive_orb(0)); // o is the resulted CDD
	    int j2 = op2q.get_s().getirrep(), j1 = op1q.get_s().getirrep(), j21 = oq.get_s().getirrep();
	    int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = oq.get_symm().getirrep();
	    double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
	    factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());


      double parity = 1.0;
      if (otherblock == b->get_leftBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), q);
      //if (loopblock == b->get_leftBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), q);

      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v[ilock], q, factor*parity);

    }
  }
}

//**********************************************************************************************************

void SpinAdapted::opxop::ccd_dxcccomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, SparseMatrix* o)
{
  int ilock = omp_get_thread_num();
  int numthrds = 1;//MAX_THRD;
  const SpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  for (int opind1=0; opind1<opvec1.size(); opind1++) {
    boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind1)->getworkingrepresentation(loopblock); // DES_i
    int i = op1->get_orbs(0);
    if (!otherblock->get_op_array(CCD_CRE_CRECOMP).has_local_index(i))
      return;

    const std::vector<boost::shared_ptr<SparseMatrix> >& opvec2 = otherblock->get_op_array(CCD_CRE_CRECOMP).get_element(i); // P_{ki}
    for (int opind2=0; opind2<opvec2.size(); opind2++) {
      boost::shared_ptr<SparseMatrix> op2 = opvec2.at(opind2)->getworkingrepresentation(otherblock);
      
      double scale = 1.0;
      double parity = 1.0;
      if (otherblock == b->get_rightBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0));
      
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], scale*parity, numthrds);	    
    }

  }
}

void SpinAdapted::opxop::ccd_cxcdcomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, SparseMatrix* o)
{
  int ilock = omp_get_thread_num();
  int numthrds = 1;//MAX_THRD;
  const SpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  for (int opind1=0; opind1<opvec1.size(); opind1++) {
    boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind1)->getworkingrepresentation(loopblock); // CRE_i
    int i = op1->get_orbs(0);
    if (!otherblock->get_op_array(CCD_CRE_DESCOMP).has_local_index(i))
      return;

    const std::vector<boost::shared_ptr<SparseMatrix> >& opvec2 = otherblock->get_op_array(CCD_CRE_DESCOMP).get_element(i); // P_{ki}
    for (int opind2=0; opind2<opvec2.size(); opind2++) {
      boost::shared_ptr<SparseMatrix> op2 = opvec2.at(opind2)->getworkingrepresentation(otherblock);
    
      double scale = 1.0;
      double parity = 1.0;
      if (otherblock == b->get_rightBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0));
      
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], scale*parity, numthrds);	    
    }
  }
}

void SpinAdapted::opxop::ccd_dxcccomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q)
{
  int ilock = omp_get_thread_num();
  int numthrds = 1;//MAX_THRD;
  const SpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
    
  for (int opind1=0; opind1<opvec1.size(); opind1++) {
    boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind1)->getworkingrepresentation(loopblock);
    int i = op1->get_orbs(0);
    if (!otherblock->get_op_array(CCD_CRE_CRECOMP).has_local_index(i))
      return;
    const std::vector<boost::shared_ptr<SparseMatrix> >& opvec2 = otherblock->get_op_array(CCD_CRE_CRECOMP).get_element(i); // P_{ki}
    for (int opind2=0; opind2<opvec2.size(); opind2++) {
      boost::shared_ptr<SparseMatrix> op2 = opvec2.at(opind2)->getworkingrepresentation(otherblock);

	    SpinQuantum op1q = op1->get_deltaQuantum()[0], op2q = op2->get_deltaQuantum()[0], oq = getSpinQuantum(b->nonactive_orb(0)); // o is the resulted CCD
	    int j2 = op2q.get_s().getirrep(), j1 = op1q.get_s().getirrep(), j21 = oq.get_s().getirrep();
	    int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = oq.get_symm().getirrep();
	    double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
	    factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());

      double parity = 1.0;
      if (otherblock == b->get_rightBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), q);

      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v[ilock], q, factor*parity);
    }
  }
}

void SpinAdapted::opxop::ccd_cxcdcomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q)
{
  int ilock = omp_get_thread_num();
  int numthrds = 1;//MAX_THRD;
  const SpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
    
  for (int opind1=0; opind1<opvec1.size(); opind1++) {
    boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind1)->getworkingrepresentation(loopblock);
    int i = op1->get_orbs(0);
    if (!otherblock->get_op_array(CCD_CRE_DESCOMP).has_local_index(i))
      return;
    const std::vector<boost::shared_ptr<SparseMatrix> >& opvec2 = otherblock->get_op_array(CCD_CRE_DESCOMP).get_element(i); // P_{ki}
    for (int opind2=0; opind2<opvec2.size(); opind2++) {
      boost::shared_ptr<SparseMatrix> op2 = opvec2.at(opind2)->getworkingrepresentation(otherblock);

	    SpinQuantum op1q = op1->get_deltaQuantum()[0], op2q = op2->get_deltaQuantum()[0], oq = -getSpinQuantum(b->nonactive_orb(0)); // o is the resulted CDD
	    int j2 = op2q.get_s().getirrep(), j1 = op1q.get_s().getirrep(), j21 = oq.get_s().getirrep();
	    int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = oq.get_symm().getirrep();
	    double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
	    factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());


      double parity = 1.0;
      if (otherblock == b->get_rightBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), q);
      //if (loopblock == b->get_leftBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), q);

      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v[ilock], q, factor*parity);

    }
  }
}

