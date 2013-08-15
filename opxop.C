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

//using namespace operatorfunctions;


/********************************************
Formulas for making hamiltonian matrix while blocking a block with a dot block
********************************************/


void SpinAdapted::opxop::cdxcdcomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, SparseMatrix* o)
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

    boost::shared_ptr<SparseMatrix> op3 = otherblock->get_op_array(CRE_DESCOMP).get_element(i, j).at(opind)->getworkingrepresentation(otherblock);
    double factor = 1.0;
    if (otherblock == b->get_leftBlock())
      factor = getCommuteParity(op1->get_deltaQuantum(), op3->get_deltaQuantum(), o->get_deltaQuantum());

    SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op3, *op1, b, &(b->get_stateInfo()), o[ilock], factor, numthrds);
    if (i != j) {
      factor = 1.0;
      if (otherblock == b->get_rightBlock())
	factor = getCommuteParity(-op1->get_deltaQuantum(), -op3->get_deltaQuantum(), o->get_deltaQuantum());
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, Transposeview(*op3), Transposeview(*op1), b, &(b->get_stateInfo()), o[ilock], factor, numthrds);
    }
  }
}


void SpinAdapted::opxop::ddxcccomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, SparseMatrix* o)
{
  int ilock = omp_get_thread_num();
  int numthrds = 1;//MAX_THRD;
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
      parity = getCommuteParity(op1->get_deltaQuantum(), op2->get_deltaQuantum(), o->get_deltaQuantum());

    SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], parity*factor, numthrds);

    Transposeview top1 = Transposeview(*op1);
    Transposeview top2 = Transposeview(*op2);
    parity = 1.0;
    if (otherblock == b->get_rightBlock())
      parity = getCommuteParity(-op1->get_deltaQuantum(), -op2->get_deltaQuantum(), o->get_deltaQuantum());

    SpinAdapted::operatorfunctions::TensorProduct(otherblock, top2, top1, b, &(b->get_stateInfo()), o[ilock], parity*factor, numthrds);  
  }
}



void SpinAdapted::opxop::cxcddcomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, SparseMatrix* o)
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

    if (otherblock == b->get_rightBlock())
      parity = getCommuteParity(-op1->get_deltaQuantum(), op2->get_deltaQuantum(), o->get_deltaQuantum());

    SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, top1, b, &(b->get_stateInfo()), o[ilock], scale*parity, numthrds);
    Transposeview top2 = Transposeview(op2);

    if (otherblock == b->get_leftBlock())
      parity = getCommuteParity(-op2->get_deltaQuantum(), op1->get_deltaQuantum(), o->get_deltaQuantum());

    SpinAdapted::operatorfunctions::TensorProduct(otherblock, top2, *op1, b, &(b->get_stateInfo()), o[ilock], scale*parity, numthrds);	    
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
  SpinQuantum hq(0,0,IrrepSpace(0));
  const SpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
    

  for (int opind=0; opind<opvec1.size(); opind++) {
    boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind)->getworkingrepresentation(loopblock);
    int i = op1->get_orbs(0);
    int j = op1->get_orbs(1);
    if (!otherblock->get_op_array(CRE_DESCOMP).has_local_index(i,j))
      return;
    boost::shared_ptr<SparseMatrix> op3 = otherblock->get_op_array(CRE_DESCOMP).get_element(i, j).at(opind)->getworkingrepresentation(otherblock);
    double factor = 1.0;
    SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op3, *op1, b, c, v[ilock], hq, factor);
    if (i != j)
      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, Transposeview(*op3), Transposeview(*op1), b, c, v[ilock], hq, factor);
  }
}


void SpinAdapted::opxop::ddxcccomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q)
{
  int ilock = omp_get_thread_num();
  int numthrds = 1;//MAX_THRD;
  SpinQuantum hq(0,0,IrrepSpace(0));
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
      parity = getCommuteParity(op1->get_deltaQuantum(), op2->get_deltaQuantum(), hq);
    SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v[ilock], hq, factor*parity);
    
    Transposeview top1 = Transposeview(*op1);
    Transposeview top2 = Transposeview(*op2);
    parity = 1.0;
    if (otherblock == b->get_rightBlock())
      parity = getCommuteParity(-op1->get_deltaQuantum(), -op2->get_deltaQuantum(), hq);
    SpinAdapted::operatorfunctions::TensorMultiply(otherblock, top2, top1, b, c, v[ilock], hq, factor*parity);  
  }
}



void SpinAdapted::opxop::cxcddcomp(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, Wavefunction& c, Wavefunction* v, const SpinQuantum& q)
{
  int ilock = omp_get_thread_num();
  int numthrds = 1;//MAX_THRD;
  SpinQuantum hq(0,0,IrrepSpace(0));
  const SpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  for (int opind=0; opind<opvec1.size(); opind++) {
    boost::shared_ptr<SparseMatrix> o1 = opvec1.at(opind);
    int i = o1->get_orbs(0);
    if (!otherblock->get_op_array(CRE_CRE_DESCOMP).has_local_index(i))
      return;
    boost::shared_ptr<SparseMatrix> op1 = o1->getworkingrepresentation(loopblock);
    Transposeview top1 = Transposeview(*op1);
    boost::shared_ptr<SparseMatrix> op2 = otherblock->get_op_array(CRE_CRE_DESCOMP).get_element(i).at(opind)->getworkingrepresentation(otherblock);

    double scale = 1.0;
    double parity = 1.0;
    if (otherblock == b->get_rightBlock())
      parity = getCommuteParity(op1->get_deltaQuantum(), op2->get_deltaQuantum(), hq);

    SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, top1, b, c, v[ilock], hq, scale*parity);
    Transposeview top2 = Transposeview(*op2);    
    //if (otherblock == b->get_leftBlock())
    //parity = getCommuteParity(op1->get_deltaQuantum(), op2->get_deltaQuantum(), hq);

    SpinAdapted::operatorfunctions::TensorMultiply(otherblock, top2, *op1, b, c, v[ilock], hq, scale*parity);	    
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
  SpinQuantum q(0,0,IrrepSpace(0));
  const SpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
    

  for (int opind=0; opind<opvec1.size(); opind++) {
    boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind)->getworkingrepresentation(loopblock);
    int i = op1->get_orbs(0);
    int j = op1->get_orbs(1);
    if (!otherblock->get_op_array(CRE_DESCOMP).has_local_index(i,j))
      return;

    boost::shared_ptr<SparseMatrix> op3 = otherblock->get_op_array(CRE_DESCOMP).get_element(i, j).at(opind)->getworkingrepresentation(otherblock);
    double factor = 1.0;
    SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op3, *op1, b, &(b->get_stateInfo()), e[ilock], factor);
    if (i != j)
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, Transposeview(*op3), Transposeview(*op1), b, &(b->get_stateInfo()), e[ilock], factor);
  }
}


void SpinAdapted::opxop::ddxcccomp_d(const SpinBlock* otherblock, std::vector<boost::shared_ptr<SparseMatrix> >& opvec1, const SpinBlock* b, DiagonalMatrix* e)
{
  int ilock = omp_get_thread_num();
  int numthrds = 1;//MAX_THRD;
  SpinQuantum q(0,0,IrrepSpace(0));
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
  SpinQuantum q(0,0,IrrepSpace(0));
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

  if (opvec1[0]->get_orbs(0) >= I)
  {
    for (int opind=0; opind<opvec1.size(); opind++) {    
      boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind)->getworkingrepresentation(loopblock);
      if (!otherBlock->get_op_array(CRE_DESCOMP).has_local_index(op1->get_orbs(0), I))
	return;

      const std::vector<boost::shared_ptr<SparseMatrix> >& opvec2 = otherBlock->get_op_array(CRE_DESCOMP).get_element(op1->get_orbs(0), I);
      for (int opind2 = 0; opind2<opvec2.size(); opind2++) {
	boost::shared_ptr<SparseMatrix> op2 = opvec2.at(opind2)->getworkingrepresentation(otherBlock);

	SpinQuantum op2q = op2->get_deltaQuantum(), op1q = op1->get_deltaQuantum(), oq = o->get_deltaQuantum();
	int j2 = op2q.get_s(), j1 = op1q.get_s(), j21 = oq.get_s();
	int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = (-SymmetryOfSpatialOrb(I)).getirrep();
	double factor = pow(-1.0, static_cast<int>((2+op2q.get_s())/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1));
	factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());

	double parity = 1.0;
	if (otherBlock == b->get_rightBlock())
	  parity *= getCommuteParity(op1->get_deltaQuantum(), op2->get_deltaQuantum(), o->get_deltaQuantum());
	factor*= parity;

	SpinAdapted::operatorfunctions::TensorProduct(otherBlock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], factor*scale, numthrds);
      }
    }
  }
  else
  {
    for (int opind=0; opind<opvec1.size(); opind++) {    
      boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind)->getworkingrepresentation(loopblock);
      if (!otherBlock->get_op_array(CRE_DESCOMP).has_local_index(I, op1->get_orbs(0)))
	return;
      const std::vector<boost::shared_ptr<SparseMatrix> >& opvec2 = otherBlock->get_op_array(CRE_DESCOMP).get_element(I, op1->get_orbs(0));

      for (int opind2 = 0; opind2<opvec2.size(); opind2++) {
	boost::shared_ptr<SparseMatrix> op2 = opvec2.at(opind2)->getworkingrepresentation(otherBlock);

	SpinQuantum op2q = op2->get_deltaQuantum(), op1q = op1->get_deltaQuantum(), oq = o->get_deltaQuantum();
	int j2 = op2q.get_s(), j1 = op1q.get_s(), j21 = oq.get_s();
	int l2 = (-op2q.get_symm()).getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = (-SymmetryOfSpatialOrb(I)).getirrep();
	double factor = pow(-1.0, static_cast<int>((1+1+0+op2q.get_s())/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1));
	factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());


	double parity = 1.0;
	if (otherBlock == b->get_rightBlock())
	  parity *= getCommuteParity(op1->get_deltaQuantum(), -op2->get_deltaQuantum(), o->get_deltaQuantum());

	SpinQuantum iq = getSpinQuantum(I), kq = getSpinQuantum(opvec1[0]->get_orbs(0));
	parity *= getCommuteParity(iq, -kq, -op2q); 

	factor*= parity;

	SpinAdapted::operatorfunctions::TensorProduct(otherBlock, Transposeview(*op2), *op1, b, &(b->get_stateInfo()), o[ilock], factor*scale, numthrds);
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
    boost::shared_ptr<SparseMatrix> op1 = opvec1.at(opind)->getworkingrepresentation(loopblock);

    bool transpose = true;
    int k = K, i = op1->get_orbs(0);
    if (k < i) 
      { k=i; i=K; transpose = false;}
    SpinQuantum iq = getSpinQuantum(i), kq = getSpinQuantum(k);
    
    if (!otherBlock->get_op_array(DES_DESCOMP).has_local_index(k,i))
	return;

    const std::vector<boost::shared_ptr<SparseMatrix> >& opvec2 = otherBlock->get_op_array(DES_DESCOMP).get_element(k,i);

    for (int opind2 = 0; opind2<opvec2.size(); opind2++) {
      
      Transposeview top = Transposeview(opvec2.at(opind2)->getworkingrepresentation(otherBlock));

      SpinQuantum op2q = top.get_deltaQuantum(), op1q = -op1->get_deltaQuantum(), oq = o->get_deltaQuantum();
      int j2 = op2q.get_s(), j1 = op1q.get_s(), j21 = oq.get_s();
      int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = (-SymmetryOfSpatialOrb(K)).getirrep();
      double factor = pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1));
      factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());


      double parity = 1.0;
      if (transpose)
	parity*= getCommuteParity(iq, kq, top.get_deltaQuantum()); 

      if (loopblock == b->get_leftBlock()) 
	parity*= getCommuteParity(-iq, top.get_deltaQuantum(), kq); 
	

      SpinAdapted::operatorfunctions::TensorProduct(otherBlock, top, Transposeview(op1), b, &(b->get_stateInfo()), o[ilock], parity*factor*scale, numthrds);

    }
  }
}
//**********************************************************************************************************

