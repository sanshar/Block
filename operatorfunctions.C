/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "operatorfunctions.h"
#include "wavefunction.h"
#include "couplingCoeffs.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "pario.h"


void SpinAdapted::operatorfunctions::TensorTrace(const SpinBlock *ablock, const Baseoperator<Matrix>& a, const SpinBlock *cblock, const StateInfo *cstateinfo, Baseoperator<Matrix>& c, double scale, int num_thrds)
{

  if (fabs(scale) < TINY) return;
  assert (a.get_initialised() && c.get_initialised());
#ifdef _OPENMP
#pragma omp parallel default(shared) num_threads(num_thrds) 
#endif
  {
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
  for (int cq = 0; cq < c.nrows(); ++cq)
    for (int cqprime = 0; cqprime < c.ncols(); ++cqprime)
      if (c.allowed(cq, cqprime)) 
	TensorTraceElement(ablock, a, cblock, cstateinfo, c, c.operator_element(cq, cqprime), cq, cqprime, scale);
  }
}


void SpinAdapted::operatorfunctions::TensorTraceElement(const SpinBlock *ablock, const Baseoperator<Matrix>& a, const SpinBlock *cblock, const StateInfo *cstateinfo, Baseoperator<Matrix>& c, Matrix& cel, int cq, int cqprime)
{
  TensorTraceElement(ablock, a, cblock, cstateinfo, c, cel, cq, cqprime, 1.);
}

void SpinAdapted::operatorfunctions::TensorTraceElement(const SpinBlock *ablock, const Baseoperator<Matrix>& a, const SpinBlock *cblock, const StateInfo *cstateinfo, Baseoperator<Matrix>& c, Matrix& cel, int cq, int cqprime, double scale)
{
  if (fabs(scale) < TINY) return;
  assert(c.allowed(cq, cqprime));
    
  int aq, aqprime, bq, bqprime, bstates;
  const char conjC = (ablock == cblock->get_leftBlock()) ? 'n' : 't';

  const std::vector<int> oldToNewI = cstateinfo->oldToNewState.at(cq);
  const std::vector<int> oldToNewJ = cstateinfo->oldToNewState.at(cqprime);

  const StateInfo* rS = cstateinfo->rightStateInfo, *lS = cstateinfo->leftStateInfo;
  int rowstride =0, colstride = 0;

  for (int oldi =0; oldi < oldToNewI.size(); oldi++) {
    colstride = 0;
    for (int oldj = 0; oldj < oldToNewJ.size(); oldj++)
    {
      if (conjC == 'n')
      {
	aq = cstateinfo->leftUnMapQuanta[ oldToNewI[oldi] ];
	aqprime = cstateinfo->leftUnMapQuanta[ oldToNewJ[oldj] ];
	bq = cstateinfo->rightUnMapQuanta[ oldToNewI[oldi] ];
	bqprime = cstateinfo->rightUnMapQuanta[ oldToNewJ[oldj] ];
	bstates = cstateinfo->rightStateInfo->getquantastates(bq);
      }
      else 
      {
	aq = cstateinfo->rightUnMapQuanta[ oldToNewI[oldi] ];
	aqprime = cstateinfo->rightUnMapQuanta[ oldToNewJ[oldj] ];
	bq = cstateinfo->leftUnMapQuanta[ oldToNewI[oldi] ];
	bqprime = cstateinfo->leftUnMapQuanta[ oldToNewJ[oldj] ];
	bstates = cstateinfo->leftStateInfo->getquantastates(bq);
      }
      
      if (a.allowed(aq, aqprime) && (bq == bqprime))
      {
	DiagonalMatrix unitMatrix(bstates);
	unitMatrix = 1.;
	Matrix unity(bstates, bstates);
	unity = unitMatrix;

	if (conjC == 'n')
	{
	  double scaleb = dmrginp.get_ninej()(lS->quanta[aqprime].get_s().getirrep() , rS->quanta[bqprime].get_s().getirrep(), cstateinfo->quanta[cqprime].get_s().getirrep(), 
				a.get_spin().getirrep(), 0, c.get_spin().getirrep(),
				lS->quanta[aq].get_s().getirrep() , rS->quanta[bq].get_s().getirrep(), cstateinfo->quanta[cq].get_s().getirrep());

	  scaleb *= Symmetry::spatial_ninej(lS->quanta[aqprime].get_symm().getirrep() , rS->quanta[bqprime].get_symm().getirrep(), cstateinfo->quanta[cqprime].get_symm().getirrep(), 
			       a.get_symm().getirrep(), 0, c.get_symm().getirrep(),
			       lS->quanta[aq].get_symm().getirrep() , rS->quanta[bq].get_symm().getirrep(), cstateinfo->quanta[cq].get_symm().getirrep());

	  MatrixTensorProduct (a.operator_element(aq, aqprime), a.conjugacy(), scale, unity, 'n', scaleb, 
	  	       cel, rowstride, colstride);
	}
	else {
	  double scaleb = dmrginp.get_ninej()(lS->quanta[bqprime].get_s().getirrep(), rS->quanta[aqprime].get_s().getirrep() , cstateinfo->quanta[cqprime].get_s().getirrep(), 
				0, a.get_spin().getirrep(), c.get_spin().getirrep(),
				lS->quanta[bq].get_s().getirrep(), rS->quanta[aq].get_s().getirrep() , cstateinfo->quanta[cq].get_s().getirrep());
	  scaleb *= Symmetry::spatial_ninej(lS->quanta[bqprime].get_symm().getirrep() , rS->quanta[aqprime].get_symm().getirrep(), cstateinfo->quanta[cqprime].get_symm().getirrep(), 
			       0, a.get_symm().getirrep(), c.get_symm().getirrep(),
			       lS->quanta[bq].get_symm().getirrep() , rS->quanta[aq].get_symm().getirrep(), cstateinfo->quanta[cq].get_symm().getirrep());
	  if (a.get_fermion() && IsFermion (cstateinfo->leftStateInfo->quanta[bqprime]) ) scaleb *= -1.;
	  MatrixTensorProduct (unity, 'n', scaleb, a.operator_element(aq, aqprime), a.conjugacy(), scale, 
	  		       cel, rowstride, colstride);
	}
      }
      colstride += cstateinfo->unCollectedStateInfo->quantaStates[ oldToNewJ[oldj] ];

    }
    rowstride += cstateinfo->unCollectedStateInfo->quantaStates[ oldToNewI[oldi] ];
    
  }
}


void SpinAdapted::operatorfunctions::Product (const SpinBlock *ablock, const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, Baseoperator<Matrix>& c, double scale)
{
  const StateInfo* astate = &ablock->get_stateInfo(); 
  if (fabs(scale) < TINY) return;
  int rows = c.nrows();
  for (int cq = 0; cq < rows; ++cq)
    for (int cqprime = 0; cqprime < rows; ++cqprime)
      if (c.allowed(cq, cqprime))
	for (int aprime = 0; aprime < rows; aprime++)
	  if (a.allowed(cq, aprime) && b.allowed(aprime, cqprime))
	  {
	    int apj = astate->quanta[aprime].get_s().getirrep(), cqj = astate->quanta[cq].get_s().getirrep(), cqpj = astate->quanta[cqprime].get_s().getirrep();
	    double factor = a.get_scaling(astate->quanta[cq], astate->quanta[aprime]);
	    factor *= b.get_scaling(astate->quanta[aprime], astate->quanta[cqprime]);
      if(dmrginp.spinAdapted()){

	    factor *= racah(cqpj, b.get_spin().getirrep(), cqj, a.get_spin().getirrep(), apj, c.get_spin().getirrep()) * pow( (1.0*c.get_spin().getirrep()+1.0)*(1.0*apj+1.0), 0.5 )
	            *pow(-1.0, static_cast<int>((b.get_spin().getirrep()+a.get_spin().getirrep()-c.get_spin().getirrep())/2.0));
      }
	    MatrixMultiply(a.operator_element(cq, aprime), a.conjugacy(), b.operator_element(aprime, cqprime), b.conjugacy(),
			   c.operator_element(cq, cqprime), scale*factor, 1.0);

	  }
}


void SpinAdapted::operatorfunctions::TensorProduct (const SpinBlock *ablock, const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, const SpinBlock *cblock, const StateInfo *cstateinfo, Baseoperator<Matrix>& c, double scale, int num_thrds)
{
  if (fabs(scale) < TINY) return;
  int rows = c.nrows();
  int cols = c.ncols();
#ifdef _OPENMP
#pragma omp parallel default(shared) num_threads(num_thrds) 
#endif
  {
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
  for (int cq = 0; cq < rows; ++cq)
    for (int cqprime = 0; cqprime < cols; ++cqprime)
      if (c.allowed(cq, cqprime)) {
	TensorProductElement(ablock, a, b, cblock, cstateinfo, c, c.operator_element(cq, cqprime), cq, cqprime, scale);
      }
  }
}


void SpinAdapted::operatorfunctions::TensorProductElement(const SpinBlock *ablock, const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, const SpinBlock *cblock, const StateInfo *cstateinfo, Baseoperator<Matrix>& c, Matrix& cel, int cq, int cqprime, double scale)
{
  //cstateinfo is not used
  //This function can be used for different bra and ket stateinfo.
  if (fabs(scale) < TINY) return;
  assert (a.get_initialised());
  assert (b.get_initialised());
  assert (c.get_initialised());
 
  const StateInfo *ketstateinfo = &cblock->get_ketStateInfo(), 
    *brastateinfo = &cblock->get_braStateInfo();

  const SpinBlock* bblock = (cblock->get_leftBlock() == ablock) ? cblock->get_rightBlock() : cblock->get_leftBlock();
 
  const std::vector<int>& oldToNewI = brastateinfo->oldToNewState.at(cq);
  const std::vector<int>& oldToNewJ = ketstateinfo->oldToNewState.at(cqprime);
 
  const char conjC = (cblock->get_leftBlock() == ablock) ? 'n' : 't';
  
  const StateInfo* lbraS = brastateinfo->leftStateInfo, *rbraS = brastateinfo->rightStateInfo;
  const StateInfo* lketS = ketstateinfo->leftStateInfo, *rketS = ketstateinfo->rightStateInfo;
  int rowstride = 0, colstride = 0;

  int aq, aqprime, bq, bqprime;

  //pout << "old to new size "<<oldToNewI.size()<<" "<<oldToNewJ.size()<<endl;
  for (int oldi =0; oldi < oldToNewI.size(); oldi++) {
    colstride = 0;
    for (int oldj = 0; oldj < oldToNewJ.size(); oldj++)
    {
      if (conjC == 'n')
      {
	aq = brastateinfo->leftUnMapQuanta[ oldToNewI[oldi] ];
	aqprime = ketstateinfo->leftUnMapQuanta[ oldToNewJ[oldj] ];
	bq = brastateinfo->rightUnMapQuanta[ oldToNewI[oldi] ];
	bqprime = ketstateinfo->rightUnMapQuanta[ oldToNewJ[oldj] ];
      }
      else 
      {
	aq = brastateinfo->rightUnMapQuanta[ oldToNewI[oldi] ];
	aqprime = ketstateinfo->rightUnMapQuanta[ oldToNewJ[oldj] ];
	bq = brastateinfo->leftUnMapQuanta[ oldToNewI[oldi] ];
	bqprime = ketstateinfo->leftUnMapQuanta[ oldToNewJ[oldj] ];
      }
  
      Real scaleA = scale;
      Real scaleB = 1.0;
      if (a.allowed(aq, aqprime) && b.allowed(bq, bqprime))
      {
	if (conjC == 'n')
	{
	  scaleB = dmrginp.get_ninej()(lketS->quanta[aqprime].get_s().getirrep() , rketS->quanta[bqprime].get_s().getirrep(), ketstateinfo->quanta[cqprime].get_s().getirrep(), 
			 a.get_spin().getirrep(), b.get_spin().getirrep(), c.get_spin().getirrep(),
			 lbraS->quanta[aq].get_s().getirrep() , rbraS->quanta[bq].get_s().getirrep(), brastateinfo->quanta[cq].get_s().getirrep());
	  scaleB *= Symmetry::spatial_ninej(lketS->quanta[aqprime].get_symm().getirrep() , rketS->quanta[bqprime].get_symm().getirrep(), ketstateinfo->quanta[cqprime].get_symm().getirrep(), 
			       a.get_symm().getirrep(), b.get_symm().getirrep(), c.get_symm().getirrep(),
			       lbraS->quanta[aq].get_symm().getirrep() , rbraS->quanta[bq].get_symm().getirrep(), brastateinfo->quanta[cq].get_symm().getirrep());
	  scaleB *= b.get_scaling(rbraS->quanta[bq], rketS->quanta[bqprime]);
	  scaleA *= a.get_scaling(lbraS->quanta[aq], lketS->quanta[aqprime]);
	  if (b.get_fermion() && IsFermion (ketstateinfo->leftStateInfo->quanta [aqprime])) scaleB *= -1;
	  MatrixTensorProduct (a.operator_element(aq, aqprime), a.conjugacy(), scaleA, 
			       b.operator_element(bq, bqprime), b.conjugacy(), scaleB, cel,rowstride, colstride);
	}
	else
	{
	  scaleB = dmrginp.get_ninej()(lketS->quanta[bqprime].get_s().getirrep(), rketS->quanta[aqprime].get_s().getirrep() , ketstateinfo->quanta[cqprime].get_s().getirrep(), 
			 b.get_spin().getirrep(), a.get_spin().getirrep(), c.get_spin().getirrep(),
			 lbraS->quanta[bq].get_s().getirrep(), rbraS->quanta[aq].get_s().getirrep() , brastateinfo->quanta[cq].get_s().getirrep());
	  scaleB *= Symmetry::spatial_ninej(lketS->quanta[bqprime].get_symm().getirrep() , rketS->quanta[aqprime].get_symm().getirrep(), ketstateinfo->quanta[cqprime].get_symm().getirrep(), 
			       b.get_symm().getirrep(), a.get_symm().getirrep(), c.get_symm().getirrep(),
			       lbraS->quanta[bq].get_symm().getirrep() , rbraS->quanta[aq].get_symm().getirrep(), brastateinfo->quanta[cq].get_symm().getirrep());
	  scaleB *= b.get_scaling(lbraS->quanta[bq], lketS->quanta[bqprime]);
	  scaleA *= a.get_scaling(rbraS->quanta[aq], rketS->quanta[aqprime]);
	  if (a.get_fermion() && IsFermion (ketstateinfo->leftStateInfo->quanta[bqprime]) ) scaleB *= -1.;
	  
	  MatrixTensorProduct (b.operator_element(bq, bqprime), b.conjugacy(), scaleB, 
			       a.operator_element(aq, aqprime), a.conjugacy(), scaleA, cel, rowstride, colstride);
	}
      }
      colstride += ketstateinfo->unCollectedStateInfo->quantaStates[ oldToNewJ[oldj] ];

    }
    rowstride += brastateinfo->unCollectedStateInfo->quantaStates[ oldToNewI[oldi] ];
    
  }

  
}

/*
void SpinAdapted::operatorfunctions::TensorMultiply(const SpinBlock *ablock, const Baseoperator<Matrix>& a, const SpinBlock *cblock, Wavefunction& c, Wavefunction& v, const SpinQuantum dQ, double scale, int num_thrds)
{
  // cannot be used for situation with different bra and ket
  const int leftOpSz = cblock->get_leftBlock()->get_stateInfo().quanta.size ();
  const int rightOpSz = cblock->get_rightBlock()->get_stateInfo().quanta.size ();

  const StateInfo* rS = cblock->get_stateInfo().rightStateInfo, *lS = cblock->get_stateInfo().leftStateInfo;

  assert (cblock->get_leftBlock() == ablock || cblock->get_rightBlock() == ablock);
  if (cblock->get_leftBlock() == ablock)
    {
      //#pragma omp parallel default(shared)  num_threads(num_thrds)
      {
	//#pragma omp for schedule(dynamic)
      for (int lQ = 0; lQ < leftOpSz; ++lQ) {
	for (int lQPrime = 0; lQPrime < leftOpSz; ++lQPrime)
	  {
	    if (a.allowed(lQ, lQPrime))
              {
		const Matrix& aop = a.operator_element(lQ, lQPrime);
		for (int rQ = 0; rQ < rightOpSz; ++rQ)
		  if (c.allowed(lQPrime, rQ) && v.allowed(lQ, rQ))
		    {
                      double fac = scale;
		      fac *= dmrginp.get_ninej()(lS->quanta[lQPrime].get_s().getirrep(), rS->quanta[rQ].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
						   a.get_spin().getirrep(), 0, a.get_spin().getirrep(),
						   lS->quanta[lQ].get_s().getirrep(), rS->quanta[rQ].get_s().getirrep() , v.get_deltaQuantum(0).get_s().getirrep());
		      fac *= Symmetry::spatial_ninej(lS->quanta[lQPrime].get_symm().getirrep() , rS->quanta[rQ].get_symm().getirrep(), c.get_symm().getirrep(), 
					   a.get_symm().getirrep(), 0, a.get_symm().getirrep(),
					   lS->quanta[lQ].get_symm().getirrep() , rS->quanta[rQ].get_symm().getirrep(), v.get_symm().getirrep());
		      fac *= a.get_scaling(lS->quanta[lQ], lS->quanta[lQPrime]);
		      MatrixMultiply (aop, a.conjugacy(), c.operator_element(lQPrime, rQ), c.conjugacy(),
				      v.operator_element(lQ, rQ), fac);
		    }

              }
	  }
      }
      }
    }
  else
    {
      //#pragma omp parallel default(shared)  num_threads(num_thrds)
      {
	//#pragma omp for schedule(dynamic)
      for (int rQ = 0; rQ < rightOpSz; ++rQ) {
	for (int rQPrime = 0; rQPrime < rightOpSz; ++rQPrime)
	  if (a.allowed(rQ, rQPrime))
	    {
	      const Matrix& aop = a.operator_element(rQ, rQPrime);
	      for (int lQ = 0; lQ < leftOpSz; ++lQ) 
		if (v.allowed(lQ, rQ) && c.allowed(lQ, rQPrime)) {
                  double fac = scale;
                  fac *= dmrginp.get_ninej()(lS->quanta[lQ].get_s().getirrep(), rS->quanta[rQPrime].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
					       0, a.get_spin().getirrep(), a.get_spin().getirrep(),
					       lS->quanta[lQ].get_s().getirrep(), rS->quanta[rQ].get_s().getirrep() , v.get_deltaQuantum(0).get_s().getirrep());
		  fac *= Symmetry::spatial_ninej(lS->quanta[lQ].get_symm().getirrep() , rS->quanta[rQPrime].get_symm().getirrep(), c.get_symm().getirrep(), 
				      0, a.get_symm().getirrep(), a.get_symm().getirrep(),
				      lS->quanta[lQ].get_symm().getirrep() , rS->quanta[rQ].get_symm().getirrep(), v.get_symm().getirrep());
		  fac *= a.get_scaling(rS->quanta[rQ], rS->quanta[rQPrime]);
		  double parity = a.get_fermion() && IsFermion(lS->quanta[lQ]) ? -1 : 1;

		  MatrixMultiply (c.operator_element(lQ, rQPrime), c.conjugacy(),
				  aop, TransposeOf(a.conjugacy()), v.operator_element(lQ, rQ), fac*parity);
		}

	    }
      }
      }
    }
}
*/


void SpinAdapted::operatorfunctions::TensorMultiply(const SpinBlock *ablock, const Baseoperator<Matrix>& a, const SpinBlock *cblock, Wavefunction& c, Wavefunction& v, const SpinQuantum dQ, double scale, int num_thrds)
{
  // cannot be used for situation with different bra and ket
  const int leftBraOpSz = cblock->get_leftBlock()->get_braStateInfo().quanta.size ();
  const int leftKetOpSz = cblock->get_leftBlock()->get_ketStateInfo().quanta.size ();
  const int rightBraOpSz = cblock->get_rightBlock()->get_braStateInfo().quanta.size ();
  const int rightKetOpSz = cblock->get_rightBlock()->get_ketStateInfo().quanta.size ();

  const StateInfo* lbraS = cblock->get_braStateInfo().leftStateInfo, *lketS = cblock->get_ketStateInfo().leftStateInfo;
  const StateInfo* rbraS = cblock->get_braStateInfo().rightStateInfo, *rketS = cblock->get_ketStateInfo().rightStateInfo;

  assert (cblock->get_leftBlock() == ablock || cblock->get_rightBlock() == ablock);
  if (cblock->get_leftBlock() == ablock)
    {
      //#pragma omp parallel default(shared)  num_threads(num_thrds)
      {
	//#pragma omp for schedule(dynamic)
      for (int lQ = 0; lQ < leftBraOpSz; ++lQ) {
	for (int lQPrime = 0; lQPrime < leftKetOpSz; ++lQPrime)
	  {
	    if (a.allowed(lQ, lQPrime))
              {
		const Matrix& aop = a.operator_element(lQ, lQPrime);
		  for (int rQ = 0; rQ < rightKetOpSz; ++rQ) 
		    if (c.allowed(lQPrime, rQ) && v.allowed(lQ, rQ))
		    {
                      double fac=scale;
		      fac *= dmrginp.get_ninej()(lketS->quanta[lQPrime].get_s().getirrep(), rketS->quanta[rQ].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
						   a.get_spin().getirrep(), 0, a.get_spin().getirrep(),
						   lbraS->quanta[lQ].get_s().getirrep(), rketS->quanta[rQ].get_s().getirrep() , v.get_deltaQuantum(0).get_s().getirrep());
		      fac *= Symmetry::spatial_ninej(lketS->quanta[lQPrime].get_symm().getirrep() , rketS->quanta[rQ].get_symm().getirrep(), c.get_symm().getirrep(), 
					   a.get_symm().getirrep(), 0, a.get_symm().getirrep(),
					   lbraS->quanta[lQ].get_symm().getirrep() , rketS->quanta[rQ].get_symm().getirrep(), v.get_symm().getirrep());
		      fac *= a.get_scaling(lbraS->quanta[lQ], lketS->quanta[lQPrime]);
		      MatrixMultiply (aop, a.conjugacy(), c.operator_element(lQPrime, rQ), c.conjugacy(),
				      v.operator_element(lQ, rQ), fac);
		    }

              }
	  }
      }
      }
    }
  else
    {
      //#pragma omp parallel default(shared)  num_threads(num_thrds)
      {
	//#pragma omp for schedule(dynamic)
      for (int rQ = 0; rQ < rightBraOpSz; ++rQ) {
	for (int rQPrime = 0; rQPrime < rightKetOpSz; ++rQPrime)
	  if (a.allowed(rQ, rQPrime))
	    {
	      const Matrix& aop = a.operator_element(rQ, rQPrime);
	      for (int lQPrime = 0; lQPrime < leftKetOpSz; ++lQPrime) 
		if (v.allowed(lQPrime, rQ) && c.allowed(lQPrime, rQPrime)) {
                  double fac = scale;
		  fac *= dmrginp.get_ninej()(lketS->quanta[lQPrime].get_s().getirrep(), rketS->quanta[rQPrime].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
					       0, a.get_spin().getirrep(), a.get_spin().getirrep(),
					       lketS->quanta[lQPrime].get_s().getirrep(), rbraS->quanta[rQ].get_s().getirrep() , v.get_deltaQuantum(0).get_s().getirrep());
		  fac *= Symmetry::spatial_ninej(lketS->quanta[lQPrime].get_symm().getirrep() , rketS->quanta[rQPrime].get_symm().getirrep(), c.get_symm().getirrep(), 
				      0, a.get_symm().getirrep(), a.get_symm().getirrep(),
				      lketS->quanta[lQPrime].get_symm().getirrep() , rbraS->quanta[rQ].get_symm().getirrep(), v.get_symm().getirrep());
		  fac *= a.get_scaling(rbraS->quanta[rQ], rketS->quanta[rQPrime]);
		  double parity = a.get_fermion() && IsFermion(lketS->quanta[lQPrime]) ? -1 : 1;

		  MatrixMultiply (c.operator_element(lQPrime, rQPrime), c.conjugacy(),
				  aop, TransposeOf(a.conjugacy()), v.operator_element(lQPrime, rQ), fac*parity);
		}

	    }
      }
      }
    }
}

void SpinAdapted::operatorfunctions::braTensorMultiply(const SpinBlock *ablock, const Baseoperator<Matrix>& a, const SpinBlock *cblock, Wavefunction& c, Wavefunction& v, double scale, int num_thrds)
{
  //It get a result of <\Psi|O. 
  //It is similar to Transposeview(O)|\Psi>
  //However, spin coupling coefficients is different for transpose.
  //It is convenient for npdm with intermediate.
  const int leftBraOpSz = cblock->get_leftBlock()->get_braStateInfo().quanta.size ();
  const int leftKetOpSz = cblock->get_leftBlock()->get_ketStateInfo().quanta.size ();
  const int rightBraOpSz = cblock->get_rightBlock()->get_braStateInfo().quanta.size ();
  const int rightKetOpSz = cblock->get_rightBlock()->get_ketStateInfo().quanta.size ();

  const StateInfo* lbraS = cblock->get_braStateInfo().leftStateInfo, *lketS = cblock->get_ketStateInfo().leftStateInfo;
  const StateInfo* rbraS = cblock->get_braStateInfo().rightStateInfo, *rketS = cblock->get_ketStateInfo().rightStateInfo;

  assert (cblock->get_leftBlock() == ablock || cblock->get_rightBlock() == ablock);
  if (cblock->get_leftBlock() == ablock)
    {
      //#pragma omp parallel default(shared)  num_threads(num_thrds)
      {
	//#pragma omp for schedule(dynamic)
      for (int lQ = 0; lQ < leftKetOpSz; ++lQ) {
	for (int lQPrime = 0; lQPrime < leftBraOpSz; ++lQPrime)
	  {
	    if (a.allowed(lQPrime, lQ))
              {
		const Matrix& aop = a.operator_element(lQPrime, lQ);
		  for (int rQ = 0; rQ < rightBraOpSz; ++rQ) 
		    if (c.allowed(lQPrime, rQ) && v.allowed(lQ, rQ))
		    {
                      double fac=scale;
		      fac *= dmrginp.get_ninej()(lbraS->quanta[lQPrime].get_s().getirrep(), rbraS->quanta[rQ].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
						   (-a.get_spin()).getirrep(), 0, (-a.get_spin()).getirrep(),
						   lketS->quanta[lQ].get_s().getirrep(), rbraS->quanta[rQ].get_s().getirrep() , v.get_deltaQuantum(0).get_s().getirrep());
		      fac *= Symmetry::spatial_ninej(lbraS->quanta[lQPrime].get_symm().getirrep() , rbraS->quanta[rQ].get_symm().getirrep(), c.get_symm().getirrep(), 
					   (-a.get_symm()).getirrep(), 0, (-a.get_symm()).getirrep(),
					   lketS->quanta[lQ].get_symm().getirrep() , rbraS->quanta[rQ].get_symm().getirrep(), v.get_symm().getirrep());
		      fac *= a.get_scaling(lbraS->quanta[lQPrime], lketS->quanta[lQ]);
		      MatrixMultiply (aop, TransposeOf(a.conjugacy()), c.operator_element(lQPrime, rQ), c.conjugacy(),
				      v.operator_element(lQ, rQ), fac);
		    }

              }
	  }
      }
      }
    }
  else
    {
      //#pragma omp parallel default(shared)  num_threads(num_thrds)
      {
	//#pragma omp for schedule(dynamic)
      for (int rQ = 0; rQ < rightKetOpSz; ++rQ) {
	for (int rQPrime = 0; rQPrime < rightBraOpSz; ++rQPrime)
	  if (a.allowed(rQPrime, rQ))
	    {
	      const Matrix& aop = a.operator_element(rQ, rQPrime);
	      for (int lQPrime = 0; lQPrime < leftBraOpSz; ++lQPrime) 
		if (v.allowed(lQPrime, rQ) && c.allowed(lQPrime, rQPrime)) {
                  double fac = scale;
		  fac *= dmrginp.get_ninej()(lbraS->quanta[lQPrime].get_s().getirrep(), rbraS->quanta[rQPrime].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
					       0, (-a.get_spin()).getirrep(), (-a.get_spin()).getirrep(),
					       lbraS->quanta[lQPrime].get_s().getirrep(), rketS->quanta[rQ].get_s().getirrep() , v.get_deltaQuantum(0).get_s().getirrep());
		  fac *= Symmetry::spatial_ninej(lbraS->quanta[lQPrime].get_symm().getirrep() , rbraS->quanta[rQPrime].get_symm().getirrep(), c.get_symm().getirrep(), 
				      0, (-a.get_symm()).getirrep(), (-a.get_symm()).getirrep(),
				      lbraS->quanta[lQPrime].get_symm().getirrep() , rketS->quanta[rQ].get_symm().getirrep(), v.get_symm().getirrep());
		  fac *= a.get_scaling(rbraS->quanta[rQPrime], rketS->quanta[rQ]);
		  double parity = a.get_fermion() && IsFermion(lbraS->quanta[lQPrime]) ? -1 : 1;

		  MatrixMultiply (c.operator_element(lQPrime, rQPrime), c.conjugacy(),
				  aop, a.conjugacy(), v.operator_element(lQPrime, rQ), fac*parity);
		}

	    }
      }
      }
    }
}


void SpinAdapted::operatorfunctions::TensorMultiply(const SpinBlock *ablock, const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, const SpinBlock *cblock, Wavefunction& c, Wavefunction& v, const SpinQuantum opQ, double scale)
{
  // can be used for situation with different bra and ket
  const int leftBraOpSz = cblock->get_leftBlock()->get_braStateInfo().quanta.size ();
  const int leftKetOpSz = cblock->get_leftBlock()->get_ketStateInfo().quanta.size ();
  const int rightBraOpSz = cblock->get_rightBlock()->get_braStateInfo().quanta.size ();
  const int rightKetOpSz = cblock->get_rightBlock()->get_ketStateInfo().quanta.size ();

  const StateInfo* lbraS = cblock->get_braStateInfo().leftStateInfo, *rbraS = cblock->get_braStateInfo().rightStateInfo;
  const StateInfo* lketS = cblock->get_ketStateInfo().leftStateInfo, *rketS = cblock->get_ketStateInfo().rightStateInfo;

  const char conjC = (cblock->get_leftBlock() == ablock) ? 'n' : 't';

  const Baseoperator<Matrix>& leftOp = (conjC == 'n') ? a : b; // an ugly hack to support the release memory optimisation
  const Baseoperator<Matrix>& rightOp = (conjC == 'n') ? b : a;
  const char leftConj = (conjC == 'n') ? a.conjugacy() : b.conjugacy();
  const char rightConj = (conjC == 'n') ? b.conjugacy() : a.conjugacy();


  int totalmem =0;

  for (int lQrQPrime = 0; lQrQPrime<leftBraOpSz*rightKetOpSz; ++lQrQPrime)
  {
    int rQPrime = lQrQPrime%rightKetOpSz, lQ = lQrQPrime/rightKetOpSz;
    for (int lQPrime = 0; lQPrime < leftKetOpSz; lQPrime++)
      if (leftOp.allowed(lQ, lQPrime) && c.allowed(lQPrime, rQPrime))
      {	    
	Matrix m; m.ReSize(lbraS->getquantastates(lQ), rketS->getquantastates(rQPrime));
        
	double factor = leftOp.get_scaling(lbraS->quanta[lQ], lketS->quanta[lQPrime]);
	MatrixMultiply (leftOp.operator_element(lQ, lQPrime), leftConj, c.operator_element(lQPrime, rQPrime), 'n',
			m, factor, 0.);	      
	
	for (int rQ = 0; rQ<rightBraOpSz; rQ++) {
	  if (v.allowed(lQ, rQ) && rightOp.allowed(rQ, rQPrime)) {
	    double factor = scale;
	    
	    factor *= dmrginp.get_ninej()(lketS->quanta[lQPrime].get_s().getirrep(), rketS->quanta[rQPrime].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
					  leftOp.get_spin().getirrep(), rightOp.get_spin().getirrep(), opQ.get_s().getirrep(),
					  lbraS->quanta[lQ].get_s().getirrep(), rbraS->quanta[rQ].get_s().getirrep() , v.get_deltaQuantum(0).get_s().getirrep());
	    factor *= Symmetry::spatial_ninej(lketS->quanta[lQPrime].get_symm().getirrep() , rketS->quanta[rQPrime].get_symm().getirrep(), c.get_symm().getirrep(), 
					      leftOp.get_symm().getirrep(), rightOp.get_symm().getirrep(), opQ.get_symm().getirrep(),
					      lbraS->quanta[lQ].get_symm().getirrep() , rbraS->quanta[rQ].get_symm().getirrep(), v.get_symm().getirrep());
	    int parity = rightOp.get_fermion() && IsFermion(lketS->quanta[lQPrime]) ? -1 : 1;
	    factor *=  rightOp.get_scaling(rbraS->quanta[rQ], rketS->quanta[rQPrime]);
	    MatrixMultiply (m, 'n', rightOp(rQ, rQPrime), TransposeOf(rightOp.conjugacy()), v.operator_element(lQ, rQ), factor*parity);
	  }
	}
	
      }
  }

}

/*
void SpinAdapted::operatorfunctions::TensorMultiply(const SpinBlock *ablock, const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, const SpinBlock *cblock, Wavefunction& c, Wavefunction& v, const SpinQuantum opQ, double scale)
{
  // can be used for situation with different bra and ket
  const int leftBraOpSz = cblock->get_leftBlock()->get_braStateInfo().quanta.size ();
  const int leftKetOpSz = cblock->get_leftBlock()->get_ketStateInfo().quanta.size ();
  const int rightBraOpSz = cblock->get_rightBlock()->get_braStateInfo().quanta.size ();
  const int rightKetOpSz = cblock->get_rightBlock()->get_ketStateInfo().quanta.size ();

  const StateInfo* lbraS = cblock->get_braStateInfo().leftStateInfo, *rbraS = cblock->get_braStateInfo().rightStateInfo;
  const StateInfo* lketS = cblock->get_ketStateInfo().leftStateInfo, *rketS = cblock->get_ketStateInfo().rightStateInfo;

  const char conjC = (cblock->get_leftBlock() == ablock) ? 'n' : 't';

  const Baseoperator<Matrix>& leftOp = (conjC == 'n') ? a : b; // an ugly hack to support the release memory optimisation
  const Baseoperator<Matrix>& rightOp = (conjC == 'n') ? b : a;
  const char leftConj = (conjC == 'n') ? a.conjugacy() : b.conjugacy();
  const char rightConj = (conjC == 'n') ? b.conjugacy() : a.conjugacy();

  Wavefunction u;
  u.resize(leftBraOpSz*leftKetOpSz, rightKetOpSz);

  int totalmem =0;

  {
    for (int lQrQPrime = 0; lQrQPrime<leftBraOpSz*rightKetOpSz; ++lQrQPrime)
    {
      int rQPrime = lQrQPrime%rightKetOpSz, lQ = lQrQPrime/rightKetOpSz;
	for (int lQPrime = 0; lQPrime < leftKetOpSz; lQPrime++)
	  if (leftOp.allowed(lQ, lQPrime) && c.allowed(lQPrime, rQPrime))
	  {
	    int lindex = lQ*leftKetOpSz+lQPrime;
	    u.allowed(lindex, rQPrime) = true;
            
	    u(lindex,rQPrime).ReSize(lbraS->getquantastates(lQ), rketS->getquantastates(rQPrime));
	    double factor = leftOp.get_scaling(lbraS->quanta[lQ], lketS->quanta[lQPrime]);
	    MatrixMultiply (leftOp.operator_element(lQ, lQPrime), leftConj, c.operator_element(lQPrime, rQPrime), 'n',
			    u.operator_element(lindex, rQPrime), factor, 0.);	      

	  }
    }
  }

  pout << "after first step in tensormultiply"<<endl;
      mcheck("before davidson but after all blocks are built");

  {
    for (int lQrQ = 0; lQrQ<leftBraOpSz*rightBraOpSz; ++lQrQ)
    {
      int rQ = lQrQ%rightBraOpSz, lQ=lQrQ/rightBraOpSz;
	if (v.allowed(lQ, rQ))
	  for (int rQPrime = 0; rQPrime < rightKetOpSz; rQPrime++)
	    if (rightOp.allowed(rQ, rQPrime))
	      for (int lQPrime = 0; lQPrime < leftKetOpSz; lQPrime++)
		if (leftOp.allowed(lQ, lQPrime) && u.allowed(lQ*leftKetOpSz+lQPrime, rQPrime))
		{
		  int lindex = lQ*leftKetOpSz+lQPrime;
		  double factor = scale;

		  factor *= dmrginp.get_ninej()(lketS->quanta[lQPrime].get_s().getirrep(), rketS->quanta[rQPrime].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
						leftOp.get_spin().getirrep(), rightOp.get_spin().getirrep(), opQ.get_s().getirrep(),
						lbraS->quanta[lQ].get_s().getirrep(), rbraS->quanta[rQ].get_s().getirrep() , v.get_deltaQuantum(0).get_s().getirrep());
		  factor *= Symmetry::spatial_ninej(lketS->quanta[lQPrime].get_symm().getirrep() , rketS->quanta[rQPrime].get_symm().getirrep(), c.get_symm().getirrep(), 
				       leftOp.get_symm().getirrep(), rightOp.get_symm().getirrep(), opQ.get_symm().getirrep(),
				       lbraS->quanta[lQ].get_symm().getirrep() , rbraS->quanta[rQ].get_symm().getirrep(), v.get_symm().getirrep());
		  int parity = rightOp.get_fermion() && IsFermion(lketS->quanta[lQPrime]) ? -1 : 1;
		  factor *=  rightOp.get_scaling(rbraS->quanta[rQ], rketS->quanta[rQPrime]);
		  MatrixMultiply (u.operator_element(lindex, rQPrime), 'n',
				  rightOp(rQ, rQPrime), TransposeOf(rightOp.conjugacy()), v.operator_element(lQ, rQ), factor*parity);
		}
    }
  }
	      
}
*/
void SpinAdapted::operatorfunctions::OperatorScaleAdd(double scaleV, const SpinBlock& b, const Baseoperator<Matrix>& op1, Baseoperator<Matrix>& op2)
{
  const StateInfo& s = b.get_stateInfo();
  for (int lQ = 0; lQ< op2.nrows(); lQ++)
    for (int rQ = 0; rQ<op2.ncols(); rQ++)
      if (op2.allowed(lQ, rQ) && op1.allowed(lQ,rQ))
      {
	double factor = op1.get_scaling(s.quanta[lQ], s.quanta[rQ]);
	if (op1.conjugacy() == 't')
	  MatrixScaleAdd(scaleV*factor, op1.operator_element(lQ,rQ).t(), op2.operator_element(lQ,rQ));
	else
	  MatrixScaleAdd(scaleV*factor, op1.operator_element(lQ,rQ), op2.operator_element(lQ,rQ));
      }

}

void SpinAdapted::operatorfunctions::MultiplyProduct(const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, Baseoperator<Matrix>& c, Real scale)
{
  if (fabs(scale) < TINY) return;
  const int aSz = a.nrows();
  const int aSzPrime = a.ncols();
  const int bSzPrime = b.ncols();

  assert (a.ncols() == b.nrows() && c.nrows() == a.nrows() &&
          c.ncols() == b.ncols());

  for (int aQ = 0; aQ < aSz; ++aQ)
    for (int aQPrime = 0; aQPrime < aSzPrime; ++aQPrime)
      for (int bQPrime = 0; bQPrime < bSzPrime; ++bQPrime)
        {
          if (a.allowed(aQ, aQPrime) && b.allowed(aQPrime, bQPrime) && c.allowed(aQ, bQPrime) ) {

	    MatrixMultiply (a.operator_element(aQ, aQPrime), a.conjugacy(), b.operator_element(aQPrime, bQPrime), b.conjugacy(),
			    c.operator_element(aQ, bQPrime), scale);
	  }
        }
}


void SpinAdapted::operatorfunctions::TensorTrace (const SpinBlock *ablock, const Baseoperator<Matrix>& a, const SpinBlock* cblock, const StateInfo* cstateinfo, DiagonalMatrix& cDiagonal, Real scale)
{
  if (fabs(scale) < TINY) return;
  try
    {
      assert (a.get_initialised());
      const char conjC = (ablock == cblock->get_leftBlock()) ? 'n' : 't';

      const int aSz = ablock->get_stateInfo().quanta.size ();
      const int bSz = (conjC == 'n') ? cblock->get_stateInfo().rightStateInfo->quanta.size () : cblock->get_stateInfo().leftStateInfo->quanta.size ();
      const StateInfo& s = cblock->get_stateInfo();

      const StateInfo* lS = s.leftStateInfo, *rS = s.rightStateInfo;

      for (int aQ = 0; aQ < aSz; ++aQ)
	if (a.allowed(aQ, aQ))
	  for (int bQ = 0; bQ < bSz; ++bQ)
	    if (s.allowedQuanta (aQ, bQ, conjC))
	      {
		int cQ = s.quantaMap (aQ, bQ, conjC)[0];
		for (int cQState = 0; cQState < s.quantaStates [cQ]; ++cQState)
		  {
		    Real scaleB = 1.0;
		    int aQState;
		    int bQState;
		    
		    if (conjC == 'n')
		      {
			s.UnMapQuantumState (cQState, s.rightStateInfo->quantaStates [bQ], aQState, bQState);
			scaleB *= dmrginp.get_ninej()(lS->quanta[aQ].get_s().getirrep() , rS->quanta[bQ].get_s().getirrep(), cstateinfo->quanta[cQ].get_s().getirrep(), 
				a.get_spin().getirrep(), 0, 0,
				lS->quanta[aQ].get_s().getirrep() , rS->quanta[bQ].get_s().getirrep(), cstateinfo->quanta[cQ].get_s().getirrep());
			scaleB *= Symmetry::spatial_ninej(lS->quanta[aQ].get_symm().getirrep() , rS->quanta[bQ].get_symm().getirrep(), cstateinfo->quanta[cQ].get_symm().getirrep(), 
					   a.get_symm().getirrep(), 0, 0,
					   lS->quanta[aQ].get_symm().getirrep() , rS->quanta[bQ].get_symm().getirrep(), cstateinfo->quanta[cQ].get_symm().getirrep());
		      }
		    else
		      {
			scaleB *= dmrginp.get_ninej()(lS->quanta[bQ].get_s().getirrep() , rS->quanta[aQ].get_s().getirrep(), cstateinfo->quanta[cQ].get_s().getirrep(), 
					0, a.get_spin().getirrep(), 0,
				lS->quanta[bQ].get_s().getirrep() , rS->quanta[aQ].get_s().getirrep(), cstateinfo->quanta[cQ].get_s().getirrep());
			scaleB *= Symmetry::spatial_ninej(lS->quanta[bQ].get_symm().getirrep() , rS->quanta[aQ].get_symm().getirrep(), cstateinfo->quanta[cQ].get_symm().getirrep(), 
					     0, a.get_symm().getirrep(), 0,
					     lS->quanta[bQ].get_symm().getirrep() , rS->quanta[aQ].get_symm().getirrep(), cstateinfo->quanta[cQ].get_symm().getirrep());
			if (a.get_fermion()&& IsFermion(lS->quanta[bQ])) scaleB *= -1.0;
			s.UnMapQuantumState (cQState, s.rightStateInfo->quantaStates [aQ], bQState, aQState);
		      }
		    cDiagonal(s.unBlockedIndex [cQ] + cQState + 1) += scale * scaleB * a.operator_element(aQ, aQ)(aQState + 1, aQState + 1); 

		  }
    }
    }
  catch (Exception)
    {
      pout << Exception::what () << endl;
      abort ();
    }
}

  
void SpinAdapted::operatorfunctions::TensorProduct (const SpinBlock *ablock, const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, const SpinBlock* cblock, const StateInfo* cstateinfo, DiagonalMatrix& cDiagonal, double scale)
{
  if (fabs(scale) < TINY) return;
  const int aSz = a.nrows();
  const int bSz = b.nrows();
  const char conjC = (cblock->get_leftBlock() == ablock) ? 'n' : 't';
  const SpinBlock* bblock = (cblock->get_leftBlock() == ablock) ? cblock->get_rightBlock() : cblock->get_leftBlock();
  const StateInfo& s = cblock->get_stateInfo();
  const StateInfo* lS = s.leftStateInfo, *rS = s.rightStateInfo;

  for (int aQ = 0; aQ < aSz; ++aQ)
    if (a.allowed(aQ, aQ))
      for (int bQ = 0; bQ < bSz; ++bQ)
	if (b.allowed(bQ, bQ))
	  if (s.allowedQuanta (aQ, bQ, conjC))
	  {
	    int cQ = s.quantaMap (aQ, bQ, conjC)[0];
	    Real scaleA = scale;
	    Real scaleB = 1;
	    if (conjC == 'n')
	      {
		scaleB *= dmrginp.get_ninej()(lS->quanta[aQ].get_s().getirrep() , rS->quanta[bQ].get_s().getirrep(), cstateinfo->quanta[cQ].get_s().getirrep(), 
					      a.get_spin().getirrep(), b.get_spin().getirrep(), 0,
					      lS->quanta[aQ].get_s().getirrep() , rS->quanta[bQ].get_s().getirrep(), cstateinfo->quanta[cQ].get_s().getirrep());
		scaleB *= Symmetry::spatial_ninej(lS->quanta[aQ].get_symm().getirrep() , rS->quanta[bQ].get_symm().getirrep(), cstateinfo->quanta[cQ].get_symm().getirrep(), 
				     a.get_symm().getirrep(), b.get_symm().getirrep(), 0,
				     lS->quanta[aQ].get_symm().getirrep() , rS->quanta[bQ].get_symm().getirrep(), cstateinfo->quanta[cQ].get_symm().getirrep());
		
		if (b.get_fermion() && IsFermion (lS->quanta [aQ])) scaleB *= -1.0;
		for (int aQState = 0; aQState < lS->quantaStates[aQ] ; aQState++)
		  MatrixDiagonalScale(a.operator_element(aQ, aQ)(aQState+1, aQState+1)*scaleA*scaleB, b.operator_element(bQ, bQ), 
				      cDiagonal.Store()+s.unBlockedIndex[cQ]+aQState*rS->quantaStates[bQ]);

	      }
	    else
	      {
		scaleB *= dmrginp.get_ninej()(lS->quanta[bQ].get_s().getirrep() , rS->quanta[aQ].get_s().getirrep(), cstateinfo->quanta[cQ].get_s().getirrep(), 
					      b.get_spin().getirrep(), a.get_spin().getirrep(), 0,
					      lS->quanta[bQ].get_s().getirrep() , rS->quanta[aQ].get_s().getirrep(), cstateinfo->quanta[cQ].get_s().getirrep());
		scaleB *= Symmetry::spatial_ninej(lS->quanta[bQ].get_symm().getirrep() , rS->quanta[aQ].get_symm().getirrep(), cstateinfo->quanta[cQ].get_symm().getirrep(), 
				     b.get_symm().getirrep(), a.get_symm().getirrep(), 0,
				     lS->quanta[bQ].get_symm().getirrep() , rS->quanta[aQ].get_symm().getirrep(), cstateinfo->quanta[cQ].get_symm().getirrep());
		
		if (a.get_fermion()&& IsFermion(lS->quanta[bQ])) scaleB *= -1.0;
		for (int bQState = 0; bQState < lS->quantaStates[bQ] ; bQState++)
		  MatrixDiagonalScale(b.operator_element(bQ, bQ)(bQState+1, bQState+1)*scaleA*scaleB, a.operator_element(aQ, aQ), 
				      cDiagonal.Store()+s.unBlockedIndex[cQ]+bQState*rS->quantaStates[aQ]);
	      }
	  }
}


void SpinAdapted::operatorfunctions::TensorProduct (const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, const StateInfo *brastateinfo, const StateInfo *ketstateinfo, Baseoperator<Matrix>& c, double scale, bool aIsLeftOp, int num_thrds)
{
  if (fabs(scale) < TINY) return;
  int rows = c.nrows();
  int cols = c.ncols();
#ifdef _OPENMP
#pragma omp parallel default(shared) num_threads(num_thrds) 
#endif
  {
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
  for (int cq = 0; cq < rows; ++cq)
    for (int cqprime = 0; cqprime < cols; ++cqprime)
      if (c.allowed(cq, cqprime)) {
	TensorProductElement(a, b, brastateinfo, ketstateinfo, c, c.operator_element(cq, cqprime), cq, cqprime, aIsLeftOp, scale);
      }
  }
}


void SpinAdapted::operatorfunctions::TensorProductElement(const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, const StateInfo *brastateinfo, const StateInfo *ketstateinfo, Baseoperator<Matrix>& c, Matrix& cel, int cq, int cqprime, bool aIsLeftOp, double scale)
{
  if (fabs(scale) < TINY) return;
  assert (a.get_initialised());
  assert (b.get_initialised());
  assert (c.get_initialised());
 
  const std::vector<int>& oldToNewI = brastateinfo->oldToNewState.at(cq);
  const std::vector<int>& oldToNewJ = ketstateinfo->oldToNewState.at(cqprime);
 
  const char conjC = aIsLeftOp ? 'n' : 't';
  
  const StateInfo* lbraS = brastateinfo->leftStateInfo, *rbraS = brastateinfo->rightStateInfo;
  const StateInfo* lketS = ketstateinfo->leftStateInfo, *rketS = ketstateinfo->rightStateInfo;
  int rowstride = 0, colstride = 0;

  int aq, aqprime, bq, bqprime;

  //pout << "old to new size "<<oldToNewI.size()<<" "<<oldToNewJ.size()<<endl;
  for (int oldi =0; oldi < oldToNewI.size(); oldi++) {
    colstride = 0;
    for (int oldj = 0; oldj < oldToNewJ.size(); oldj++)
    {
      if (conjC == 'n')
      {
	aq = brastateinfo->leftUnMapQuanta[ oldToNewI[oldi] ];
	aqprime = ketstateinfo->leftUnMapQuanta[ oldToNewJ[oldj] ];
	bq = brastateinfo->rightUnMapQuanta[ oldToNewI[oldi] ];
	bqprime = ketstateinfo->rightUnMapQuanta[ oldToNewJ[oldj] ];
      }
      else 
      {
	aq = brastateinfo->rightUnMapQuanta[ oldToNewI[oldi] ];
	aqprime = ketstateinfo->rightUnMapQuanta[ oldToNewJ[oldj] ];
	bq = brastateinfo->leftUnMapQuanta[ oldToNewI[oldi] ];
	bqprime = ketstateinfo->leftUnMapQuanta[ oldToNewJ[oldj] ];
      }
  
      Real scaleA = scale;
      Real scaleB = 1.0;
      if (a.allowed(aq, aqprime) && b.allowed(bq, bqprime))
      {
	if (conjC == 'n')
	{
	  scaleB = dmrginp.get_ninej()(lketS->quanta[aqprime].get_s().getirrep() , rketS->quanta[bqprime].get_s().getirrep(), ketstateinfo->quanta[cqprime].get_s().getirrep(), 
			 a.get_spin().getirrep(), b.get_spin().getirrep(), c.get_spin().getirrep(),
			 lbraS->quanta[aq].get_s().getirrep() , rbraS->quanta[bq].get_s().getirrep(), brastateinfo->quanta[cq].get_s().getirrep());
	  scaleB *= Symmetry::spatial_ninej(lketS->quanta[aqprime].get_symm().getirrep() , rketS->quanta[bqprime].get_symm().getirrep(), ketstateinfo->quanta[cqprime].get_symm().getirrep(), 
			       a.get_symm().getirrep(), b.get_symm().getirrep(), c.get_symm().getirrep(),
			       lbraS->quanta[aq].get_symm().getirrep() , rbraS->quanta[bq].get_symm().getirrep(), brastateinfo->quanta[cq].get_symm().getirrep());
	  scaleB *= b.get_scaling(rbraS->quanta[bq], rketS->quanta[bqprime]);
	  scaleA *= a.get_scaling(lbraS->quanta[aq], lketS->quanta[aqprime]);
	  if (b.get_fermion() && IsFermion (ketstateinfo->leftStateInfo->quanta [aqprime])) scaleB *= -1;
	  MatrixTensorProduct (a.operator_element(aq, aqprime), a.conjugacy(), scaleA, 
			       b.operator_element(bq, bqprime), b.conjugacy(), scaleB, cel,rowstride, colstride);
	}
	else
	{
	  scaleB = dmrginp.get_ninej()(lketS->quanta[bqprime].get_s().getirrep(), rketS->quanta[aqprime].get_s().getirrep() , ketstateinfo->quanta[cqprime].get_s().getirrep(), 
			 b.get_spin().getirrep(), a.get_spin().getirrep(), c.get_spin().getirrep(),
			 lbraS->quanta[bq].get_s().getirrep(), rbraS->quanta[aq].get_s().getirrep() , brastateinfo->quanta[cq].get_s().getirrep());
	  scaleB *= Symmetry::spatial_ninej(lketS->quanta[bqprime].get_symm().getirrep() , rketS->quanta[aqprime].get_symm().getirrep(), ketstateinfo->quanta[cqprime].get_symm().getirrep(), 
			       b.get_symm().getirrep(), a.get_symm().getirrep(), c.get_symm().getirrep(),
			       lbraS->quanta[bq].get_symm().getirrep() , rbraS->quanta[aq].get_symm().getirrep(), brastateinfo->quanta[cq].get_symm().getirrep());
	  scaleB *= b.get_scaling(lbraS->quanta[bq], lketS->quanta[bqprime]);
	  scaleA *= a.get_scaling(rbraS->quanta[aq], rketS->quanta[aqprime]);
	  if (a.get_fermion() && IsFermion (ketstateinfo->leftStateInfo->quanta[bqprime]) ) scaleB *= -1.;
	  
	  MatrixTensorProduct (b.operator_element(bq, bqprime), b.conjugacy(), scaleB, 
			       a.operator_element(aq, aqprime), a.conjugacy(), scaleA, cel, rowstride, colstride);
	}
      }
      colstride += ketstateinfo->unCollectedStateInfo->quantaStates[ oldToNewJ[oldj] ];

    }
    rowstride += brastateinfo->unCollectedStateInfo->quantaStates[ oldToNewI[oldi] ];
    
  }

  
}


void SpinAdapted::operatorfunctions::TensorMultiply(const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, const StateInfo *brastateinfo, const StateInfo *ketstateinfo, const Wavefunction& c, Wavefunction& v, const SpinQuantum opQ, bool aIsLeftOp, double scale)
{
  const int leftBraOpSz = brastateinfo->leftStateInfo->quanta.size ();
  const int leftKetOpSz = ketstateinfo->leftStateInfo->quanta.size ();
  const int rightBraOpSz = brastateinfo->rightStateInfo->quanta.size ();
  const int rightKetOpSz = ketstateinfo->rightStateInfo->quanta.size ();

  const StateInfo* lbraS = brastateinfo->leftStateInfo, *rbraS = brastateinfo->rightStateInfo;
  const StateInfo* lketS = ketstateinfo->leftStateInfo, *rketS = ketstateinfo->rightStateInfo;

  const char conjC = (aIsLeftOp) ? 'n' : 't';

  const Baseoperator<Matrix>& leftOp = (conjC == 'n') ? a : b; // an ugly hack to support the release memory optimisation
  const Baseoperator<Matrix>& rightOp = (conjC == 'n') ? b : a;
  const char leftConj = (conjC == 'n') ? a.conjugacy() : b.conjugacy();
  const char rightConj = (conjC == 'n') ? b.conjugacy() : a.conjugacy();

  Wavefunction u;
  u.resize(leftBraOpSz*leftKetOpSz, rightKetOpSz);

  int totalmem =0;

  {
    for (int lQrQPrime = 0; lQrQPrime<leftBraOpSz*rightKetOpSz; ++lQrQPrime)
    {
      int rQPrime = lQrQPrime%rightKetOpSz, lQ = lQrQPrime/rightKetOpSz;
	for (int lQPrime = 0; lQPrime < leftKetOpSz; lQPrime++)
	  if (leftOp.allowed(lQ, lQPrime) && c.allowed(lQPrime, rQPrime))
	  {
	    int lindex = lQ*leftKetOpSz+lQPrime;
	    u.allowed(lindex, rQPrime) = true;

	    u(lindex,rQPrime).ReSize(lbraS->getquantastates(lQ), rketS->getquantastates(rQPrime));
	    double factor = leftOp.get_scaling(lbraS->quanta[lQ], lketS->quanta[lQPrime]);
	    MatrixMultiply (leftOp.operator_element(lQ, lQPrime), leftConj, c.operator_element(lQPrime, rQPrime), 'n',
			    u.operator_element(lindex, rQPrime), factor, 0.);	      

	  }
    }
  }

  {
    for (int lQrQ = 0; lQrQ<leftBraOpSz*rightBraOpSz; ++lQrQ)
    {
      int rQ = lQrQ%rightBraOpSz, lQ=lQrQ/rightBraOpSz;
	if (v.allowed(lQ, rQ))
	  for (int rQPrime = 0; rQPrime < rightKetOpSz; rQPrime++)
	    if (rightOp.allowed(rQ, rQPrime))
	      for (int lQPrime = 0; lQPrime < leftKetOpSz; lQPrime++)
		if (leftOp.allowed(lQ, lQPrime) && u.allowed(lQ*leftKetOpSz+lQPrime, rQPrime))
		{
		  int lindex = lQ*leftKetOpSz+lQPrime;
		  double factor = scale;
      //if(dmrginp.spinAdapted()){
      //ninej has already considered non spin-adapted
      //it is just 1 in nonspin-adapted

		  factor *= dmrginp.get_ninej()(lketS->quanta[lQPrime].get_s().getirrep(), rketS->quanta[rQPrime].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
						leftOp.get_spin().getirrep(), rightOp.get_spin().getirrep(), opQ.get_s().getirrep(),
						lbraS->quanta[lQ].get_s().getirrep(), rbraS->quanta[rQ].get_s().getirrep() , v.get_deltaQuantum(0).get_s().getirrep());
      //}
		  factor *= Symmetry::spatial_ninej(lketS->quanta[lQPrime].get_symm().getirrep() , rketS->quanta[rQPrime].get_symm().getirrep(), c.get_symm().getirrep(), 
				       leftOp.get_symm().getirrep(), rightOp.get_symm().getirrep(), opQ.get_symm().getirrep(),
				       lbraS->quanta[lQ].get_symm().getirrep() , rbraS->quanta[rQ].get_symm().getirrep(), v.get_symm().getirrep());
		  int parity = rightOp.get_fermion() && IsFermion(lketS->quanta[lQPrime]) ? -1 : 1;
		  factor *=  rightOp.get_scaling(rbraS->quanta[rQ], rketS->quanta[rQPrime]);
		  MatrixMultiply (u.operator_element(lindex, rQPrime), 'n',
				  rightOp(rQ, rQPrime), TransposeOf(rightOp.conjugacy()), v.operator_element(lQ, rQ), factor*parity);
		}
    }
  }
	      
}

void SpinAdapted::operatorfunctions::TensorMultiply(const Baseoperator<Matrix>& a, const StateInfo *brastateinfo, const StateInfo *ketstateinfo, const Wavefunction& c, Wavefunction& v, const SpinQuantum dQ, bool left, double scale,  int num_thrds)
{
  //Calculate O_{l or r} |\Psi> without building big block.
  const StateInfo* lbraS = brastateinfo->leftStateInfo, *lketS = ketstateinfo->leftStateInfo;
  const StateInfo* rbraS = brastateinfo->rightStateInfo, *rketS = ketstateinfo->rightStateInfo;
  const int leftBraOpSz = brastateinfo->leftStateInfo->quanta.size ();
  const int leftKetOpSz = ketstateinfo->leftStateInfo->quanta.size ();
  const int rightBraOpSz = brastateinfo->rightStateInfo->quanta.size ();
  const int rightKetOpSz = ketstateinfo->rightStateInfo->quanta.size ();

  if (left)
    {
      //#pragma omp parallel default(shared)  num_threads(num_thrds)
      {
	//#pragma omp for schedule(dynamic)
      for (int lQ = 0; lQ < leftBraOpSz; ++lQ) {
	for (int lQPrime = 0; lQPrime < leftKetOpSz; ++lQPrime)
	  {
	    if (a.allowed(lQ, lQPrime))
              {
		const Matrix& aop = a.operator_element(lQ, lQPrime);
		  for (int rQ = 0; rQ < rightKetOpSz; ++rQ) 
		    if (c.allowed(lQPrime, rQ) && v.allowed(lQ, rQ))
		    {
                      double fac=scale;
		      fac *= dmrginp.get_ninej()(lketS->quanta[lQPrime].get_s().getirrep(), rketS->quanta[rQ].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
						   a.get_spin().getirrep(), 0, a.get_spin().getirrep(),
						   lbraS->quanta[lQ].get_s().getirrep(), rbraS->quanta[rQ].get_s().getirrep() , v.get_deltaQuantum(0).get_s().getirrep());
		      fac *= Symmetry::spatial_ninej(lketS->quanta[lQPrime].get_symm().getirrep() , rketS->quanta[rQ].get_symm().getirrep(), c.get_symm().getirrep(), 
					   a.get_symm().getirrep(), 0, a.get_symm().getirrep(),
					   lbraS->quanta[lQ].get_symm().getirrep() , rbraS->quanta[rQ].get_symm().getirrep(), v.get_symm().getirrep());
		      fac *= a.get_scaling(lbraS->quanta[lQ], lketS->quanta[lQPrime]);
		      MatrixMultiply (aop, a.conjugacy(), c.operator_element(lQPrime, rQ), c.conjugacy(),
				      v.operator_element(lQ, rQ), fac);
		    }

              }
	  }
      }
      }
    }
  else
    {
      //#pragma omp parallel default(shared)  num_threads(num_thrds)
      {
	//#pragma omp for schedule(dynamic)
      for (int rQ = 0; rQ < rightBraOpSz; ++rQ) {
	for (int rQPrime = 0; rQPrime < rightKetOpSz; ++rQPrime)
	  if (a.allowed(rQ, rQPrime))
	    {
	      const Matrix& aop = a.operator_element(rQ, rQPrime);
	      for (int lQPrime = 0; lQPrime < leftKetOpSz; ++lQPrime) 
		if (v.allowed(lQPrime, rQ) && c.allowed(lQPrime, rQPrime)) {
                  double fac = scale;
		  fac *= dmrginp.get_ninej()(lketS->quanta[lQPrime].get_s().getirrep(), rketS->quanta[rQPrime].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
					       0, a.get_spin().getirrep(), a.get_spin().getirrep(),
					       lbraS->quanta[lQPrime].get_s().getirrep(), rbraS->quanta[rQ].get_s().getirrep() , v.get_deltaQuantum(0).get_s().getirrep());
		  fac *= Symmetry::spatial_ninej(lketS->quanta[lQPrime].get_symm().getirrep() , rketS->quanta[rQPrime].get_symm().getirrep(), c.get_symm().getirrep(), 
				      0, a.get_symm().getirrep(), a.get_symm().getirrep(),
				      lbraS->quanta[lQPrime].get_symm().getirrep() , rbraS->quanta[rQ].get_symm().getirrep(), v.get_symm().getirrep());
		  fac *= a.get_scaling(rbraS->quanta[rQ], rketS->quanta[rQPrime]);
		  double parity = a.get_fermion() && IsFermion(lketS->quanta[lQPrime]) ? -1 : 1;

		  MatrixMultiply (c.operator_element(lQPrime, rQPrime), c.conjugacy(),
				  aop, TransposeOf(a.conjugacy()), v.operator_element(lQPrime, rQ), fac*parity);
		}

	    }
      }
      }
    }
}


