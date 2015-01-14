/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "twopdm.h"
#include "MatrixBLAS.h"
#include <boost/format.hpp>
#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include "operatorfunctions.h"
#include "execinfo.h"
#include "pario.h"

namespace SpinAdapted{
  void compute_twopdm_sweep(std::vector<Wavefunction>& wavefunctions, const SpinBlock& system, const SpinBlock& systemDot, const SpinBlock& newSystem, const SpinBlock& newEnvironment, const SpinBlock& big, const int numprocs, int state)
{
  const int nroots = wavefunctions.size();
  array_4d<double> twopdm(big.size()*2, big.size()*2, big.size()*2, big.size()*2);
    //for(int j=0;j<=i;++j)
      {
	int i = state;
	int j = i;
	Wavefunction &wavefunction1 = wavefunctions[0];
	Wavefunction &wavefunction2 = wavefunctions[0];
	load_twopdm_binary(twopdm, i ,j);

	const std::vector<int> distribute_work = distribute_procs(numprocs,4);

	p2out << "\t\t\t compute 1_3_0"<<endl;
	if(mpigetrank() == distribute_work[0])
	  compute_two_pdm_1_3_0(wavefunction1, wavefunction2, big, twopdm);

	p2out << "\t\t\t compute 0_4_0"<<endl;
	if(mpigetrank() == distribute_work[1])
	  compute_two_pdm_0_4_0(wavefunction1, wavefunction2, big, twopdm);

	p2out << "\t\t\t compute 0_3_1"<<endl;
	compute_two_pdm_0_3_1(wavefunction1, wavefunction2, big, twopdm);

	p2out << "\t\t\t compute 1_2_1"<<endl;
	compute_two_pdm_1_2_1(wavefunction1, wavefunction2, big, twopdm);

	p2out << "\t\t\t compute 0_2_2"<<endl;
	compute_two_pdm_0_2_2(wavefunction1, wavefunction2, big, twopdm);

	p2out << "\t\t\t compute 1_1_2"<<endl;
	compute_two_pdm_1_1_2(wavefunction1, wavefunction2, big, twopdm);

	accumulate_twopdm(twopdm);
	save_twopdm_binary(twopdm, i ,j);
      }
}

void compute_twopdm_final(std::vector<Wavefunction>& wavefunctions, const SpinBlock& system, const SpinBlock& systemDot, const SpinBlock& newSystem, const SpinBlock& newEnvironment, const SpinBlock& big, const int numprocs, int state)
{
  const int nroots = wavefunctions.size();
  array_4d<double> twopdm(big.size()*2, big.size()*2, big.size()*2, big.size()*2);
  //for(int i=0;i<nroots;++i)
    //for(int j=0;j<=i;++j)
      {
	int i = state;
	int j = i;
	Wavefunction &wavefunction1 = wavefunctions[0];
	Wavefunction &wavefunction2 = wavefunctions[0];
	load_twopdm_binary(twopdm, i ,j);
	pout <<"\t\t\t Performing sweep calculation "<<endl;
	const std::vector<int> distribute_work = distribute_procs(numprocs,2);

	p2out << "\t\t\t compute 0_0_4"<<endl;	
	if(mpigetrank() == distribute_work[0])
	  compute_two_pdm_0_0_4(wavefunction1, wavefunction2, big, twopdm);

	p2out << "\t\t\t compute 1_3"<<endl;
	if(mpigetrank() == distribute_work[1])
	  compute_two_pdm_1_3(wavefunction1, wavefunction2, big, twopdm);

	accumulate_twopdm(twopdm);
	save_twopdm_binary(twopdm, i ,j);
      }
}


void compute_twopdm_initial(std::vector<Wavefunction>& wavefunctions, const SpinBlock& system, const SpinBlock& systemDot, const SpinBlock& newSystem, const SpinBlock& newEnvironment, const SpinBlock& big, const int numprocs, int state)
{
  const int nroots = wavefunctions.size();
  array_4d<double> twopdm(big.size()*2, big.size()*2, big.size()*2, big.size()*2);
  //for(int i=0;i<nroots;++i)
    //for(int j=0;j<=i;++j)
      {
	int i = state;
	int j = i;
	Wavefunction &wavefunction1 = wavefunctions[0];//there is only one wavefunction in the vector the 
	Wavefunction &wavefunction2 = wavefunctions[0];
	load_twopdm_binary(twopdm, i ,j);
	const std::vector<int> distribute_work = distribute_procs(numprocs,3);

	pout <<"\t\t\t Performing sweep calculation: 2PDM "<<endl;

	p2out << "\t\t\t compute 4_0_0"<<endl;	
	if(mpigetrank() == distribute_work[0])
	  compute_two_pdm_4_0_0(wavefunction1, wavefunction2, big, twopdm);

	p2out << "\t\t\t compute 3_1_0"<<endl;
	if(mpigetrank() == distribute_work[1])
	  compute_two_pdm_3_1_0(wavefunction1, wavefunction2, big, twopdm);

	p2out << "\t\t\t compute 2_2_0"<<endl;
	if(mpigetrank() == distribute_work[2])
	  compute_two_pdm_2_2_0(wavefunction1, wavefunction2, big, twopdm);

	p2out << "\t\t\t compute 3_0_1"<<endl;
	compute_two_pdm_3_0_1(wavefunction1, wavefunction2, big, twopdm);

	p2out << "\t\t\t compute 2_0_2"<<endl;
	compute_two_pdm_2_0_2(wavefunction1, wavefunction2, big, twopdm);

	p2out << "\t\t\t compute 2_1_1"<<endl;
	compute_two_pdm_2_1_1(wavefunction1, wavefunction2, big, twopdm);

	accumulate_twopdm(twopdm);
	save_twopdm_binary(twopdm, i ,j);
      }
}


void compute_two_pdm_0_4_0(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm)
{
  SpinBlock* leftBlock = big.get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();
  SpinBlock* dotBlock = leftBlock->get_rightBlock();
  bool start = leftBlock->get_sites().size() ==2 ? true : false;
  int dotindex = dotBlock->get_sites()[0];

  SpinQuantum sq0 = (getSpinQuantum(dotindex)+getSpinQuantum(dotindex))[0];//(SpinQuantum(1,1,SymmetryOfSpatialOrb(dotindex))+SpinQuantum(1,1,SymmetryOfSpatialOrb(dotindex)))[0];
  boost::shared_ptr<SparseMatrix> dotop0 = dotBlock->get_op_rep(CRE_CRE, sq0, dotindex, dotindex);
  Cre dotop;
  dotop.set_orbs() = dotop0->get_orbs(); dotop.set_orbs().push_back(dotindex); dotop.set_orbs().push_back(dotindex);
  dotop.set_initialised() = true;
  dotop.set_fermion() = false;
  dotop.set_deltaQuantum(1, (dotop0->get_deltaQuantum(0) - dotop0->get_deltaQuantum(0))[0]);
  dotop.allocate(dotBlock->get_stateInfo());
  operatorfunctions::Product(dotBlock, *dotop0, Transposeview(*dotop0), dotop, 1.0);

  SparseMatrix* leftOp = 0;
  SparseMatrix* rightOp = 0;
  vector<double> expectations;
  spinExpectation(wave1, wave2, *leftOp, dotop, *rightOp, big, expectations, false);
  int ix = 2*dotindex;
  assign_antisymmetric(twopdm, ix, ix+1, ix+1, ix, 0.5*expectations[0]);
}

void compute_two_pdm_4_0_0(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm)
{
  SpinBlock* leftBlock = big.get_leftBlock()->get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();
  SpinBlock* dotBlock = leftBlock->get_rightBlock();
  int leftindex = leftBlock->get_sites()[0];

  SpinQuantum sq0 = (getSpinQuantum(leftindex)+getSpinQuantum(leftindex))[0];//(SpinQuantum(1,1,SymmetryOfSpatialOrb(leftindex))+SpinQuantum(1,1,SymmetryOfSpatialOrb(leftindex)))[0];
  boost::shared_ptr<SparseMatrix> leftop0 = leftBlock->get_op_rep(CRE_CRE, sq0, leftindex, leftindex);

  Cre leftop;
  leftop.set_orbs() = leftop0->get_orbs(); leftop.set_orbs().push_back(leftindex); leftop.set_orbs().push_back(leftindex);
  leftop.set_initialised() = true;
  leftop.set_fermion() = false;
  leftop.set_deltaQuantum(1, (leftop0->get_deltaQuantum(0) - leftop0->get_deltaQuantum(0))[0]);
  leftop.allocate(leftBlock->get_stateInfo());
  operatorfunctions::Product(leftBlock, *leftop0, Transposeview(*leftop0), leftop, 1.0);

  SparseMatrix* dotOp = 0;
  SparseMatrix* rightOp = 0;
  vector<double> expectations;
  spinExpectation(wave1, wave2, leftop, *dotOp, *rightOp, big, expectations, false);
  int ix = 2*leftindex;
  assign_antisymmetric(twopdm, ix, ix+1, ix+1, ix, 0.5*expectations[0]);
}

void compute_two_pdm_0_0_4(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm)
{
  //only on the final sweep iteration
  SpinBlock* leftBlock = big.get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();
  SpinBlock* dotBlock = leftBlock->get_rightBlock();
  bool start = leftBlock->get_sites().size() ==2 ? true : false;
  int rightindex = rightBlock->get_sites()[0];

  SpinQuantum sq0 = (getSpinQuantum(rightindex)+getSpinQuantum(rightindex))[0];//(SpinQuantum(1,1,SymmetryOfSpatialOrb(rightindex))+SpinQuantum(1,1,SymmetryOfSpatialOrb(rightindex)))[0];
  boost::shared_ptr<SparseMatrix> rightop0 = rightBlock->get_op_rep(CRE_CRE, sq0, rightindex, rightindex);
  Cre rightop;
  rightop.set_orbs() = rightop0->get_orbs(); rightop.set_orbs().push_back(rightindex); rightop.set_orbs().push_back(rightindex);
  rightop.set_initialised() = true;
  rightop.set_fermion() = false;
  rightop.set_deltaQuantum(1, (rightop0->get_deltaQuantum(0) - rightop0->get_deltaQuantum(0))[0]);
  rightop.allocate(rightBlock->get_stateInfo());
  operatorfunctions::Product(rightBlock, *rightop0, Transposeview(*rightop0), rightop, 1.0);

  SparseMatrix* leftOp = 0;
  SparseMatrix* dotop = 0;
  vector<double> expectations;
  spinExpectation(wave1, wave2, *leftOp, *dotop, rightop, big, expectations, false);
  int ix = 2*rightindex;
  assign_antisymmetric(twopdm, ix, ix+1, ix+1, ix, 0.5*expectations[0]);
}

void compute_two_pdm_0_2_2(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm)
{
  SpinBlock* leftBlock = big.get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();
  SpinBlock* dotBlock = leftBlock->get_rightBlock();
  bool start = leftBlock->get_sites().size() ==2 ? true : false;
  int dotindex = dotBlock->get_sites()[0];

  for (int ij = 0; ij < dotBlock->get_op_array(CRE_CRE).get_size(); ++ij)
  {
    boost::shared_ptr<SparseMatrix> dotop2 = dotBlock->get_op_array(CRE_CRE).get_local_element(ij)[1]->getworkingrepresentation(leftBlock);
    int ix = dotop2->get_orbs(0);
    int jx = dotop2->get_orbs(1);
    SpinQuantum sq0 = (getSpinQuantum(ix)+getSpinQuantum(jx))[0];//(SpinQuantum(1,1,SymmetryOfSpatialOrb(ix))+SpinQuantum(1,1,SymmetryOfSpatialOrb(jx)))[0];
    boost::shared_ptr<SparseMatrix> dotop0 = dotBlock->get_op_rep(CRE_CRE, sq0, ix, jx);
    SparseMatrix* leftOp = 0;
        
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
    {
      int ccsize = rightBlock->get_op_array(CRE_CRE).get_size();
#ifdef _OPENMP
#pragma omp for schedule(guided) nowait
#endif
    for (int kl =0; kl <ccsize; kl++)
    {
      boost::shared_ptr<SparseMatrix> rightop2 = rightBlock->get_op_array(CRE_CRE).get_local_element(kl)[1]->getworkingrepresentation(rightBlock);
      int kx = rightop2->get_orbs(0);
      int lx = rightop2->get_orbs(1);

      SpinQuantum sq0 = (getSpinQuantum(kx)+getSpinQuantum(lx))[0];//(SpinQuantum(1,1,SymmetryOfSpatialOrb(kx))+SpinQuantum(1,1,SymmetryOfSpatialOrb(lx)))[0];
      boost::shared_ptr<SparseMatrix> rightop0 = rightBlock->get_op_rep(CRE_CRE, sq0, kx, lx);
      vector<double> expectations;
      Transposeview rop0 = Transposeview(*rightop0), rop2 = Transposeview(*rightop2); 
      spinExpectation(wave1, wave2, *leftOp, *dotop0, rop0, big, expectations, false);
      spinExpectation(wave1, wave2, *leftOp, *dotop2, rop2, big, expectations, false);
      vector<int> indices(4,0);
      indices[0] = ix; indices[1] = jx; indices[2] = lx; indices[3] = kx;
      expectations[0]*=-1; expectations[1]*=-1;
      spin_to_nonspin(indices, expectations, twopdm, CC_DD, true);
    }
    }
  }      


  for (int ij = 0; ij < dotBlock->get_op_array(CRE_DES).get_size(); ++ij)
  {
    boost::shared_ptr<SparseMatrix> dotop2 = dotBlock->get_op_array(CRE_DES).get_local_element(ij)[1]->getworkingrepresentation(leftBlock);
    int ix = dotop2->get_orbs(0);
    int jx = dotop2->get_orbs(1);
    boost::shared_ptr<SparseMatrix> dotop0 = dotBlock->get_op_array(CRE_DES).get_local_element(ij)[0]->getworkingrepresentation(leftBlock);//dotBlock->get_op_rep(CRE_DES_S0, ix, jx);
    SparseMatrix* leftOp = 0;
        
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
    {
      int cdsize = rightBlock->get_op_array(CRE_DES).get_size();
#ifdef _OPENMP
#pragma omp for schedule(guided) nowait
#endif
    for (int kl =0; kl <cdsize; kl++)
    {
      boost::shared_ptr<SparseMatrix> rightop2 = rightBlock->get_op_array(CRE_DES).get_local_element(kl)[1]->getworkingrepresentation(rightBlock);
      int kx = rightop2->get_orbs(0);
      int lx = rightop2->get_orbs(1);
      boost::shared_ptr<SparseMatrix> rightop0 = rightBlock->get_op_array(CRE_DES).get_local_element(kl)[0]->getworkingrepresentation(rightBlock);//rightBlock->get_op_rep(CRE_DES_S0, kx, lx);
      vector<double> expectations;
      spinExpectation(wave1, wave2, *leftOp, *dotop0, *rightop0, big, expectations, false);
      spinExpectation(wave1, wave2, *leftOp, *dotop2, *rightop2, big, expectations, false);
    
      vector<int> indices(4,0);
      indices[0] = ix; indices[1] = kx; indices[2] = jx; indices[3] = lx;
      spin_to_nonspin(indices, expectations, twopdm, CD_CD, true);
    }
    }
  }      
}

void compute_two_pdm_2_0_2(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm)
{
  SpinBlock* leftBlock = big.get_leftBlock()->get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();
  SpinBlock* dotBlock = leftBlock->get_rightBlock();
  bool start = leftBlock->get_sites().size() ==2 ? true : false;
  int leftindex = leftBlock->get_sites()[0];

#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
  {
  int ccsize = leftBlock->get_op_array(CRE_CRE).get_size();
#ifdef _OPENMP
#pragma omp for schedule(guided) nowait
#endif
  for (int ij = 0; ij < ccsize; ++ij)
  {
    boost::shared_ptr<SparseMatrix> leftop2 = leftBlock->get_op_array(CRE_CRE).get_local_element(ij)[1]->getworkingrepresentation(leftBlock);
    int ix = leftop2->get_orbs(0);
    int jx = leftop2->get_orbs(1);
    boost::shared_ptr<SparseMatrix> leftop0 = leftBlock->get_op_array(CRE_CRE).get_local_element(ij)[0]->getworkingrepresentation(leftBlock);//leftBlock->get_op_rep(CRE_CRE_S0, ix, jx);
    SparseMatrix* dotop = 0;
        
    for (int kl =0; kl <rightBlock->get_op_array(CRE_CRE).get_size(); kl++)
    {
      boost::shared_ptr<SparseMatrix> rightop2 = rightBlock->get_op_array(CRE_CRE).get_local_element(kl)[1]->getworkingrepresentation(rightBlock);
      int kx = rightop2->get_orbs(0);
      int lx = rightop2->get_orbs(1);
      boost::shared_ptr<SparseMatrix> rightop0 = rightBlock->get_op_array(CRE_CRE).get_local_element(kl)[0]->getworkingrepresentation(rightBlock);
      vector<double> expectations;
      Transposeview rop0 = Transposeview(*rightop0), rop2 = Transposeview(*rightop2); 
      spinExpectation(wave1, wave2, *leftop0, *dotop, rop0, big, expectations, false);
      spinExpectation(wave1, wave2, *leftop2, *dotop, rop2, big, expectations, false);
      vector<int> indices(4,0);
      indices[0] = ix; indices[1] = jx; indices[2] = lx; indices[3] = kx;
      expectations[0]*=-1; expectations[1]*=-1;
      spin_to_nonspin(indices, expectations, twopdm, CC_DD, true);
    }
  }      
  }

#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
  {
    int cdsize = leftBlock->get_op_array(CRE_DES).get_size();
#ifdef _OPENMP
#pragma omp for schedule(guided) nowait
#endif
  for (int ij = 0; ij < cdsize; ++ij)
  {
    boost::shared_ptr<SparseMatrix> leftop2 = leftBlock->get_op_array(CRE_DES).get_local_element(ij)[1]->getworkingrepresentation(leftBlock);
    int ix = leftop2->get_orbs(0);
    int jx = leftop2->get_orbs(1);
    boost::shared_ptr<SparseMatrix> leftop0 = leftBlock->get_op_array(CRE_DES).get_local_element(ij)[0]->getworkingrepresentation(leftBlock);//leftBlock->get_op_rep(CRE_DES_S0, ix, jx);
    SparseMatrix* dotop = 0;
        
    for (int kl =0; kl <rightBlock->get_op_array(CRE_DES).get_size(); kl++)
    {
      boost::shared_ptr<SparseMatrix> rightop2 = rightBlock->get_op_array(CRE_DES).get_local_element(kl)[1]->getworkingrepresentation(rightBlock);
      int kx = rightop2->get_orbs(0);
      int lx = rightop2->get_orbs(1);
      boost::shared_ptr<SparseMatrix> rightop0 = rightBlock->get_op_array(CRE_DES).get_local_element(kl)[0]->getworkingrepresentation(rightBlock);//rightBlock->get_op_rep(CRE_DES_S0, kx, lx);
      vector<double> expectations;
      spinExpectation(wave1, wave2, *leftop0, *dotop, *rightop0, big, expectations, false);
      spinExpectation(wave1, wave2, *leftop2, *dotop, *rightop2, big, expectations, false);
    
      vector<int> indices(4,0);
      indices[0] = ix; indices[1] = kx; indices[2] = jx; indices[3] = lx;
      spin_to_nonspin(indices, expectations, twopdm, CD_CD, true);
    }
  }
  }      
}

void compute_two_pdm_2_2_0(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm)
{
  SpinBlock* leftBlock = big.get_leftBlock()->get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();
  SpinBlock* dotBlock = big.get_leftBlock()->get_rightBlock();
  bool start = leftBlock->get_sites().size() ==2 ? true : false;
  int leftindex = leftBlock->get_sites()[0];

#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
  {
    int ccsize = leftBlock->get_op_array(CRE_CRE).get_size();
#ifdef _OPENMP
#pragma omp for schedule(guided) nowait
#endif
  for (int ij = 0; ij < ccsize; ++ij)
  {
    boost::shared_ptr<SparseMatrix> leftop2 = leftBlock->get_op_array(CRE_CRE).get_local_element(ij)[1]->getworkingrepresentation(leftBlock);
    int ix = leftop2->get_orbs(0);
    int jx = leftop2->get_orbs(1);
    boost::shared_ptr<SparseMatrix> leftop0 = leftBlock->get_op_array(CRE_CRE).get_local_element(ij)[0]->getworkingrepresentation(leftBlock);//leftBlock->get_op_rep(CRE_CRE_S0, ix, jx);
    SparseMatrix* rightop = 0;
        
    for (int kl =0; kl <dotBlock->get_op_array(CRE_CRE).get_size(); kl++)
    {
      boost::shared_ptr<SparseMatrix> dotop2 = dotBlock->get_op_array(CRE_CRE).get_local_element(kl)[1]->getworkingrepresentation(dotBlock);
      int kx = dotop2->get_orbs(0);
      int lx = dotop2->get_orbs(1);
      boost::shared_ptr<SparseMatrix> dotop0 = dotBlock->get_op_array(CRE_CRE).get_local_element(kl)[0]->getworkingrepresentation(dotBlock);//dotBlock->get_op_rep(CRE_CRE_S0, kx, lx);
      vector<double> expectations;
      Transposeview dop0 = Transposeview(*dotop0), dop2 = Transposeview(*dotop2); 
      spinExpectation(wave1, wave2, *leftop0, dop0, *rightop, big, expectations, false);
      spinExpectation(wave1, wave2, *leftop2, dop2, *rightop, big, expectations, false);
      vector<int> indices(4,0);
      indices[0] = ix; indices[1] = jx; indices[2] = lx; indices[3] = kx;
      expectations[0]*=-1; expectations[1]*=-1;
      spin_to_nonspin(indices, expectations, twopdm, CC_DD, true);
    }
  }      
  }


#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
  {
    int cdsize = leftBlock->get_op_array(CRE_DES).get_size();
#ifdef _OPENMP
#pragma omp for schedule(guided) nowait
#endif
  for (int ij = 0; ij < cdsize; ++ij)
  {
    boost::shared_ptr<SparseMatrix> leftop2 = leftBlock->get_op_array(CRE_DES).get_local_element(ij)[1]->getworkingrepresentation(leftBlock);
    int ix = leftop2->get_orbs(0);
    int jx = leftop2->get_orbs(1);
    boost::shared_ptr<SparseMatrix> leftop0 = leftBlock->get_op_array(CRE_DES).get_local_element(ij)[0]->getworkingrepresentation(leftBlock);//leftBlock->get_op_rep(CRE_DES_S0, ix, jx);
    SparseMatrix* rightop = 0;
        
    for (int kl =0; kl <dotBlock->get_op_array(CRE_DES).get_size(); kl++)
    {
      boost::shared_ptr<SparseMatrix> dotop2 = dotBlock->get_op_array(CRE_DES).get_local_element(kl)[1]->getworkingrepresentation(dotBlock);
      int kx = dotop2->get_orbs(0);
      int lx = dotop2->get_orbs(1);
      boost::shared_ptr<SparseMatrix> dotop0 = dotBlock->get_op_array(CRE_DES).get_local_element(kl)[0]->getworkingrepresentation(dotBlock);//dotBlock->get_op_rep(CRE_DES_S0, kx, lx);
      vector<double> expectations;
      spinExpectation(wave1, wave2, *leftop0, *dotop0, *rightop, big, expectations, false);
      spinExpectation(wave1, wave2, *leftop2, *dotop2, *rightop, big, expectations, false);
    
      vector<int> indices(4,0);
      indices[0] = ix; indices[1] = kx; indices[2] = jx; indices[3] = lx;
      spin_to_nonspin(indices, expectations, twopdm, CD_CD, true);
    }
  } 
  }     
}

void compute_two_pdm_1_1_2(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm)
{
  SpinBlock* leftBlock = big.get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();
  SpinBlock* dotBlock = leftBlock->get_rightBlock();

  int dotindex = dotBlock->get_sites()[0];

  for (int j = 0; j < leftBlock->get_leftBlock()->get_op_array(CRE).get_size(); ++j)
  {
    boost::shared_ptr<SparseMatrix> leftop = leftBlock->get_leftBlock()->get_op_array(CRE).get_local_element(j)[0];
    boost::shared_ptr<SparseMatrix> dotop  = dotBlock->get_op_rep(CRE, getSpinQuantum(dotindex), dotindex);
    int ix = leftop->get_orbs(0);
    int jx = dotindex;
        
    //#pragma omp parallel default(shared)
    {
      int ccsize = rightBlock->get_op_array(CRE_CRE).get_size();
      //#pragma omp for schedule(guided) nowait
    for (int kl =0; kl <ccsize; kl++)
    {
      boost::shared_ptr<SparseMatrix> rightop2 = rightBlock->get_op_array(CRE_CRE).get_local_element(kl)[1]->getworkingrepresentation(rightBlock);
      int kx = rightop2->get_orbs(0);
      int lx = rightop2->get_orbs(1);
      boost::shared_ptr<SparseMatrix> rightop0 = rightBlock->get_op_array(CRE_CRE).get_local_element(kl)[0]->getworkingrepresentation(rightBlock);//rightBlock->get_op_rep(CRE_CRE_S0, kx, lx);
      vector<double> expectations;
      Transposeview rop0 = Transposeview(*rightop0), rop2 = Transposeview(*rightop2); 
      spinExpectation(wave1, wave2, *leftop, *dotop, rop0, big, expectations, false);
      spinExpectation(wave1, wave2, *leftop, *dotop, rop2, big, expectations, false);
      vector<int> indices(4,0);
      indices[0] = ix; indices[1] = jx; indices[2] = kx; indices[3] = lx;
      expectations[0]*=-1; expectations[1]*=-1;
      spin_to_nonspin(indices, expectations, twopdm, CC_DD, true);
    }
    }
  }      


  for (int j = 0; j < leftBlock->get_leftBlock()->get_op_array(CRE).get_size(); ++j)
  {
    boost::shared_ptr<SparseMatrix> leftop = leftBlock->get_leftBlock()->get_op_array(CRE).get_local_element(j)[0];
    boost::shared_ptr<SparseMatrix> dotop  = dotBlock->get_op_rep(CRE, getSpinQuantum(dotindex) , dotindex);
    int ix = leftop->get_orbs(0);
    int jx = dotindex;
        
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
    {
      int cdsize = rightBlock->get_op_array(CRE_DES).get_size();
#ifdef _OPENMP
#pragma omp for schedule(guided) nowait
#endif
    for (int kl =0; kl <cdsize; kl++)
    {
      boost::shared_ptr<SparseMatrix> rightop2 = rightBlock->get_op_array(CRE_DES).get_local_element(kl)[1]->getworkingrepresentation(rightBlock);
      int kx = rightop2->get_orbs(0);
      int lx = rightop2->get_orbs(1);
      boost::shared_ptr<SparseMatrix> rightop0 = rightBlock->get_op_array(CRE_DES).get_local_element(kl)[0]->getworkingrepresentation(rightBlock);//rightBlock->get_op_rep(CRE_DES_S0, kx, lx);
      vector<double> expectations;
      Transposeview tdot = Transposeview(*dotop);
      spinExpectation(wave1, wave2, *leftop, tdot, *rightop0, big, expectations, true);
      spinExpectation(wave1, wave2, *leftop, tdot, *rightop2, big, expectations, true);
      vector<int> indices(4,0);
      vector<double> expect1(2), expect2(2);
      expect1[0] = expectations[0]; expect1[1] = expectations[2];
      indices[0] = ix; indices[1] = kx; indices[2] = jx; indices[3] = lx;
      spin_to_nonspin(indices, expect1, twopdm, CD_CD, true);

      expect2[0] = expectations[1]; expect2[1] = -expectations[3];
      indices[0] = jx; indices[1] = kx; indices[2] = ix; indices[3] = lx;
      spin_to_nonspin(indices, expect2, twopdm, CD_CD, true);
    }
    }
  }      
  
}



void compute_two_pdm_1_2_1(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm)
{
  SpinBlock* leftBlock = big.get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();
  SpinBlock* dotBlock = leftBlock->get_rightBlock();

  int dotindex = dotBlock->get_sites()[0];
  int jx = dotindex;

#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
  {
    int csize = leftBlock->get_leftBlock()->get_op_array(CRE).get_size();
#ifdef _OPENMP
#pragma omp for schedule(guided) nowait
#endif
  for (int j = 0; j < csize; ++j)
  {
    boost::shared_ptr<SparseMatrix> leftop = leftBlock->get_leftBlock()->get_op_array(CRE).get_local_element(j)[0];
    int ix = leftop->get_orbs(0);
    if (ix == jx ) continue;
    boost::shared_ptr<SparseMatrix> dotop0  = dotBlock->get_op_array(CRE_CRE).get_local_element(0)[0];
    boost::shared_ptr<SparseMatrix> dotop2  = dotBlock->get_op_array(CRE_CRE).get_local_element(0)[1];//dotBlock->get_op_rep(CRE_CRE_S2, dotindex, dotindex);
        
    for (int k =0; k <rightBlock->get_op_array(CRE).get_size(); k++)
    {
      SparseMatrix& rightop = *rightBlock->get_op_array(CRE).get_local_element(k)[0];
      int kx = rightop.get_orbs(0);
      Transposeview trop = Transposeview(rightop);
      Transposeview tlop = Transposeview(*leftop);

      vector<double> expectations;
      spinExpectation(wave1, wave2, tlop, *dotop0, trop, big, expectations, false);
      spinExpectation(wave1, wave2, tlop, *dotop2, trop, big, expectations, false);
      vector<int> indices(4,0);
      indices[0] = jx; indices[1] = jx; indices[2] = ix; indices[3] = kx;
      spin_to_nonspin(indices, expectations, twopdm, D_CC_D, true);
    }

    dotop0  = dotBlock->get_op_array(CRE_DES).get_local_element(0)[0];//dotBlock->get_op_rep(CRE_DES_S0, dotindex, dotindex);
    dotop2  = dotBlock->get_op_array(CRE_DES).get_local_element(0)[1];//dotBlock->get_op_rep(CRE_DES_S2, dotindex, dotindex);        
    for (int k =0; k <rightBlock->get_op_array(CRE).get_size(); k++)
    {
      SparseMatrix& rightop = *rightBlock->get_op_array(CRE).get_local_element(k)[0];
      int kx = rightop.get_orbs(0);
      Transposeview trop = Transposeview(rightop);
      Transposeview tlop = Transposeview(*leftop);

      vector<double> expectations;
      spinExpectation(wave1, wave2, *leftop, *dotop0, trop, big, expectations, false);
      spinExpectation(wave1, wave2, *leftop, *dotop2, trop, big, expectations, false);
      vector<int> indices(4,0);
      indices[0] = ix; indices[1] = jx; indices[2] = jx; indices[3] = kx;
      spin_to_nonspin(indices, expectations, twopdm, C_CD_D, true);

      expectations.resize(0);
      spinExpectation(wave1, wave2, tlop, *dotop0, rightop, big, expectations, false);
      spinExpectation(wave1, wave2, tlop, *dotop2, rightop, big, expectations, false);
      indices[0] = kx; indices[1] = jx; indices[2] = jx; indices[3] = ix;
      spin_to_nonspin(indices, expectations, twopdm, D_CD_C, true);

    }
  }      
  }
}

void compute_two_pdm_2_1_1(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm)
{
  SpinBlock* leftBlock = big.get_leftBlock()->get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();
  SpinBlock* dotBlock = big.get_leftBlock()->get_rightBlock();

  int dotindex = dotBlock->get_sites()[0];
  int jx = dotindex;
  int leftindex = leftBlock->get_sites()[0];
  int ix = leftindex;

  SparseMatrix& dotop = *dotBlock->get_op_array(CRE).get_local_element(0)[0];

  boost::shared_ptr<SparseMatrix> leftop0  = leftBlock->get_op_array(CRE_CRE).get_local_element(0)[0];
  boost::shared_ptr<SparseMatrix> leftop2  = leftBlock->get_op_array(CRE_CRE).get_local_element(0)[1];//leftBlock->get_op_rep(CRE_CRE_S2, leftindex, leftindex);
        
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
  {
    int csize = rightBlock->get_op_array(CRE).get_size();
#ifdef _OPENMP
#pragma omp for schedule(guided) nowait
#endif
  for (int k =0; k <csize; k++)
    {
      SparseMatrix& rightop = *rightBlock->get_op_array(CRE).get_local_element(k)[0];
      int kx = rightop.get_orbs(0);
      Transposeview trop = Transposeview(rightop);
      Transposeview tdop = Transposeview(dotop);
      
      vector<double> expectations;
      spinExpectation(wave1, wave2, *leftop0, tdop, trop, big, expectations, false);
      spinExpectation(wave1, wave2, *leftop2, tdop, trop, big, expectations, false);
      vector<int> indices(4,0);
      indices[0] = ix; indices[1] = ix; indices[2] = jx; indices[3] = kx;
      spin_to_nonspin(indices, expectations, twopdm, CC_D_D, true);
    }
  }

  leftop0  = leftBlock->get_op_array(CRE_DES).get_local_element(0)[0];//leftBlock->get_op_rep(CRE_DES_S0, leftindex, leftindex);
  leftop2  = leftBlock->get_op_array(CRE_DES).get_local_element(0)[1];//leftBlock->get_op_rep(CRE_DES_S2, leftindex, leftindex);        
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
  {
    int csize = rightBlock->get_op_array(CRE).get_size();
#ifdef _OPENMP
#pragma omp for schedule(guided) nowait
#endif
  for (int k =0; k <csize; k++)
    {
      SparseMatrix& rightop = *rightBlock->get_op_array(CRE).get_local_element(k)[0];
      int kx = rightop.get_orbs(0);
      Transposeview trop = Transposeview(rightop);
      Transposeview tdop = Transposeview(dotop);

      vector<double> expectations;
      spinExpectation(wave1, wave2, *leftop0, dotop, trop, big, expectations, false);
      spinExpectation(wave1, wave2, *leftop2, dotop, trop, big, expectations, false);
      vector<int> indices(4,0);
      indices[0] = ix; indices[1] = jx; indices[2] = ix; indices[3] = kx;
      spin_to_nonspin(indices, expectations, twopdm, CD_CD, true);

    }
  }
}

void compute_two_pdm_0_3_1(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm)
{
  SpinBlock* leftBlock = big.get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();
  SpinBlock* dotBlock = leftBlock->get_rightBlock();

  int dotindex = dotBlock->get_sites()[0];
  int jx = dotindex;

  boost::shared_ptr<SparseMatrix> dotop0 = dotBlock->get_op_array(CRE).get_local_element(0)[0];//dotBlock->get_op_rep(CRE_S1, dotindex);
  boost::shared_ptr<SparseMatrix> dotop1  = dotBlock->get_op_array(CRE_CRE).get_local_element(0)[0];//dotBlock->get_op_rep(CRE_CRE_S0, dotindex, dotindex);
  boost::shared_ptr<SparseMatrix> dotop2  = dotBlock->get_op_array(CRE_CRE).get_local_element(0)[1];//dotBlock->get_op_rep(CRE_CRE_S2, dotindex, dotindex);

  Cre Dotop1, Dotop2;
  Dotop1.set_orbs() = dotop0->get_orbs(); Dotop1.set_orbs().push_back(dotindex); Dotop1.set_orbs().push_back(dotindex);
  Dotop1.set_initialised() = true;
  Dotop1.set_fermion() = true;
  Dotop1.set_deltaQuantum(1, (dotop1->get_deltaQuantum(0) - dotop0->get_deltaQuantum(0))[0]);
  Dotop1.allocate(dotBlock->get_stateInfo());
  operatorfunctions::Product(dotBlock, *dotop1, Transposeview(*dotop0), Dotop1, 1.0);

  Dotop2.set_orbs() = dotop0->get_orbs(); Dotop2.set_orbs().push_back(dotindex); Dotop2.set_orbs().push_back(dotindex);
  Dotop2.set_initialised() = true;
  Dotop2.set_fermion() = true;
  Dotop2.set_deltaQuantum(1, (dotop2->get_deltaQuantum(0) - dotop0->get_deltaQuantum(0))[0]);
  Dotop2.allocate(dotBlock->get_stateInfo());
  operatorfunctions::Product(dotBlock, *dotop2, Transposeview(*dotop0), Dotop2, 1.0);
        
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
  {
  SparseMatrix *leftop = 0;
  int csize = rightBlock->get_op_array(CRE).get_size();
#ifdef _OPENMP
#pragma omp for schedule(guided) nowait
#endif
  for (int k =0; k <csize; k++)
    {
      SparseMatrix& rightop = *rightBlock->get_op_array(CRE).get_local_element(k)[0];
      int kx = rightop.get_orbs(0);
      Transposeview trop = Transposeview(rightop);

      vector<double> expectations;
      spinExpectation(wave1, wave2, *leftop, Dotop1, trop, big, expectations, false);
      spinExpectation(wave1, wave2, *leftop, Dotop2, trop, big, expectations, false);
      vector<int> indices(4,0);
      indices[0] = jx; indices[1] = jx; indices[2] = jx; indices[3] = kx;
      spin_to_nonspin(indices, expectations, twopdm, CC_D_D, true);
    }
  }
  /*
  dotop1  = dotBlock->get_op_rep(CRE_DES_S0, dotindex, dotindex);
  dotop2  = dotBlock->get_op_rep(CRE_DES_S2, dotindex, dotindex);        
  Dotop1.set_orbs() = dotop0->get_orbs(); Dotop1.set_orbs().push_back(dotindex); Dotop1.set_orbs().push_back(dotindex);
  Dotop1.set_initialised() = true;
  Dotop1.set_fermion() = true;
  Dotop1.set_deltaQuantum() = ( dotop1->get_deltaQuantum() - dotop0->get_deltaQuantum() )[0];
  Dotop1.allocate(dotBlock->get_stateInfo());
  operatorfunctions::Product(dotBlock, Transposeview(*dotop0), *dotop1, Dotop1, 1.0);

  Dotop2.set_orbs() = dotop0->get_orbs(); Dotop2.set_orbs().push_back(dotindex); Dotop2.set_orbs().push_back(dotindex);
  Dotop2.set_initialised() = true;
  Dotop2.set_fermion() = true;
  Dotop2.set_deltaQuantum() = ( dotop2->get_deltaQuantum() - dotop0->get_deltaQuantum() )[0];
  Dotop2.allocate(dotBlock->get_stateInfo());
  operatorfunctions::Product(dotBlock, Transposeview(*dotop0), *dotop2, Dotop2, 1.0);

  for (int k =0; k <rightBlock->get_op_array(CRE_S1).get_size(); k++)
    {
      SparseMatrix& rightop = rightBlock->get_op_array(CRE_S1).get_local_element(k);
      int kx = rightop.get_orbs(0);
      Transposeview trop = Transposeview(rightop);
      Transposeview tlop = Transposeview(*leftop);
      
      vector<double> expectations;
      spinExpectation(wave1, wave2, *leftop, Dotop1, rightop, big, expectations, false);
      spinExpectation(wave1, wave2, *leftop, Dotop2, rightop, big, expectations, false);
      vector<int> indices(4,0);
      indices[0] = kx; indices[1] = jx; indices[2] = jx; indices[3] = jx;
      spin_to_nonspin(indices, expectations, twopdm, D_CD_C, false);

    }
  */
}


void compute_two_pdm_0_3_1_notranspose(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm)
{
  SpinBlock* leftBlock = big.get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();
  SpinBlock* dotBlock = leftBlock->get_rightBlock();

  int dotindex = dotBlock->get_sites()[0];
  int jx = dotindex;

  boost::shared_ptr<SparseMatrix> dotop0 = dotBlock->get_op_array(CRE).get_local_element(0)[0];//dotBlock->get_op_rep(CRE_S1, dotindex);
  boost::shared_ptr<SparseMatrix> dotop0d = dotBlock->get_op_array(DES).get_local_element(0)[0];//dotBlock->get_op_rep(CRE_S1, dotindex);
  boost::shared_ptr<SparseMatrix> dotop1  = dotBlock->get_op_array(CRE_CRE).get_local_element(0)[0];//dotBlock->get_op_rep(CRE_CRE_S0, dotindex, dotindex);
  boost::shared_ptr<SparseMatrix> dotop2  = dotBlock->get_op_array(CRE_CRE).get_local_element(0)[1];//dotBlock->get_op_rep(CRE_CRE_S2, dotindex, dotindex);

  Cre Dotop1, Dotop2;
  Dotop1.set_orbs() = dotop0->get_orbs(); Dotop1.set_orbs().push_back(dotindex); Dotop1.set_orbs().push_back(dotindex);
  Dotop1.set_initialised() = true;
  Dotop1.set_fermion() = true;
  Dotop1.set_deltaQuantum(1, (dotop1->get_deltaQuantum(0) - dotop0->get_deltaQuantum(0))[0]);
  Dotop1.allocate(dotBlock->get_stateInfo());
  operatorfunctions::Product(dotBlock, *dotop1, *dotop0d, Dotop1, 1.0);

  Dotop2.set_orbs() = dotop0->get_orbs(); Dotop2.set_orbs().push_back(dotindex); Dotop2.set_orbs().push_back(dotindex);
  Dotop2.set_initialised() = true;
  Dotop2.set_fermion() = true;
  Dotop2.set_deltaQuantum(1, (dotop2->get_deltaQuantum(0) - dotop0->get_deltaQuantum(0))[0]);
  Dotop2.allocate(dotBlock->get_stateInfo());
  operatorfunctions::Product(dotBlock, *dotop2, *dotop0d, Dotop2, 1.0);
        
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
  {
  SparseMatrix *leftop = 0;
  int dsize = rightBlock->get_op_array(DES).get_size();
#ifdef _OPENMP
#pragma omp for schedule(guided) nowait
#endif
  for (int k =0; k <dsize; k++)
    {
      SparseMatrix& rightop = *rightBlock->get_op_array(DES).get_local_element(k)[0];
      int kx = rightop.get_orbs(0);

      vector<double> expectations;
      spinExpectation(wave1, wave2, *leftop, Dotop1, rightop, big, expectations, false);
      spinExpectation(wave1, wave2, *leftop, Dotop2, rightop, big, expectations, false);
      vector<int> indices(4,0);
      indices[0] = jx; indices[1] = jx; indices[2] = jx; indices[3] = kx;
      spin_to_nonspin(indices, expectations, twopdm, CC_D_D, false);
    }
  }


  //******
  {
    dotop0 = dotBlock->get_op_array(CRE).get_local_element(0)[0];//dotBlock->get_op_rep(CRE_S1, dotindex);
    dotop0d = dotBlock->get_op_array(DES).get_local_element(0)[0];//dotBlock->get_op_rep(CRE_S1, dotindex);
    dotop1  = dotBlock->get_op_array(DES_DES).get_local_element(0)[0];//dotBlock->get_op_rep(CRE_CRE_S0, dotindex, dotindex);
    dotop2  = dotBlock->get_op_array(DES_DES).get_local_element(0)[1];//dotBlock->get_op_rep(CRE_CRE_S2, dotindex, dotindex);
    
    Cre Dotop1, Dotop2;
    Dotop1.set_orbs() = dotop0->get_orbs(); Dotop1.set_orbs().push_back(dotindex); Dotop1.set_orbs().push_back(dotindex);
    Dotop1.set_initialised() = true;
    Dotop1.set_fermion() = true;
    Dotop1.set_deltaQuantum(1, (dotop1->get_deltaQuantum(0) + dotop0->get_deltaQuantum(0))[0]);
    Dotop1.allocate(dotBlock->get_stateInfo());
    operatorfunctions::Product(dotBlock, *dotop0, *dotop1, Dotop1, 1.0);
    
    Dotop2.set_orbs() = dotop0->get_orbs(); Dotop2.set_orbs().push_back(dotindex); Dotop2.set_orbs().push_back(dotindex);
    Dotop2.set_initialised() = true;
    Dotop2.set_fermion() = true;
    Dotop2.set_deltaQuantum(1, (dotop2->get_deltaQuantum(0) + dotop0->get_deltaQuantum(0))[0]);
    Dotop2.allocate(dotBlock->get_stateInfo());
    operatorfunctions::Product(dotBlock, *dotop0, *dotop2, Dotop2, 1.0);
    
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
    {
      SparseMatrix *leftop = 0;
      int dsize = rightBlock->get_op_array(CRE).get_size();
#ifdef _OPENMP
#pragma omp for schedule(guided) nowait
#endif
      for (int k =0; k <dsize; k++)
	{
	  SparseMatrix& rightop = *rightBlock->get_op_array(CRE).get_local_element(k)[0];
	  int kx = rightop.get_orbs(0);
	  
	  vector<double> expectations;
	  spinExpectation(wave1, wave2, *leftop, Dotop1, rightop, big, expectations, false);
	  spinExpectation(wave1, wave2, *leftop, Dotop2, rightop, big, expectations, false);
	  vector<int> indices(4,0);
	  indices[0] = kx; indices[1] = jx; indices[2] = jx; indices[3] = jx;
	  spin_to_nonspin(indices, expectations, twopdm, CC_D_D, false);
	}
    }
  }

}


void compute_two_pdm_3_0_1(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm)
{
  SpinBlock* leftBlock = big.get_leftBlock()->get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();
  SpinBlock* dotBlock = big.get_leftBlock()->get_rightBlock();

  int leftindex = leftBlock->get_sites()[0];
  int jx = leftindex;

  boost::shared_ptr<SparseMatrix> leftop0 = leftBlock->get_op_array(CRE).get_local_element(0)[0];//leftBlock->get_op_rep(CRE_S1, leftindex);
  boost::shared_ptr<SparseMatrix> leftop1  = leftBlock->get_op_array(CRE_CRE).get_local_element(0)[0];//leftBlock->get_op_rep(CRE_CRE_S0, leftindex, leftindex);
  boost::shared_ptr<SparseMatrix> leftop2  = leftBlock->get_op_array(CRE_CRE).get_local_element(0)[1];//leftBlock->get_op_rep(CRE_CRE_S2, leftindex, leftindex);

  Cre Leftop1, Leftop2;
  Leftop1.set_orbs() = leftop0->get_orbs(); Leftop1.set_orbs().push_back(leftindex); Leftop1.set_orbs().push_back(leftindex);
  Leftop1.set_initialised() = true;
  Leftop1.set_fermion() = true;
  Leftop1.set_deltaQuantum(1, (leftop1->get_deltaQuantum(0) - leftop0->get_deltaQuantum(0))[0]);
  Leftop1.allocate(leftBlock->get_stateInfo());
  operatorfunctions::Product(leftBlock, *leftop1, Transposeview(*leftop0), Leftop1, 1.0);

  Leftop2.set_orbs() = leftop0->get_orbs(); Leftop2.set_orbs().push_back(leftindex); Leftop2.set_orbs().push_back(leftindex);
  Leftop2.set_initialised() = true;
  Leftop2.set_fermion() = true;
  Leftop2.set_deltaQuantum(1, (leftop2->get_deltaQuantum(0) - leftop0->get_deltaQuantum(0))[0]);
  Leftop2.allocate(leftBlock->get_stateInfo());
  operatorfunctions::Product(leftBlock, *leftop2, Transposeview(*leftop0), Leftop2, 1.0);
        
#ifdef _OPENMP
#pragma omp parallel default(shared) 
#endif
  {
  SparseMatrix *dotop;
  int cresize = rightBlock->get_op_array(CRE).get_size();
#ifdef _OPENMP
#pragma omp for schedule(guided) nowait private(dotop)
#endif
  for (int k =0; k <cresize; k++)
    {
      SparseMatrix& rightop = *rightBlock->get_op_array(CRE).get_local_element(k)[0];
      int kx = rightop.get_orbs(0);
      Transposeview trop = Transposeview(rightop);

      vector<double> expectations;
      dotop = 0;
      spinExpectation(wave1, wave2, Leftop1, *dotop, trop, big, expectations, false);
      spinExpectation(wave1, wave2, Leftop2, *dotop, trop, big, expectations, false);
      vector<int> indices(4,0);
      indices[0] = jx; indices[1] = jx; indices[2] = jx; indices[3] = kx;
      spin_to_nonspin(indices, expectations, twopdm, CC_D_D, true);
    }
  }
}

void compute_two_pdm_3_1_0(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm)
{
  bool w1_eq_w2 = &wave1==&wave2 ? true : false;
  SpinBlock* leftBlock = big.get_leftBlock()->get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();
  SpinBlock* dotBlock = big.get_leftBlock()->get_rightBlock();

  int leftindex = leftBlock->get_sites()[0];
  int jx = leftindex;

  boost::shared_ptr<SparseMatrix> leftop0 = leftBlock->get_op_array(CRE).get_local_element(0)[0];//leftBlock->get_op_rep(CRE_S1, leftindex);
  boost::shared_ptr<SparseMatrix> leftop1  = leftBlock->get_op_array(CRE_CRE).get_local_element(0)[0];//leftBlock->get_op_rep(CRE_CRE_S0, leftindex, leftindex);
  boost::shared_ptr<SparseMatrix> leftop2  = leftBlock->get_op_array(CRE_CRE).get_local_element(0)[1];//leftBlock->get_op_rep(CRE_CRE_S2, leftindex, leftindex);

  Cre Leftop1, Leftop2;
  Leftop1.set_orbs() = leftop0->get_orbs(); Leftop1.set_orbs().push_back(leftindex); Leftop1.set_orbs().push_back(leftindex);
  Leftop1.set_initialised() = true;
  Leftop1.set_fermion() = true;
  Leftop1.set_deltaQuantum(1, (leftop1->get_deltaQuantum(0) - leftop0->get_deltaQuantum(0))[0]);
  Leftop1.allocate(leftBlock->get_stateInfo());
  operatorfunctions::Product(leftBlock, *leftop1, Transposeview(*leftop0), Leftop1, 1.0);

  Leftop2.set_orbs() = leftop0->get_orbs(); Leftop2.set_orbs().push_back(leftindex); Leftop2.set_orbs().push_back(leftindex);
  Leftop2.set_initialised() = true;
  Leftop2.set_fermion() = true;
  Leftop2.set_deltaQuantum(1, (leftop2->get_deltaQuantum(0) - leftop0->get_deltaQuantum(0))[0]);
  Leftop2.allocate(leftBlock->get_stateInfo());
  operatorfunctions::Product(leftBlock, *leftop2, Transposeview(*leftop0), Leftop2, 1.0);

  SparseMatrix *rightop = 0;
  {
    SparseMatrix& dotop = *dotBlock->get_op_array(CRE).get_local_element(0)[0];
    int kx = dotop.get_orbs(0);
    Transposeview tdop = Transposeview(dotop);
    
    vector<double> expectations;

    spinExpectation(wave1, wave2, Leftop1, tdop, *rightop, big, expectations, false);
    spinExpectation(wave1, wave2, Leftop2, tdop, *rightop, big, expectations, false);
    vector<int> indices(4,0);
    indices[0] = jx; indices[1] = jx; indices[2] = jx; indices[3] = kx;

    spin_to_nonspin(indices, expectations, twopdm, CC_D_D, true);
  }

  /*
  leftop1  = leftBlock->get_op_rep(CRE_DES_S0, leftindex, leftindex);
  leftop2  = leftBlock->get_op_rep(CRE_DES_S2, leftindex, leftindex);        
  Leftop1.set_orbs() = leftop0->get_orbs(); Leftop1.set_orbs().push_back(leftindex); Leftop1.set_orbs().push_back(leftindex);
  Leftop1.set_initialised() = true;
  Leftop1.set_fermion() = true;
  Leftop1.set_deltaQuantum() = ( leftop1->get_deltaQuantum() - leftop0->get_deltaQuantum() )[0];
  Leftop1.allocate(leftBlock->get_stateInfo());
  operatorfunctions::Product(leftBlock, *leftop1, Transposeview(*leftop0), Leftop1, 1.0);

  Leftop2.set_orbs() = leftop0->get_orbs(); Leftop2.set_orbs().push_back(leftindex); Leftop2.set_orbs().push_back(leftindex);
  Leftop2.set_initialised() = true;
  Leftop2.set_fermion() = true;
  Leftop2.set_deltaQuantum() = ( leftop2->get_deltaQuantum() - leftop0->get_deltaQuantum() )[0];
  Leftop2.allocate(leftBlock->get_stateInfo());
  operatorfunctions::Product(leftBlock, *leftop2, Transposeview(*leftop0), Leftop2, 1.0);


  {
    SparseMatrix& dotop = dotBlock->get_op_array(CRE_S1).get_local_element(0);
    int kx = dotop.get_orbs(0);
    
    vector<double> expectations;
    spinExpectation(wave1, wave2, Leftop1, dotop, *rightop, big, expectations, false);
    spinExpectation(wave1, wave2, Leftop2, dotop, *rightop, big, expectations, false);
    vector<int> indices(4,0);
    indices[0] = jx; indices[1] = kx; indices[2] = jx; indices[3] = jx;
    spin_to_nonspin(indices, expectations, twopdm, CD_D_C, false);
    
  }
  */
}


void compute_two_pdm_1_3_0(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm)
{
  SpinBlock* leftBlock = big.get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();
  SpinBlock* dotBlock = leftBlock->get_rightBlock();

  int dotindex = dotBlock->get_sites()[0];
  int jx = dotindex;

  boost::shared_ptr<SparseMatrix> dotop0 = dotBlock->get_op_array(CRE).get_local_element(0)[0];//dotBlock->get_op_rep(CRE_S1, dotindex);
  boost::shared_ptr<SparseMatrix> dotop1  = dotBlock->get_op_array(CRE_CRE).get_local_element(0)[0];//dotBlock->get_op_rep(CRE_CRE_S0, dotindex, dotindex);
  boost::shared_ptr<SparseMatrix> dotop2  = dotBlock->get_op_array(CRE_CRE).get_local_element(0)[1];//dotBlock->get_op_rep(CRE_CRE_S2, dotindex, dotindex);

  Cre Dotop1, Dotop2;
  Dotop1.set_orbs() = dotop0->get_orbs(); Dotop1.set_orbs().push_back(dotindex); Dotop1.set_orbs().push_back(dotindex);
  Dotop1.set_initialised() = true;
  Dotop1.set_fermion() = true;
  Dotop1.set_deltaQuantum(1, (dotop1->get_deltaQuantum(0) - dotop0->get_deltaQuantum(0))[0]);
  Dotop1.allocate(dotBlock->get_stateInfo());
  operatorfunctions::Product(dotBlock, *dotop1, Transposeview(*dotop0), Dotop1, 1.0);

  Dotop2.set_orbs() = dotop0->get_orbs(); Dotop2.set_orbs().push_back(dotindex); Dotop2.set_orbs().push_back(dotindex);
  Dotop2.set_initialised() = true;
  Dotop2.set_fermion() = true;
  Dotop2.set_deltaQuantum(1, (dotop2->get_deltaQuantum(0) - dotop0->get_deltaQuantum(0))[0]);
  Dotop2.allocate(dotBlock->get_stateInfo());
  operatorfunctions::Product(dotBlock, *dotop2, Transposeview(*dotop0), Dotop2, 1.0);
        
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
  {
  SparseMatrix *rightop = 0;
  int csize = leftBlock->get_leftBlock()->get_op_array(CRE).get_size();
#ifdef _OPENMP
#pragma omp for schedule(guided) nowait
#endif
  for (int k =0; k <csize; k++)
    {
      SparseMatrix& leftop = *leftBlock->get_leftBlock()->get_op_array(CRE).get_local_element(k)[0];
      int kx = leftop.get_orbs(0);
      Transposeview tlop = Transposeview(leftop);

      vector<double> expectations;
      spinExpectation(wave1, wave2, tlop, Dotop1, *rightop, big, expectations, false);
      spinExpectation(wave1, wave2, tlop, Dotop2, *rightop, big, expectations, false);

      vector<int> indices(4,0);
      indices[0] = jx; indices[1] = jx; indices[2] = jx; indices[3] = kx;
      spin_to_nonspin(indices, expectations, twopdm, CC_D_D, true);
    }
  }
  /*
  dotop1  = dotBlock->get_op_rep(CRE_DES_S0, dotindex, dotindex);
  dotop2  = dotBlock->get_op_rep(CRE_DES_S2, dotindex, dotindex);
  Dotop1.set_deltaQuantum() = ( dotop1->get_deltaQuantum() - dotop0->get_deltaQuantum() )[0];
  Dotop1.allocate(dotBlock->get_stateInfo());
  operatorfunctions::Product(dotBlock, *dotop1, Transposeview(*dotop0), Dotop1, 1.0);

  Dotop2.set_deltaQuantum() = ( dotop2->get_deltaQuantum() - dotop0->get_deltaQuantum() )[0];
  Dotop2.allocate(dotBlock->get_stateInfo());
  operatorfunctions::Product(dotBlock, *dotop2, Transposeview(*dotop0), Dotop2, 1.0);
        
  for (int k =0; k <leftBlock->get_leftBlock()->get_op_array(CRE_S1).get_size(); k++)
    {
      SparseMatrix& leftop = leftBlock->get_leftBlock()->get_op_array(CRE_S1).get_local_element(k);
      int kx = leftop.get_orbs(0);

      vector<double> expectations;
      spinExpectation(wave1, wave2, leftop, Dotop1, *rightop, big, expectations, false);
      spinExpectation(wave1, wave2, leftop, Dotop2, *rightop, big, expectations, false);

      vector<int> indices(4,0);
      indices[0] = kx; indices[1] = jx; indices[2] = jx; indices[3] = jx;
      spin_to_nonspin(indices, expectations, twopdm, C_CD_D, false);
    }
  */
}

void compute_two_pdm_1_3(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, array_4d<double>& twopdm)
{
  SpinBlock* leftBlock = big.get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();
  SpinBlock* dotBlock = leftBlock->get_rightBlock();

  int rightindex = rightBlock->get_sites()[0];
  int jx = rightindex;

  boost::shared_ptr<SparseMatrix> rightop0 = rightBlock->get_op_array(CRE).get_local_element(0)[0];//rightBlock->get_op_rep(CRE_S1, rightindex);
  boost::shared_ptr<SparseMatrix> rightop1  = rightBlock->get_op_array(CRE_DES).get_local_element(0)[0];//rightBlock->get_op_rep(CRE_DES_S0, rightindex, rightindex);
  boost::shared_ptr<SparseMatrix> rightop2  = rightBlock->get_op_array(CRE_DES).get_local_element(0)[1];//rightBlock->get_op_rep(CRE_DES_S2, rightindex, rightindex);

  Cre Rightop1, Rightop2;
  Rightop1.set_orbs() = rightop0->get_orbs(); Rightop1.set_orbs().push_back(rightindex); Rightop1.set_orbs().push_back(rightindex);
  Rightop1.set_initialised() = true;
  Rightop1.set_fermion() = true;
  Rightop1.set_deltaQuantum(1, (rightop1->get_deltaQuantum(0) - rightop0->get_deltaQuantum(0))[0]);
  Rightop1.allocate(rightBlock->get_stateInfo());
  operatorfunctions::Product(rightBlock, *rightop1, Transposeview(*rightop0), Rightop1, 1.0);

  Rightop2.set_orbs() = rightop0->get_orbs(); Rightop2.set_orbs().push_back(rightindex); Rightop2.set_orbs().push_back(rightindex);
  Rightop2.set_initialised() = true;
  Rightop2.set_fermion() = true;
  Rightop2.set_deltaQuantum(1, (rightop2->get_deltaQuantum(0) - rightop0->get_deltaQuantum(0))[0]);
  Rightop2.allocate(rightBlock->get_stateInfo());
  operatorfunctions::Product(rightBlock, *rightop2, Transposeview(*rightop0), Rightop2, 1.0);
        
  for (int k =0; k <leftBlock->get_leftBlock()->get_op_array(CRE).get_size(); k++)
    {
      SparseMatrix& leftop = *leftBlock->get_leftBlock()->get_op_array(CRE).get_local_element(k)[0];
      int kx = leftop.get_orbs(0);
      Transposeview tlop = Transposeview(leftop);
      SparseMatrix* dotop = 0;

      vector<double> expectations;
      spinExpectation(wave1, wave2, leftop, *dotop, Rightop1, big, expectations, false);
      spinExpectation(wave1, wave2, leftop, *dotop, Rightop2, big, expectations, false);

      vector<int> indices(4,0);
      indices[0] = kx; indices[1] = jx; indices[2] = jx; indices[3] = jx;
      spin_to_nonspin(indices, expectations, twopdm, C_CD_D, true);
    }

  SparseMatrix& dotop = *dotBlock->get_op_array(CRE).get_local_element(0)[0];
  int kx = dotop.get_orbs(0);
  Transposeview tdop = Transposeview(dotop);
  SparseMatrix* leftop = 0;
  
  vector<double> expectations;
  spinExpectation(wave1, wave2, *leftop, dotop, Rightop1, big, expectations, false);
  spinExpectation(wave1, wave2, *leftop, dotop, Rightop2, big, expectations, false);
  
  vector<int> indices(4,0);
  indices[0] = kx; indices[1] = jx; indices[2] = jx; indices[3] = jx;
  spin_to_nonspin(indices, expectations, twopdm, C_CD_D, true);

}



std::vector<int> distribute_procs(const int numprocs, const int numjobs) {
  // In this we are going to distribute the jobs between the procs.
  std::vector<int> retval(numjobs);

  const int residual = numjobs%numprocs;
  const int jobsperproc = numjobs/numprocs;

  // bad method... will think abt later
  for(int i = 0; i < jobsperproc; ++i)
    for(int j = 0; j < numprocs; ++j)
      retval[i*numprocs+j] = j;

  for(int i = 0; i < residual; ++i)
    retval[(jobsperproc-1)*numprocs+numprocs + i]= i;

  return retval;
}

}
