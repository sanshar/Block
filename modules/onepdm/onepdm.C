/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "onepdm.h"
#include "MatrixBLAS.h"
#include <boost/format.hpp>
#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include "operatorfunctions.h"
#include "execinfo.h"
#include "include/newmatutils.h"
#include "pario.h"

namespace SpinAdapted{
void compute_onepdm(std::vector<Wavefunction>& wavefunctions, const SpinBlock& system, const SpinBlock& systemDot, const SpinBlock& newSystem, const SpinBlock& newEnvironment, const SpinBlock& big, const int numprocs)
{

  const int nroots = wavefunctions.size();
  for(int i=0;i<nroots;++i)
    for(int j=0;j<=i;++j)
      {
	Matrix onepdm;
	load_onepdm_binary(onepdm, i ,j);
	Wavefunction &wavefunction1 = wavefunctions[i];
	Wavefunction &wavefunction2 = wavefunctions[j];

	pout <<"\t\t\t Performing sweep calculation: 1PDM "<<endl;

	//if (big.get_leftBlock()->size() == 2) {
	p2out << "\t\t\t compute 2_0 "<<mpigetrank()<<endl;
	compute_one_pdm_2_0_0(wavefunction1, wavefunction2, big, onepdm);
	compute_one_pdm_0_2_0(wavefunction1, wavefunction2, big, onepdm);
	compute_one_pdm_1_1_0(wavefunction1, wavefunction2, big, onepdm);
	//}
	
	//if (big.get_rightBlock()->size() == 1) {
	p2out << "\t\t\t compute 0_2 "<<mpigetrank()<<endl;
	compute_one_pdm_0_2(wavefunction1, wavefunction2, big, onepdm);
	//}
	
	p2out << "\t\t\t compute 1_1 "<<mpigetrank()<<endl;
	compute_one_pdm_1_1(wavefunction1, wavefunction2, big, onepdm);
	accumulate_onepdm(onepdm);
	save_onepdm_binary(onepdm, i, j);
      }
}

void compute_one_pdm_2_0_0(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, Matrix& onepdm)
{
  SpinBlock* leftBlock = big.get_leftBlock()->get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();
  SpinBlock* dotBlock = big.get_leftBlock()->get_rightBlock();

  //this is function 2_0_0
  for (int ij = 0; ij < leftBlock->get_op_array(CRE_DES).get_size(); ++ij)
  {
    boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_array(CRE_DES).get_local_element(ij)[0]->getworkingrepresentation(leftBlock);//spin 0
    int ix = op->get_orbs(0);
    int jx = op->get_orbs(1);

    //now for a leftop that sits on the system block, this is done by tensortrace function
    Cre leftop1;
    leftop1.set_orbs() = op->get_orbs(); 
    leftop1.set_initialised() = true;
    leftop1.set_fermion() = false;
    leftop1.set_deltaQuantum(1, op->get_deltaQuantum(0));
    leftop1.allocate(big.get_leftBlock()->get_stateInfo());
    operatorfunctions::TensorTrace(leftBlock, *op, big.get_leftBlock(), &(big.get_leftBlock()->get_stateInfo()), leftop1);

    Wavefunction opw2;
    vector<SpinQuantum> dQ = wave1.get_deltaQuantum();
    opw2.initialise(dQ, &big, true);
    operatorfunctions::TensorMultiply(big.get_leftBlock(), leftop1, &big, wave2, opw2, dQ[0], 1.0); // FIXME dQ is actually never used
    double sum = sqrt(2.0)*DotProduct(wave1, opw2);

    if(dmrginp.spinAdapted()) {
      onepdm(2*ix+1, 2*jx+1) = (sum)/2.0;
      onepdm(2*ix+2, 2*jx+2) = (sum)/2.0;
      onepdm(2*jx+1, 2*ix+1) = (sum)/2.0;
      onepdm(2*jx+2, 2*ix+2) = (sum)/2.0;
    }
    else {
      onepdm(ix+1, jx+1) = sum/sqrt(2.0);
      onepdm(jx+1, ix+1) = sum/sqrt(2.0);
    }
  }      
}

void compute_pair_2_0_0(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, Matrix& onepdm) {
  SpinBlock* leftBlock = big.get_leftBlock()->get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();
  SpinBlock* dotBlock = big.get_leftBlock()->get_rightBlock();

  for (int ij = 0; ij < leftBlock->get_op_array(CRE_CRE).get_size(); ++ij) {
    boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_array(CRE_CRE).get_local_element(ij)[0]->getworkingrepresentation(leftBlock);//spin 0
    int ix = op->get_orbs(0);
    int jx = op->get_orbs(1);

    Transposeview top(op);
    Cre leftop1;
    leftop1.set_orbs() = top.get_orbs(); 
    leftop1.set_initialised() = true;
    leftop1.set_fermion() = false;
    leftop1.set_deltaQuantum(1, top.get_deltaQuantum(0));
    leftop1.allocate(big.get_leftBlock()->get_stateInfo());
    operatorfunctions::TensorTrace(leftBlock, top, big.get_leftBlock(), &(big.get_leftBlock()->get_stateInfo()), leftop1);

    Wavefunction opw2;
    vector<SpinQuantum> dQ = wave1.get_deltaQuantum();
    opw2.initialise(dQ, &big, true);
    operatorfunctions::TensorMultiply(big.get_leftBlock(), leftop1, &big, wave2, opw2, dQ[0], 1.0);
    double sum = sqrt(2.0)*DotProduct(wave1, opw2);

    if(dmrginp.spinAdapted()) {
      pout << "BCS with spin adaption not implemented yet." << endl;
    }
    else {
      onepdm(ix+1, jx+1) = -sum/sqrt(2.0);
      onepdm(jx+1, ix+1) = sum/sqrt(2.0);
    }
  }
}

void compute_one_pdm_0_2_0(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, Matrix& onepdm)
{
  SpinBlock* leftBlock = big.get_leftBlock()->get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();
  SpinBlock* dotBlock = big.get_leftBlock()->get_rightBlock();


  for (int ij = 0; ij < dotBlock->get_op_array(CRE_DES).get_size(); ++ij)
  {
    boost::shared_ptr<SparseMatrix> op = dotBlock->get_op_array(CRE_DES).get_local_element(ij)[0]->getworkingrepresentation(leftBlock);//spin 0
    int ix = op->get_orbs(0);
    int jx = op->get_orbs(1);

    //now for a leftop that sits on the system block, this is done by tensortrace function
    Cre leftop1;
    leftop1.set_orbs() = op->get_orbs(); 
    leftop1.set_initialised() = true;
    leftop1.set_fermion() = false;
    leftop1.set_deltaQuantum(1, op->get_deltaQuantum(0));
    leftop1.allocate(big.get_leftBlock()->get_stateInfo());
    operatorfunctions::TensorTrace(dotBlock, *op, big.get_leftBlock(), &(big.get_leftBlock()->get_stateInfo()), leftop1);

    Wavefunction opw2;
    vector<SpinQuantum> dQ = wave1.get_deltaQuantum();
    opw2.initialise(dQ, &big, true);
    operatorfunctions::TensorMultiply(big.get_leftBlock(), leftop1, &big, wave2, opw2, dQ[0], 1.0);
    double sum = sqrt(2.0)*DotProduct(wave1, opw2);

    double difference = 0.0;

    if(dmrginp.spinAdapted()) {
      onepdm(2*ix+1, 2*jx+1) = (sum+difference)/2.0;
      onepdm(2*ix+2, 2*jx+2) = (sum-difference)/2.0;
      onepdm(2*jx+1, 2*ix+1) = (sum+difference)/2.0;
      onepdm(2*jx+2, 2*ix+2) = (sum-difference)/2.0;
    }
    else {
      onepdm(ix+1, jx+1) = sum/sqrt(2.0);
      onepdm(jx+1, ix+1) = sum/sqrt(2.0);
    }
  }      
}

void compute_pair_0_2_0(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, Matrix& onepdm) {
  SpinBlock* leftBlock = big.get_leftBlock()->get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();
  SpinBlock* dotBlock = big.get_leftBlock()->get_rightBlock();

  for (int ij = 0; ij < dotBlock->get_op_array(CRE_CRE).get_size(); ++ij) {
    boost::shared_ptr<SparseMatrix> op = dotBlock->get_op_array(CRE_CRE).get_local_element(ij)[0]->getworkingrepresentation(leftBlock);//spin 0
    int ix = op->get_orbs(0);
    int jx = op->get_orbs(1);
    
    Transposeview top(op);
    Cre leftop1;
    leftop1.set_orbs() = top.get_orbs(); 
    leftop1.set_initialised() = true;
    leftop1.set_fermion() = false;
    leftop1.set_deltaQuantum(1, top.get_deltaQuantum(0));
    leftop1.allocate(big.get_leftBlock()->get_stateInfo());
    operatorfunctions::TensorTrace(dotBlock, top, big.get_leftBlock(), &(big.get_leftBlock()->get_stateInfo()), leftop1);

    Wavefunction opw2;
    vector<SpinQuantum> dQ = wave1.get_deltaQuantum();
    opw2.initialise(dQ, &big, true);
    operatorfunctions::TensorMultiply(big.get_leftBlock(), leftop1, &big, wave2, opw2, dQ[0], 1.0);
    double sum = sqrt(2.0)*DotProduct(wave1, opw2);

    if(dmrginp.spinAdapted()) {
      pout << "BCS with spin adaption not implemented yet." << endl;      
    }
    else {
      onepdm(ix+1, jx+1) = -sum/sqrt(2.0);
      onepdm(jx+1, ix+1) = sum/sqrt(2.0);
    }
  }      
}


void compute_one_pdm_1_1_0(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, Matrix& onepdm)
{
  SpinBlock* leftBlock = big.get_leftBlock()->get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();
  SpinBlock* dotBlock = big.get_leftBlock()->get_rightBlock();


  for (int i = 0; i < leftBlock->get_op_array(CRE).get_size(); ++i)
  for (int j = 0; j < dotBlock->get_op_array(CRE).get_size(); ++j)
  {
    boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_array(CRE).get_local_element(i)[0]->getworkingrepresentation(leftBlock);//spin 0
    boost::shared_ptr<SparseMatrix> op2 = dotBlock->get_op_array(CRE).get_local_element(j)[0]->getworkingrepresentation(dotBlock);//spin 0
    int ix = op1->get_orbs(0);
    int jx = op2->get_orbs(0);

    //now for a leftop that sits on the system block, this is done by tensortrace function
    Cre leftop1;
    leftop1.set_orbs() = op1->get_orbs(); leftop1.set_orbs().push_back(op2->get_orbs(0));
    leftop1.set_initialised() = true;
    leftop1.set_fermion() = false;
    leftop1.set_deltaQuantum(1, (op1->get_deltaQuantum(0)-op2->get_deltaQuantum(0))[0]);
    leftop1.allocate(big.get_leftBlock()->get_stateInfo());
    operatorfunctions::TensorProduct(leftBlock, *op1, Transposeview(op2), big.get_leftBlock(), &(big.get_leftBlock()->get_stateInfo()), leftop1, 1.0);

    Wavefunction opw2;
    vector<SpinQuantum> dQ = wave1.get_deltaQuantum();
    opw2.initialise(dQ, &big, true);
    operatorfunctions::TensorMultiply(big.get_leftBlock(), leftop1, &big, wave2, opw2, dQ[0], 1.0);
    double sum = sqrt(2.0)*DotProduct(wave1, opw2);
    double difference = 0.0;

    if(dmrginp.spinAdapted()) {
      onepdm(2*ix+1, 2*jx+1) = (sum+difference)/2.0;
      onepdm(2*ix+2, 2*jx+2) = (sum-difference)/2.0;
      onepdm(2*jx+1, 2*ix+1) = (sum+difference)/2.0;
      onepdm(2*jx+2, 2*ix+2) = (sum-difference)/2.0;
    }
    else {
      onepdm(ix+1, jx+1) = sum/sqrt(2.0);
      onepdm(jx+1, ix+1) = sum/sqrt(2.0);
    }
  }      
}

void compute_pair_1_1_0(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, Matrix& onepdm) {
  SpinBlock* leftBlock = big.get_leftBlock()->get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();
  SpinBlock* dotBlock = big.get_leftBlock()->get_rightBlock();

  for (int i = 0; i < leftBlock->get_op_array(CRE).get_size(); ++i)
  for (int j = 0; j < dotBlock->get_op_array(CRE).get_size(); ++j)
  {
    boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_array(CRE).get_local_element(i)[0]->getworkingrepresentation(leftBlock);//spin 0
    boost::shared_ptr<SparseMatrix> op2 = dotBlock->get_op_array(CRE).get_local_element(j)[0]->getworkingrepresentation(dotBlock);//spin 0
    int ix = op1->get_orbs(0);
    int jx = op2->get_orbs(0);

    //now for a leftop that sits on the system block, this is done by tensortrace function
    Cre leftop1;
    leftop1.set_orbs() = op1->get_orbs(); leftop1.set_orbs().push_back(op2->get_orbs(0));
    leftop1.set_initialised() = true;
    leftop1.set_fermion() = false;
    leftop1.set_deltaQuantum(1, (-op1->get_deltaQuantum(0)-op2->get_deltaQuantum(0))[0]);
    leftop1.allocate(big.get_leftBlock()->get_stateInfo());
    operatorfunctions::TensorProduct(leftBlock, Transposeview(op1), Transposeview(op2), big.get_leftBlock(), &(big.get_leftBlock()->get_stateInfo()), leftop1, 1.0);

    Wavefunction opw2;
    vector<SpinQuantum> dQ = wave1.get_deltaQuantum();
    opw2.initialise(dQ, &big, true);
    operatorfunctions::TensorMultiply(big.get_leftBlock(), leftop1, &big, wave2, opw2, dQ[0], 1.0);
    double sum = sqrt(2.0)*DotProduct(wave1, opw2);

    if(dmrginp.spinAdapted()) {
      pout << "BCS with spin adaption not implemented yet." << endl;
    }
    else {
      onepdm(ix+1, jx+1) = sum/sqrt(2.0);
      onepdm(jx+1, ix+1) = -sum/sqrt(2.0);
    }
  }
}

void compute_one_pdm_0_2(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, Matrix& onepdm)
{
  if(mpigetrank() == 0) {
  SpinBlock* leftBlock = big.get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();

  for (int ij = 0; ij < rightBlock->get_op_array(CRE_DES).get_size(); ++ij)
  {
    boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_array(CRE_DES).get_local_element(ij)[0]->getworkingrepresentation(rightBlock);
    int ix = op->get_orbs(0);
    int jx = op->get_orbs(1);

    Wavefunction opw2;
    vector<SpinQuantum> dQ = wave1.get_deltaQuantum();
    opw2.initialise(dQ, &big, true);
    operatorfunctions::TensorMultiply(rightBlock, *op, &big, wave2, opw2, dQ[0], 1.0);
    double sum = sqrt(2.0)*DotProduct(wave1, opw2);

    double difference = 0.0;
    
    if(dmrginp.spinAdapted()) {
      onepdm(2*ix+1, 2*jx+1) = (sum+difference)/2.0;
      onepdm(2*ix+2, 2*jx+2) = (sum-difference)/2.0;
      onepdm(2*jx+1, 2*ix+1) = (sum+difference)/2.0;
      onepdm(2*jx+2, 2*ix+2) = (sum-difference)/2.0;
    }
    else {
      onepdm(ix+1, jx+1) = sum/sqrt(2.0);
      onepdm(jx+1, ix+1) = sum/sqrt(2.0);
    }
  }      
  }
}

void compute_pair_0_2(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, Matrix& onepdm)
{
  if(mpigetrank() == 0) {
  SpinBlock* leftBlock = big.get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();

  for (int ij = 0; ij < rightBlock->get_op_array(CRE_CRE).get_size(); ++ij) {
    boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_array(CRE_CRE).get_local_element(ij)[0]->getworkingrepresentation(rightBlock);
    int ix = op->get_orbs(0);
    int jx = op->get_orbs(1);

    Transposeview top(op);
    Wavefunction opw2;
    vector<SpinQuantum> dQ = wave1.get_deltaQuantum();
    opw2.initialise(dQ, &big, true);
    operatorfunctions::TensorMultiply(rightBlock, top, &big, wave2, opw2, dQ[0], 1.0);
    double sum = sqrt(2.0)*DotProduct(wave1, opw2);

    if(dmrginp.spinAdapted()) {
      pout << "BCS with spin adaption not implemented yet." << endl;
    }
    else {
      onepdm(ix+1, jx+1) = -sum/sqrt(2.0);
      onepdm(jx+1, ix+1) = sum/sqrt(2.0);
    }
  }      
  }
}

void compute_one_pdm_1_1(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, Matrix& onepdm)
{
  SpinBlock* leftBlock = big.get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();

  for (int j = 0; j < rightBlock->get_op_array(CRE).get_size(); ++j)
  {
    boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_array(CRE).get_local_element(j)[0]->getworkingrepresentation(rightBlock);
    int jx = op2->get_orbs(0);
    for (int i = 0; i < leftBlock->get_op_array(CRE).get_size(); ++i)
    {
      boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_array(CRE).get_local_element(i)[0]->getworkingrepresentation(leftBlock);
      int ix = op1->get_orbs(0);

      vector<SpinQuantum> opQ = op1->get_deltaQuantum(0)-op2->get_deltaQuantum(0);
      Wavefunction opw2;
      vector<SpinQuantum> dQ = wave1.get_deltaQuantum();
      opw2.initialise(dQ, &big, true);
      operatorfunctions::TensorMultiply(leftBlock, *op1, Transposeview(*op2), &big, wave2, opw2, opQ[0], 1.0);
      double sum = sqrt(2.0)*DotProduct(wave1, opw2);

      double difference = 0.0;

    if(dmrginp.spinAdapted()) {
      onepdm(2*ix+1, 2*jx+1) = (sum+difference)/2.0;
      onepdm(2*ix+2, 2*jx+2) = (sum-difference)/2.0;
      onepdm(2*jx+1, 2*ix+1) = (sum+difference)/2.0;
      onepdm(2*jx+2, 2*ix+2) = (sum-difference)/2.0;
    }
    else {
      onepdm(ix+1, jx+1) = sum/sqrt(2.0);
      onepdm(jx+1, ix+1) = sum/sqrt(2.0);
    }
    }
  }      
}

void compute_pair_1_1(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, Matrix& onepdm) {
  SpinBlock* leftBlock = big.get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();

  for (int j = 0; j < rightBlock->get_op_array(CRE).get_size(); ++j) {
    boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_array(CRE).get_local_element(j)[0]->getworkingrepresentation(rightBlock);
    int jx = op2->get_orbs(0);
    for (int i = 0; i < leftBlock->get_op_array(CRE).get_size(); ++i) {
      boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_array(CRE).get_local_element(i)[0]->getworkingrepresentation(leftBlock);
      int ix = op1->get_orbs(0);

      vector<SpinQuantum> opQ = -op1->get_deltaQuantum(0)-op2->get_deltaQuantum(0);
      Wavefunction opw2;
      vector<SpinQuantum> dQ = wave1.get_deltaQuantum();
      opw2.initialise(dQ, &big, true);
      operatorfunctions::TensorMultiply(leftBlock, Transposeview(op1), Transposeview(op2), &big, wave2, opw2, opQ[0], 1.0);
      double sum = sqrt(2.0)*DotProduct(wave1, opw2);

      if(dmrginp.spinAdapted()) {
        pout << "BCS with spin adaption not implemented yet." << endl;
      }
      else {
        onepdm(ix+1, jx+1) = sum/sqrt(2.0);
        onepdm(jx+1, ix+1) = -sum/sqrt(2.0);
      }
    }
  }
}

void accumulate_onepdm(Matrix& onepdm)
{

#ifndef SERIAL
  Matrix tmp_recv;
  mpi::communicator world;
  if (world.size() == 1)
    return;
  if (!mpigetrank())
    {
      for(int i=1;i<world.size();++i)
	{
	  world.recv(i, i, tmp_recv);
	  for(int k=0;k<onepdm.Nrows();++k)
	    for(int l=0;l<onepdm.Ncols();++l)
	      if(tmp_recv(k+1,l+1) != 0.) {
		onepdm(k+1,l+1) = tmp_recv(k+1,l+1);
	      }
	}
    }
  else
    {
      world.send(0, mpigetrank(), onepdm);
    }
#endif
}


void save_onepdm_text(const Matrix& onepdm, const int &i, const int &j)
{
  //the spatial has a factor of 1/2 in front of it 
  if(!mpigetrank())
  {
    std::vector<int> reorder;
    reorder.resize(onepdm.Nrows()/2);
    for (int k=0; k<onepdm.Nrows()/2; k++) {
      reorder.at(dmrginp.reorder_vector()[k]) = k;
    }


    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/onepdm.", i, j, ".txt");
    ofstream ofs(file);
    ofs << onepdm.Nrows() << endl;
    for(int k=0;k<onepdm.Nrows()/2;++k)
      for(int l=0;l<onepdm.Ncols()/2;++l) {
	int K = reorder.at(k), L = reorder.at(l);

	double opdm = onepdm(2*K+1, 2*L+1) ;
	ofs << boost::format("%d %d %20.14e\n") % (2*k) % (2*l) % opdm;

	opdm = onepdm(2*K+2, 2*L+2);
	ofs << boost::format("%d %d %20.14e\n") % (2*k+1) % (2*l+1) % opdm;
      }

    ofs.close();
  }
}

void save_pairmat_text(const Matrix& pairmat, const int &i, const int &j) {
  if(!mpigetrank()) {
    std::vector<int> reorder;
    reorder.resize(pairmat.Nrows()/2);
    for (int k=0; k<pairmat.Nrows()/2; k++) {
      reorder.at(dmrginp.reorder_vector()[k]) = k;
    }

    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_pairmat.", i, j, ".txt");
    ofstream ofs(file);
    ofs << pairmat.Nrows()/2 << endl;
    for(int k=0;k<pairmat.Nrows()/2;++k)
      for(int l=0;l<pairmat.Ncols()/2;++l) {
	    int K = reorder.at(k), L = reorder.at(l);
	    double pair = pairmat(2*K+1, 2*L+2);
	    ofs << boost::format("%d %d %20.14e\n") % k % l % pair;
      }
    ofs.close();
  }
}

void save_onepdm_spatial_text(const Matrix& onepdm, const int &i, const int &j)
{
  //the spatial has a factor of 1/2 in front of it 
  if(!mpigetrank())
  {
    std::vector<int> reorder;
    reorder.resize(onepdm.Nrows()/2);
    for (int k=0; k<onepdm.Nrows()/2; k++) {
      reorder.at(dmrginp.reorder_vector()[k]) = k;
    }

    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_onepdm.", i, j, ".txt");
    ofstream ofs(file);
    ofs << onepdm.Nrows()/2 << endl;
    for(int k=0;k<onepdm.Nrows()/2;++k)
      for(int l=0;l<onepdm.Ncols()/2;++l) {
	int K = reorder.at(k), L = reorder.at(l);
	double opdm = (onepdm(2*K+1, 2*L+1)+onepdm(2*K+2, 2*L+2));
	ofs << boost::format("%d %d %20.14e\n") % k % l % opdm;
      }

    ofs.close();
  }
}

void save_onepdm_spatial_binary(const Matrix& onepdm, const int &i, const int &j)
{
  //the spatial has a factor of 1/2 in front of it 
  if(!mpigetrank())
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_onepdm_bin.", i, j, ".bin");
    FILE* f = fopen(file, "wb");
    int rows = onepdm.Nrows()/2;
    Matrix spatonepdm(rows, rows);
    for (int i=0; i<rows; i++)
      for (int j=0; j<rows; j++)
	spatonepdm(i+1,j+1) = onepdm(2*i+1, 2*j+1)+onepdm(2*i+2, 2*j+2);
    int result = fwrite(&rows,  1, sizeof(int), f);
    result = fwrite(&rows,  1, sizeof(int), f);
    result = fwrite(spatonepdm.Store(), rows*rows, sizeof(double), f); 
    fclose(f);
  }
}

void save_onepdm_binary(const Matrix& onepdm, const int &i, const int &j)
{
  if(!mpigetrank())
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/onepdm.", i, j);
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << onepdm;
    ofs.close();
  }
}

void load_onepdm_binary(Matrix& onepdm, const int &i, const int &j)
{
  if(!mpigetrank())
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(), "/onepdm.", i, j);
    std::ifstream ifs(file, std::ios::binary);
    boost::archive::binary_iarchive load(ifs);
    load >> onepdm;
    ifs.close();
  }
#ifndef SERIAL
  mpi::communicator world;
  mpi::broadcast(world,onepdm,0);
  if(mpigetrank())
    onepdm = 0;
#endif

}


void save_pairmat_binary(const Matrix& pairmat, const int &i, const int &j)
{
  if(!mpigetrank())
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/pairmat.", i, j);
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << pairmat;
    ofs.close();
  }
}

void load_pairmat_binary(Matrix& pairmat, const int &i, const int &j)
{
  if(!mpigetrank())
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(), "/pairmat.", i, j);
    std::ifstream ifs(file, std::ios::binary);
    boost::archive::binary_iarchive load(ifs);
    load >> pairmat;
    ifs.close();
  }
#ifndef SERIAL
  mpi::communicator world;
  mpi::broadcast(world,pairmat,0);
  if(mpigetrank())
    pairmat = 0;
#endif

}


}
