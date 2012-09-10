/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is free software: you can redistribute it and/or modify         
it under the terms of the GNU General Public License as published by         
the Free Software Foundation, either version 3 of the License, or            
(at your option) any later version.                                          
                                                                             
This program is distributed in the hope that it will be useful,              
but WITHOUT ANY WARRANTY; without even the implied warranty of               
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
GNU General Public License for more details.                                 
                                                                             
You should have received a copy of the GNU General Public License            
along with this program.  If not, see <http://www.gnu.org/licenses/>.        

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

	pout <<"Performing sweep calculation "<<endl;

	//if (big.get_leftBlock()->size() == 2) {
	cout << "compute 2_0 "<<mpigetrank()<<endl;
	  compute_one_pdm_2_0(wavefunction1, wavefunction2, big, onepdm);
	  //}

	  //if (big.get_rightBlock()->size() == 1) {
	  cout << "compute 0_2 "<<mpigetrank()<<endl;
	  compute_one_pdm_0_2(wavefunction1, wavefunction2, big, onepdm);
	  //}

	  cout << "compute 1_1 "<<mpigetrank()<<endl;
	compute_one_pdm_1_1(wavefunction1, wavefunction2, big, onepdm);
	accumulate_onepdm(onepdm);
	save_onepdm_binary(onepdm, i, j);
      }
}

void compute_one_pdm_2_0(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, Matrix& onepdm)
{
  SpinBlock* leftBlock = big.get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();

  for (int ij = 0; ij < leftBlock->get_op_array(CRE_DES).get_size(); ++ij)
  {
    boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_array(CRE_DES).get_local_element(ij)[0]->getworkingrepresentation(leftBlock);//spin 0
    int ix = op->get_orbs(0);
    int jx = op->get_orbs(1);

    Wavefunction opw2;
    SpinQuantum dQ = wave1.get_deltaQuantum();
    opw2.initialise(dQ, &big, true);
    operatorfunctions::TensorMultiply(leftBlock, *op, &big, wave2, opw2, dQ, 1.0);
    double sum = sqrt(2.0)*DotProduct(wave1, opw2);

    opw2.Clear();
    op = leftBlock->get_op_array(CRE_DES).get_local_element(ij)[1]->getworkingrepresentation(leftBlock);//spin 2
    operatorfunctions::TensorMultiply(leftBlock, *op, &big, wave2, opw2, dQ, 1.0);
    double difference = sqrt(2.0)*DotProduct(wave1, opw2)*cg(dQ.get_s(), dmrginp.Sz(), 2, 0, dQ.get_s(), dmrginp.Sz());

    onepdm(2*ix+1, 2*jx+1) = (sum+difference)/2.0;
    onepdm(2*ix+2, 2*jx+2) = (sum-difference)/2.0;

    if (ix != jx) {
      opw2.Clear();
      op = leftBlock->get_op_array(CRE_DES).get_local_element(ij)[0]->getworkingrepresentation(leftBlock);
      operatorfunctions::TensorMultiply(leftBlock, Transposeview(*op), &big, wave2, opw2, dQ, 1.0);
      sum = sqrt(2.0)*DotProduct(wave1, opw2);

      opw2.Clear();
      op = leftBlock->get_op_array(CRE_DES).get_local_element(ij)[1]->getworkingrepresentation(leftBlock);
      operatorfunctions::TensorMultiply(leftBlock, Transposeview(*op), &big, wave2, opw2, dQ, 1.0);
      difference = sqrt(2.0)*DotProduct(wave1, opw2)*cg(dQ.get_s(), dmrginp.Sz(), 2, 0, dQ.get_s(), dmrginp.Sz());
      
      onepdm(2*jx+1, 2*ix+1) = (sum+difference)/2.0;
      onepdm(2*jx+2, 2*ix+2) = (sum-difference)/2.0;
    }

  }      
}

void compute_one_pdm_0_2(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, Matrix& onepdm)
{
  SpinBlock* leftBlock = big.get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();

  for (int ij = 0; ij < rightBlock->get_op_array(CRE_DES).get_size(); ++ij)
  {
    boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_array(CRE_DES).get_local_element(ij)[0]->getworkingrepresentation(rightBlock);
    int ix = op->get_orbs(0);
    int jx = op->get_orbs(1);

    Wavefunction opw2;
    SpinQuantum dQ = wave1.get_deltaQuantum();
    opw2.initialise(dQ, &big, true);
    operatorfunctions::TensorMultiply(rightBlock, *op, &big, wave2, opw2, dQ, 1.0);
    double sum = sqrt(2.0)*DotProduct(wave1, opw2);
    opw2.Clear();
    op = rightBlock->get_op_array(CRE_DES).get_local_element(ij)[1]->getworkingrepresentation(rightBlock);
    operatorfunctions::TensorMultiply(rightBlock, *op, &big, wave2, opw2, dQ, 1.0);
    double difference = sqrt(2.0)*DotProduct(wave1, opw2)*cg(dQ.get_s(), dmrginp.Sz(), 2, 0, dQ.get_s(), dmrginp.Sz());
    
    onepdm(2*ix+1, 2*jx+1) = (sum+difference)/2.0;
    onepdm(2*ix+2, 2*jx+2) = (sum-difference)/2.0;

    if (ix != jx) {
      opw2.Clear();
      op = rightBlock->get_op_array(CRE_DES).get_local_element(ij)[0]->getworkingrepresentation(rightBlock);
      operatorfunctions::TensorMultiply(rightBlock, Transposeview(*op), &big, wave2, opw2, dQ, 1.0);
      sum = sqrt(2.0)*DotProduct(wave1, opw2);

      opw2.Clear();
      op = rightBlock->get_op_array(CRE_DES).get_local_element(ij)[1]->getworkingrepresentation(rightBlock);
      operatorfunctions::TensorMultiply(rightBlock, Transposeview(*op), &big, wave2, opw2, dQ, 1.0);
      double difference = sqrt(2.0)*DotProduct(wave1, opw2)*cg(dQ.get_s(), dmrginp.Sz(), 2, 0, dQ.get_s(), dmrginp.Sz());
    
      onepdm(2*jx+1, 2*ix+1) = (sum+difference)/2.0;
      onepdm(2*jx+2, 2*ix+2) = (sum-difference)/2.0;
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

      vector<SpinQuantum> opQ = op1->get_deltaQuantum()+op2->get_deltaQuantum();
      Wavefunction opw2;
      SpinQuantum dQ = wave1.get_deltaQuantum();
      opw2.initialise(dQ, &big, true);
      operatorfunctions::TensorMultiply(leftBlock, *op1, Transposeview(*op2), &big, wave2, opw2, opQ[0], 1.0);
      double sum = sqrt(2.0)*DotProduct(wave1, opw2);

      opw2.Clear();
      operatorfunctions::TensorMultiply(leftBlock, *op1, Transposeview(*op2), &big, wave2, opw2, opQ[1], 1.0);
      double difference = sqrt(2.0)*DotProduct(wave1, opw2)*cg(dQ.get_s(), dmrginp.Sz(), 2, 0, dQ.get_s(), dmrginp.Sz());

      onepdm(2*ix+1, 2*jx+1) = (sum+difference)/2.0;
      onepdm(2*ix+2, 2*jx+2) = (sum-difference)/2.0;

      opw2.Clear();
      operatorfunctions::TensorMultiply(leftBlock, Transposeview(*op1), *op2, &big, wave2, opw2, opQ[0], 1.0);
      sum = sqrt(2.0)*DotProduct(wave1, opw2);
      opw2.Clear();
      operatorfunctions::TensorMultiply(leftBlock, Transposeview(*op1), *op2, &big, wave2, opw2, opQ[1], 1.0);
      difference = sqrt(2.0)*DotProduct(wave1, opw2)*cg(dQ.get_s(), dmrginp.Sz(), 2, 0, dQ.get_s(), dmrginp.Sz());
    
      
      onepdm(2*jx+1, 2*ix+1) = (sum+difference)/2.0;
      onepdm(2*jx+2, 2*ix+2) = (sum-difference)/2.0;
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
	      if(tmp_recv(k+1,l+1) != 0.)
		onepdm(k+1,l+1) = tmp_recv(k+1,l+1);
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
    char file[5000];
    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/onepdm.", i, j);
    ofstream ofs(file);
    ofs << onepdm.Nrows() << endl;
    for(int k=0;k<onepdm.Nrows();++k)
      for(int l=0;l<onepdm.Ncols();++l)
	ofs << boost::format("%d %d %20.14e\n") % k % l % onepdm(k+1, l+1);

    ofs.close();
  }
}

void save_onepdm_spatial_text(const Matrix& onepdm, const int &i, const int &j)
{
  //the spatial has a factor of 1/2 in front of it 
  if(!mpigetrank())
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/spatial_onepdm.", i, j);
    ofstream ofs(file);
    ofs << onepdm.Nrows()/2 << endl;
    for(int k=0;k<onepdm.Nrows()/2;++k)
      for(int l=0;l<onepdm.Ncols()/2;++l) {
	double opdm = (onepdm(2*k+1, 2*l+1)+onepdm(2*k+2, 2*l+2));
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
    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/spatial_onepdm_bin.", i, j);
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
}
