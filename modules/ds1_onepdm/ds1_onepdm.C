/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "ds1_onepdm.h"
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
namespace ds1_onepdm{	

void compute_one_pdm_2_0(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, Matrix& onepdm)
{
  SpinBlock* leftBlock = big.get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();

  int ketS = dmrginp.total_spin_number().getirrep();
  int braS = dmrginp.bra_spin_number().getirrep();
    pout << " ketS " << ketS  << "   braS   "<< braS  <<endl;


  boost::shared_ptr<SparseMatrix> overlap = rightBlock->get_op_array(OVERLAP).get_local_element(0)[0]->getworkingrepresentation(rightBlock);//EL
  for (int ij = 0; ij < leftBlock->get_op_array(CRE_DES).get_size(); ++ij)
  {
    boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_array(CRE_DES).get_local_element(ij)[1]->getworkingrepresentation(leftBlock);//spin 0
    int ix = op->get_orbs(0);
    int jx = op->get_orbs(1);

    SpinQuantum opQ = op->get_deltaQuantum(0);

    Wavefunction opw2;
    vector<SpinQuantum> dQ = wave1.get_deltaQuantum();
    opw2.initialisebra(dQ, &big, true);

    operatorfunctions::TensorMultiply(leftBlock, *op, *overlap, &big, wave2, opw2, opQ, 1.0);
    double sum = DotProduct(wave1, opw2);

    pout << " left CREDES "  <<endl;
if ( ketS >= braS){
      onepdm(2*jx+1, 2*ix+1) = sum/2.0;
      onepdm(2*jx+2, 2*ix+2) = sum/2.0;
    pout << "onepdm(2*jx+1, 2*ix+1)   "<<  "ix   "<< ix<<"  jx    "<< jx  <<"   "<< 2*jx+1 <<"   "<< 2*ix+1<<"  "<<onepdm(2*jx+1, 2*ix+1)<<endl;
    pout << "onepdm(2*jx+2, 2*ix+2)   "<<  "ix   "<< ix<<"  jx    "<< jx  <<"   "<< 2*jx+2 <<"   "<< 2*ix+2<<"  "<<onepdm(2*jx+2, 2*ix+2)<<endl;

      }
if (braS > ketS){
      onepdm(2*ix+1, 2*jx+1) = sum/2.0;
      onepdm(2*ix+2, 2*jx+2) = sum/2.0;
    pout << "onepdm(2*ix+1, 2*jx+1)   "<<  "ix   "<< ix<<"  jx    "<< jx  <<"   "<< 2*ix+1 <<"   "<< 2*jx+1<<"  "<<onepdm(2*ix+1, 2*jx+1)<<endl;
    pout << "onepdm(2*ix+2, 2*jx+2)   "<<  "ix   "<< ix<<"  jx    "<< jx  <<"   "<< 2*ix+2 <<"   "<< 2*jx+2<<"  "<<onepdm(2*ix+2, 2*jx+2)<<endl;
      }
}

  for (int ij = 0; ij < leftBlock->get_op_array(DES_CRE).get_size(); ++ij)
  {
    boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_array(DES_CRE).get_local_element(ij)[1]->getworkingrepresentation(leftBlock);//spin 0
    int ix = op->get_orbs(0);
    int jx = op->get_orbs(1);

    if (ix!=jx) {

    SpinQuantum opQ = op->get_deltaQuantum(0);

    Wavefunction opw2;
    vector<SpinQuantum> dQ = wave1.get_deltaQuantum();
    opw2.initialisebra(dQ, &big, true);
    operatorfunctions::TensorMultiply(leftBlock, *op, *overlap, &big, wave2, opw2, opQ, 1.0);

    double sum = DotProduct(wave1, opw2);

    pout << " left DESCRE "  <<endl;

if ( ketS >= braS){
     onepdm(2*ix+1, 2*jx+1) = -sum/2.0;
     onepdm(2*ix+2, 2*jx+2) = -sum/2.0;
    pout << "onepdm(2*ix+1, 2*jx+1)   "<<  "ix   "<< ix<<"  jx    "<< jx  <<"   "<< 2*ix+1 <<"   "<< 2*jx+1<<"  "<<onepdm(2*ix+1, 2*jx+1)<<endl;
    pout << "onepdm(2*ix+2, 2*jx+2)   "<<  "ix   "<< ix<<"  jx    "<< jx  <<"   "<< 2*ix+2 <<"   "<< 2*jx+2<<"  "<<onepdm(2*ix+2, 2*jx+2)<<endl;
     }
if (braS > ketS){
     onepdm(2*jx+1, 2*ix+1) = -sum/2.0;
     onepdm(2*jx+2, 2*ix+2) = -sum/2.0;
    pout << "onepdm(2*jx+1, 2*ix+1)   "<<  "ix   "<< ix<<"  jx    "<< jx  <<"   "<< 2*jx+1 <<"   "<< 2*ix+1<<"  "<<onepdm(2*jx+1, 2*ix+1)<<endl;
    pout << "onepdm(2*jx+2, 2*ix+2)   "<<  "ix   "<< ix<<"  jx    "<< jx  <<"   "<< 2*jx+2 <<"   "<< 2*ix+2<<"  "<<onepdm(2*jx+2, 2*ix+2)<<endl;
     }     
   }
 }
}

void compute_one_pdm_0_2(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, Matrix& onepdm)
{

  SpinBlock* leftBlock = big.get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();
  int ketS = dmrginp.total_spin_number().getirrep();
  int braS = dmrginp.bra_spin_number().getirrep();

  boost::shared_ptr<SparseMatrix> overlap = leftBlock->get_op_array(OVERLAP).get_local_element(0)[0]->getworkingrepresentation(leftBlock);//EL
  for (int ij = 0; ij < rightBlock->get_op_array(CRE_DES).get_size(); ++ij)
  {
    boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_array(CRE_DES).get_local_element(ij)[1]->getworkingrepresentation(rightBlock);
    int ix = op->get_orbs(0);
    int jx = op->get_orbs(1);

    SpinQuantum opQ = op->get_deltaQuantum(0);

    Wavefunction opw2;
    vector<SpinQuantum> dQ = wave1.get_deltaQuantum();
    opw2.initialisebra(dQ, &big, true);

    operatorfunctions::TensorMultiply(rightBlock,*op, *overlap, &big, wave2, opw2, opQ, 1.0);  //EL
    double sum = DotProduct(wave1, opw2);

    pout << " right CREDES "  <<endl;

if ( ketS >= braS){
      onepdm(2*ix+1, 2*jx+1) = sum/2.0;
      onepdm(2*ix+2, 2*jx+2) = sum/2.0;
    pout << "onepdm(2*ix+1, 2*jx+1)   "<<  "ix   "<< ix<<"  jx    "<< jx  <<"   "<< 2*ix+1 <<"   "<< 2*jx+1<<"  "<<onepdm(2*ix+1, 2*jx+1)<<endl;
    pout << "onepdm(2*ix+2, 2*jx+2)   "<<  "ix   "<< ix<<"  jx    "<< jx  <<"   "<< 2*ix+2 <<"   "<< 2*jx+2<<"  "<<onepdm(2*ix+2, 2*jx+2)<<endl;
      }
if ( braS > ketS){
      onepdm(2*jx+1, 2*ix+1) = sum/2.0;
      onepdm(2*jx+2, 2*ix+2) = sum/2.0;
    pout << "onepdm(2*jx+1, 2*ix+1)   "<<  "ix   "<< ix<<"  jx    "<< jx  <<"   "<< 2*jx+1 <<"   "<< 2*ix+1<<"  "<<onepdm(2*jx+1, 2*ix+1)<<endl;
    pout << "onepdm(2*jx+2, 2*ix+2)   "<<  "ix   "<< ix<<"  jx    "<< jx  <<"   "<< 2*jx+2 <<"   "<< 2*ix+2<<"  "<<onepdm(2*jx+2, 2*ix+2)<<endl;
      }
}
//---
   for (int ij = 0; ij < rightBlock->get_op_array(DES_CRE).get_size(); ++ij)
  {
    boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_array(DES_CRE).get_local_element(ij)[1]->getworkingrepresentation(rightBlock);
    int ix = op->get_orbs(0);
    int jx = op->get_orbs(1);

    if (ix!=jx) {
    SpinQuantum opQ = op->get_deltaQuantum(0);

    Wavefunction opw2;
    vector<SpinQuantum> dQ = wave1.get_deltaQuantum();
    opw2.initialisebra(dQ, &big, true);

    operatorfunctions::TensorMultiply(rightBlock,*op, *overlap, &big, wave2, opw2, opQ, 1.0);  //EL
    double sum = DotProduct(wave1, opw2);

    pout << " right DESCRE "  <<endl;

if ( ketS >= braS){
      onepdm(2*jx+1, 2*ix+1) = -sum/2.0;
      onepdm(2*jx+2, 2*ix+2) = -sum/2.0;
    pout << "onepdm(2*jx+1, 2*ix+1)   "<<  "ix   "<< ix<<"  jx    "<< jx  <<"   "<< 2*jx+1 <<"   "<< 2*ix+1<<"  "<<onepdm(2*jx+1, 2*ix+1)<<endl;
    pout << "onepdm(2*jx+2, 2*ix+2)   "<<  "ix   "<< ix<<"  jx    "<< jx  <<"   "<< 2*jx+2 <<"   "<< 2*ix+2<<"  "<<onepdm(2*jx+2, 2*ix+2)<<endl;
      }
if ( braS > ketS){
      onepdm(2*ix+1, 2*jx+1) = -sum/2.0;
      onepdm(2*ix+2, 2*jx+2) = -sum/2.0;
    pout << "onepdm(2*ix+1, 2*jx+1)   "<<  "ix   "<< ix<<"  jx    "<< jx  <<"   "<< 2*ix+1 <<"   "<< 2*jx+1<<"  "<<onepdm(2*ix+1, 2*jx+1)<<endl;
    pout << "onepdm(2*ix+2, 2*jx+2)   "<<  "ix   "<< ix<<"  jx    "<< jx  <<"   "<< 2*ix+2 <<"   "<< 2*jx+2<<"  "<<onepdm(2*ix+2, 2*jx+2)<<endl;
   }
  }
 }
}

void compute_one_pdm_1_1(Wavefunction& wave1, Wavefunction& wave2, const SpinBlock& big, Matrix& onepdm)
{
  SpinBlock* leftBlock = big.get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();

  int ketS = dmrginp.total_spin_number().getirrep();
  int braS = dmrginp.bra_spin_number().getirrep();

  for (int j = 0; j < rightBlock->get_op_array(CRE).get_size(); ++j)	
   {
    boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_array(CRE).get_local_element(j)[0]->getworkingrepresentation(rightBlock);
    int jx = op2->get_orbs(0);
    for (int i = 0; i < leftBlock->get_op_array(DES).get_size(); ++i)
    {
      boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_array(DES).get_local_element(i)[0]->getworkingrepresentation(leftBlock);
      int ix = op1->get_orbs(0);

      vector<SpinQuantum> opQ = op2->get_deltaQuantum(0)+op1->get_deltaQuantum(0);
      Wavefunction opw2;
      vector<SpinQuantum> dQ = wave1.get_deltaQuantum();
      opw2.initialisebra(dQ, &big, true);

      operatorfunctions::TensorMultiply(rightBlock, *op2, *op1, &big, wave2, opw2, opQ[1], 1.0);
      double sum = DotProduct(wave1, opw2);

    pout << " right  CRE  left DES "  <<endl;

if ( ketS >= braS){
      onepdm(2*ix+1, 2*jx+1) = -sum/2.0;
      onepdm(2*ix+2, 2*jx+2) = -sum/2.0;
    pout << "onepdm(2*ix+1, 2*jx+1)   "<<  "ix   "<< ix<<"  jx    "<< jx  <<"   "<< 2*ix+1 <<"   "<< 2*jx+1<<"  "<<onepdm(2*ix+1, 2*jx+1)<<endl;
    pout << "onepdm(2*ix+2, 2*jx+2)   "<<  "ix   "<< ix<<"  jx    "<< jx  <<"   "<< 2*ix+2 <<"   "<< 2*jx+2<<"  "<<onepdm(2*ix+2, 2*jx+2)<<endl;
      }
if ( braS > ketS){
      onepdm(2*jx+1, 2*ix+1) = -sum/2.0;
      onepdm(2*jx+2, 2*ix+2) = -sum/2.0;
    pout << "onepdm(2*jx+1, 2*ix+1)   "<<  "ix   "<< ix<<"  jx    "<< jx  <<"   "<< 2*jx+1 <<"   "<< 2*ix+1<<"  "<<onepdm(2*jx+1, 2*ix+1)<<endl;
    pout << "onepdm(2*jx+2, 2*ix+2)   "<<  "ix   "<< ix<<"  jx    "<< jx  <<"   "<< 2*jx+2 <<"   "<< 2*ix+2<<"  "<<onepdm(2*jx+2, 2*ix+2)<<endl;
      }
    }
  }

//-----
  for (int j = 0; j < rightBlock->get_op_array(DES).get_size(); ++j)
  {
    boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_array(DES).get_local_element(j)[0]->getworkingrepresentation(rightBlock);
    int jx = op2->get_orbs(0);
    for (int i = 0; i < leftBlock->get_op_array(CRE).get_size(); ++i)
    {
      boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_array(CRE).get_local_element(i)[0]->getworkingrepresentation(leftBlock);
      int ix = op1->get_orbs(0);

      vector<SpinQuantum> opQ = op1->get_deltaQuantum(0)+op2->get_deltaQuantum(0);
      Wavefunction opw2;
      vector<SpinQuantum> dQ = wave1.get_deltaQuantum();
      opw2.initialisebra(dQ, &big, true);
      operatorfunctions::TensorMultiply(leftBlock, *op1, *op2, &big, wave2, opw2, opQ[1], 1.0);
      double sum = DotProduct(wave1, opw2);

    pout << " left  CRE  right DES "  <<endl;

if ( ketS >= braS){
      onepdm(2*jx+1, 2*ix+1) = sum/2.0;
      onepdm(2*jx+2, 2*ix+2) = sum/2.0;
    pout << "onepdm(2*jx+1, 2*ix+1)   "<<  "ix   "<< ix<<"  jx    "<< jx  <<"   "<< 2*jx+1 <<"   "<< 2*ix+1<<"  "<<onepdm(2*jx+1, 2*ix+1)<<endl;
    pout << "onepdm(2*jx+2, 2*ix+2)   "<<  "ix   "<< ix<<"  jx    "<< jx  <<"   "<< 2*jx+2 <<"   "<< 2*ix+2<<"  "<<onepdm(2*jx+2, 2*ix+2)<<endl;
     }
if ( braS > ketS){
      onepdm(2*ix+1, 2*jx+1) = sum/2.0;
      onepdm(2*ix+2, 2*jx+2) = sum/2.0;
    pout << "onepdm(2*ix+1, 2*jx+1)   "<<  "ix   "<< ix<<"  jx    "<< jx  <<"   "<< 2*ix+1 <<"   "<< 2*jx+1<<"  "<<onepdm(2*ix+1, 2*jx+1)<<endl;
    pout << "onepdm(2*ix+2, 2*jx+2)   "<<  "ix   "<< ix<<"  jx    "<< jx  <<"   "<< 2*ix+2 <<"   "<< 2*jx+2<<"  "<<onepdm(2*ix+2, 2*jx+2)<<endl;
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
    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/onepdm.", i, j);
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
    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/spatial_onepdm.", i, j);
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
}
