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

namespace SpinAdapted{

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

double spinExpectation(Wavefunction& wave1, Wavefunction& wave2, SparseMatrix& leftOp, SparseMatrix& dotOp, SparseMatrix& rightOp, const SpinBlock& big)
{

  //calculating <wave1| Oa*Ob | wave2>
  // do transpose specifies if we want  <wave1| Oa^T*Ob |wave2> separately. This can be avoided in some sitations if wave1 and wave2 are the same functions
  int leftindices=0, dotindices=0, rightindices=0;

  leftindices = &leftOp ? leftOp.get_orbs().size() : 0;
  dotindices = &dotOp ? dotOp.get_orbs().size() : 0;
  rightindices = &rightOp ? rightOp.get_orbs().size() : 0;

  int Aindices, Bindices;
  Aindices = leftindices+dotindices;
  Bindices = rightindices;

  Wavefunction opw2;
  SpinQuantum dQ = wave1.get_deltaQuantum();
  opw2.initialise(dQ, &big, true);

  SpinBlock* leftBlock = big.get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();

  Cre AOp; //This is just an example class
  int totalspin = (&rightOp) ? rightOp.get_spin() : 0;

  if (Aindices != 0)
    FormLeftOp(leftBlock, leftOp, dotOp, AOp, totalspin);
  
  //different cases
  if (Aindices == 0 && Bindices == 4)
  {
    operatorfunctions::TensorMultiply(rightBlock, rightOp, &big, wave2, opw2, dQ, 1.0);
//MAW    expectations.push_back( DotProduct(wave1, opw2, dmrginp.Sz(), big) );
  }
  else if (Aindices == 4 && Bindices == 0)
  { 
    operatorfunctions::TensorMultiply(leftBlock, AOp, &big, wave2, opw2, dQ, 1.0);
//MAW    expectations.push_back( DotProduct(wave1, opw2, dmrginp.Sz(), big) );
  }
  else if (Aindices != 0 && Bindices != 0)
  { 
    operatorfunctions::TensorMultiply(leftBlock, AOp, rightOp, &big, wave2, opw2, dQ, 1.0);
//MAW    expectations.push_back( DotProduct(wave1, opw2, dmrginp.Sz(), big) );    
//MAW    if (doTranspose)
//MAW    {
//MAW      opw2.Clear();
//MAW      operatorfunctions::TensorMultiply(leftBlock, Transposeview(AOp), rightOp, &big, wave2, opw2, dQ, 1.0);
//MAW      expectations.push_back( DotProduct(wave1, opw2, dmrginp.Sz(), big) );    
//MAW    }
  }
  return DotProduct(wave1, opw2, dmrginp.Sz(), big);

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void FormLeftOp(const SpinBlock* leftBlock, const SparseMatrix& leftOp, const SparseMatrix& dotOp, SparseMatrix& Aop, int totalspin)
{
  //Cre is just a class..it is not actually cre
  int leftindices=0, dotindices=0, rightindices=0;

  leftindices = &leftOp ? leftOp.get_orbs().size() : 0;
  dotindices = &dotOp ? dotOp.get_orbs().size() : 0;
  
  int Aindices, Bindices;
  Aindices = leftindices+dotindices;

  Aop.CleanUp();
  Aop.set_initialised() = true;
  if (dotindices == 0)
    {
      Aop.set_fermion() = false;
      Aop.set_orbs() = leftOp.get_orbs();
      Aop.set_deltaQuantum() = leftOp.get_deltaQuantum();
      Aop.allocate(leftBlock->get_stateInfo());
      operatorfunctions::TensorTrace(leftBlock->get_leftBlock(), leftOp, leftBlock, &(leftBlock->get_stateInfo()), Aop, 1.0);
    }
  else if (leftindices == 0)
    {
      Aop.set_fermion() = false;
      Aop.set_orbs() = dotOp.get_orbs();
      Aop.set_deltaQuantum() = dotOp.get_deltaQuantum();
      Aop.allocate(leftBlock->get_stateInfo());
      operatorfunctions::TensorTrace(leftBlock->get_rightBlock(), dotOp, leftBlock, &(leftBlock->get_stateInfo()), Aop, 1.0);
    }
  else
    {
      Aop.set_orbs() = leftOp.get_orbs(); copy(dotOp.get_orbs().begin(), dotOp.get_orbs().end(), back_inserter(Aop.set_orbs()));
      Aop.set_fermion() = Aop.set_orbs().size() == 2 ? true : false;
      vector<SpinQuantum> spins = (dotOp.get_deltaQuantum() + leftOp.get_deltaQuantum());
      SpinQuantum dQ;
      for (int i=0; i< spins.size(); i++) {
	if (spins[i].get_s() == totalspin) { dQ = spins[i]; break; }
      }
      Aop.set_deltaQuantum() = dQ;
      Aop.allocate(leftBlock->get_stateInfo());
      operatorfunctions::TensorProduct(leftBlock->get_leftBlock(), leftOp, dotOp, leftBlock, &(leftBlock->get_stateInfo()), Aop, 1.0);      
    }
} 

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void save_twopdm_text(const array_4d<double>& twopdm, const int &i, const int &j)
{
  if(!mpigetrank())
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/twopdm.", i, j);
    ofstream ofs(file);
    ofs << twopdm.dim1() << endl;
    for(int k=0;k<twopdm.dim1();++k)
      for(int l=0;l<twopdm.dim2();++l)
        for(int m=0;m<twopdm.dim3();++m)
          for(int n=0;n<twopdm.dim4();++n)
            ofs << boost::format("%d %d %d %d %20.14e\n") % k % l % m % n % twopdm(k,l,m,n);
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void save_spatial_twopdm_text(const array_4d<double>& twopdm, const int &i, const int &j)
{
  //the spatial has a factor of 1/2 in front of it 
  if(!mpigetrank())
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/spatial_twopdm.", i, j);
    ofstream ofs(file);
    ofs << twopdm.dim1()/2 << endl;
    for(int k=0;k<twopdm.dim1()/2;++k)
      for(int l=0;l<twopdm.dim2()/2;++l)
        for(int m=0;m<twopdm.dim3()/2;++m)
          for(int n=0;n<twopdm.dim4()/2;++n) {
	    double pdm = 0.0;
	    for (int s=0; s<2; s++)
	      for (int t =0; t<2; t++)
		pdm += twopdm(2*k+s, 2*l+t, 2*m+t, 2*n+s)*0.5;
		
            ofs << boost::format("%d %d %d %d %20.14e\n") % k % l % m % n % pdm;
	  }
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void save_spatial_twopdm_binary(const array_4d<double>& twopdm, const int &i, const int &j)
{
  //the spatial has a factor of 1/2 in front of it 
  if(!mpigetrank())
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/spatial_binary_twopdm.", i, j);
    FILE* f = fopen(file, "wb");

    int nrows = twopdm.dim1()/2;
    array_4d<double> pdm; pdm.resize(nrows, nrows, nrows, nrows);
    for(int k=0;k<twopdm.dim1()/2;++k)
      for(int l=0;l<twopdm.dim2()/2;++l)
        for(int m=0;m<twopdm.dim3()/2;++m)
          for(int n=0;n<twopdm.dim4()/2;++n) {
	    pdm(k, l, m, n) = 0.0;
	    for (int s=0; s<2; s++)
	      for (int t =0; t<2; t++)
		pdm(k, l, m, n) += twopdm(2*k+s, 2*l+t, 2*m+t, 2*n+s)*0.5;
	  }
    int result = fwrite(&nrows,  1, sizeof(int), f);
    result = fwrite(&pdm(0,0,0,0), pdm.size(), sizeof(double), f);
    fclose(f);
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void save_twopdm_binary(const array_4d<double>& twopdm, const int &i, const int &j)
{
  if(!mpigetrank())
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/twopdm.", i, j);
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << twopdm;
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void load_twopdm_binary(array_4d<double>& twopdm, const int &i, const int &j)
{
  if(!mpigetrank())
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/twopdm.", i, j);
    std::ifstream ifs(file, std::ios::binary);
    boost::archive::binary_iarchive load(ifs);
    load >> twopdm;
    ifs.close();
  }
#ifndef SERIAL
  mpi::communicator world;
  mpi::broadcast(world,twopdm,0);
  if(mpigetrank())
    twopdm.Clear();
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void save_averaged_twopdm(const int &nroots)
{
  if(!mpigetrank())
  {
    array_4d<double> twopdm;
    char file[5000];
    for(int i=0;i<nroots;++i)
    {
      sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/twopdm.", i, i);
      ifstream ifs(file);
      int size = 0;
      ifs >> size;
      if(i==0)
        twopdm.resize(size,size,size,size);
      int k,l,m,n;
      double val;
      while(ifs >> k)
      {
        ifs >> l >> m >> n >> val;
        twopdm(k,l,m,n) += dmrginp.weights()[i]*val;
      }
    }
    sprintf (file, "%s%s", dmrginp.save_prefix().c_str(),"/twopdm");
    ofstream ofs(file);
    ofs << twopdm.dim1() << endl;
    for(int k=0;k<twopdm.dim1();++k)
      for(int l=0;l<twopdm.dim2();++l)
        for(int m=0;m<twopdm.dim3();++m)
          for(int n=0;n<twopdm.dim4();++n)
            ofs << boost::format("%d %d %d %d %20.14e\n") % k % l % m % n % twopdm(k,l,m,n);

    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void accumulate_twopdm(array_4d<double>& twopdm)
{
#ifndef SERIAL
  array_4d<double> tmp_recv;
  mpi::communicator world;
  if (!mpigetrank())
    {
      for(int i=1;i<world.size();++i)
	{
	  world.recv(i, i, tmp_recv);
	  for(int k=0;k<twopdm.dim1();++k)
	    for(int l=0;l<twopdm.dim2();++l)
	      for(int m=0;m<twopdm.dim3();++m)
		for(int n=0;n<twopdm.dim4();++n)
		  if(tmp_recv(k,l,m,n) != 0.)
		    twopdm(k,l,m,n) = tmp_recv(k,l,m,n);
	}
    }
  else
    {
      world.send(0, mpigetrank(), twopdm);
    }
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void assign_antisymmetric(array_4d<double>& twopdm, const int i, const int j, const int k, const int l, const double val)
{

//MAW
if ( abs(val) > 1e-8 ) pout << "so-twopdm val: i,j,k,l = " << i << "," << j << "," << k << "," << l << "\t\t" << val << endl;
//pout << "so-twopdm val: i,j,k,l = " << i << "," << j << "," << k << "," << l << "\t\t" << val << endl;

  if ( twopdm(i, j, k, l) != 0.0 && abs(twopdm(i,j,k,l)-val) > 1e-6)
    {
      void *array[10];
      size_t size;
      size = backtrace(array, 10);
      cout << "Already calculated "<<i<<" "<<j<<" "<<k<<" "<<l<<endl;
      //backtrace_symbols_fd(array, size, 2);
      cout << "earlier value: "<<twopdm(i,j,k,l)<<endl<< "new value:     "<<val<<endl;
      assert( false );
      return;
    }

  twopdm(i, j, k, l) = val;
  twopdm(i, j, l, k) = -val;
  twopdm(j, i, k, l) = -val;
  twopdm(j, i, l, k) = val;
//  assert ( k != l );
//  assert ( i != j );
//  assert ( (i != j) && (k != l) );
//  if ( k != l ) twopdm(i, j, l, k) = -val;
//  if ( i != j ) twopdm(j, i, k, l) = -val;
//  if ( (i != j) && (k != l) ) twopdm(j, i, l, k) = val;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

double DotProduct(const Wavefunction& w1, const Wavefunction& w2, double Sz, const SpinBlock& big)
{
  int leftOpSz = big.get_leftBlock()->get_stateInfo().quanta.size ();
  int rightOpSz = big.get_rightBlock()->get_stateInfo().quanta.size ();
  const StateInfo* rS = big.get_stateInfo().rightStateInfo, *lS = big.get_stateInfo().leftStateInfo;

  double output = 0.0;
  for (int lQ =0; lQ < leftOpSz; lQ++)
    for (int rQ = 0; rQ < rightOpSz; rQ++) {
      if (w1.allowed(lQ, rQ) && w2.allowed(lQ, rQ))
      {
	double b1b2 = MatrixDotProduct(w1(lQ, rQ), w2(lQ, rQ));
	output += b1b2;
       
      }	
    }

  return output;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

}
