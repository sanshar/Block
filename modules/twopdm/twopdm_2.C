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

void spinExpectation(Wavefunction& wave1, Wavefunction& wave2, SparseMatrix& leftOp, SparseMatrix& dotOp, SparseMatrix& rightOp, const SpinBlock& big, vector<double>& expectations, bool doTranspose)
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
  vector<SpinQuantum> dQ = wave1.get_deltaQuantum();
  opw2.initialise(dQ, &big, true);

  SpinBlock* leftBlock = big.get_leftBlock();
  SpinBlock* rightBlock = big.get_rightBlock();

  Cre AOp; //This is just an example class
  int totalspin = (&rightOp) ? rightOp.get_spin().getirrep() : 0;

  if (Aindices != 0)
    FormLeftOp(leftBlock, leftOp, dotOp, AOp, totalspin);
  
  //different cases
  if (Aindices == 0 && Bindices == 4)
  {
    operatorfunctions::TensorMultiply(rightBlock, rightOp, &big, wave2, opw2, dQ[0], 1.0);
    expectations.push_back( DotProduct(wave1, opw2, dmrginp.Sz(), big) );
  }
  else if (Aindices == 4&& Bindices == 0)
  { 
    operatorfunctions::TensorMultiply(leftBlock, AOp, &big, wave2, opw2, dQ[0], 1.0);
    expectations.push_back( DotProduct(wave1, opw2, dmrginp.Sz(), big) );
  }
  else if (Aindices != 0 && Bindices != 0)
  { 
    operatorfunctions::TensorMultiply(leftBlock, AOp, rightOp, &big, wave2, opw2, dQ[0], 1.0);
    expectations.push_back( DotProduct(wave1, opw2, dmrginp.Sz(), big) );    
    if (doTranspose)
    {
      opw2.Clear();
      operatorfunctions::TensorMultiply(leftBlock, Transposeview(AOp), rightOp, &big, wave2, opw2, dQ[0], 1.0);
      expectations.push_back( DotProduct(wave1, opw2, dmrginp.Sz(), big) );    
    }
  }

}


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
      Aop.set_deltaQuantum(1, leftOp.get_deltaQuantum(0)); // FIXME does leftOp always has only one dQ?
      Aop.allocate(leftBlock->get_stateInfo());
      operatorfunctions::TensorTrace(leftBlock->get_leftBlock(), leftOp, leftBlock, &(leftBlock->get_stateInfo()), Aop, 1.0);
    }
  else if (leftindices == 0)
    {
      Aop.set_fermion() = false;
      Aop.set_orbs() = dotOp.get_orbs();
      Aop.set_deltaQuantum(1, dotOp.get_deltaQuantum(0));
      Aop.allocate(leftBlock->get_stateInfo());
      operatorfunctions::TensorTrace(leftBlock->get_rightBlock(), dotOp, leftBlock, &(leftBlock->get_stateInfo()), Aop, 1.0);
    }
  else
    {
      Aop.set_orbs() = leftOp.get_orbs(); copy(dotOp.get_orbs().begin(), dotOp.get_orbs().end(), back_inserter(Aop.set_orbs()));
      Aop.set_fermion() = Aop.set_orbs().size() == 2 ? true : false;
      vector<SpinQuantum> spins = (dotOp.get_deltaQuantum(0) + leftOp.get_deltaQuantum(0));
      SpinQuantum dQ;
      for (int i=0; i< spins.size(); i++) {
	if (spins[i].get_s().getirrep() == totalspin) { dQ = spins[i]; break; }
      }
      Aop.set_deltaQuantum(1, dQ);
      Aop.allocate(leftBlock->get_stateInfo());
      operatorfunctions::TensorProduct(leftBlock->get_leftBlock(), leftOp, dotOp, leftBlock, &(leftBlock->get_stateInfo()), Aop, 1.0);      
    }
} 


void spin_to_nonspin(vector<int>& indices, vector<double>& coeffs, array_4d<double>& twopdm, Oporder order, bool dotranspose)
{
  //indices are ix jx kx lx 
  // the six possibilities are
  // ix jx+1 kx lx+1     ix jx+1 kx+1 lx     ix+1 jx kx lx+1     ix+1 jx kx+1 lx
  // ix+1 jx+1 kx+1 lx+1   ix jx kx lx

  Matrix a(6,6);
  a  <<  1.0/sqrt(4.0)     <<  -1.0/sqrt(4.0)     <<   -1.0/sqrt(4.0)  <<  1.0/sqrt(4.0)  <<    0.0            <<  0.0           //00 **
     <<  1.0/sqrt(12.0)    <<   1.0/sqrt(12.0)    <<    1.0/sqrt(12.0) <<  1.0/sqrt(12.0) <<    1.0/sqrt(3.0)  <<  1.0/sqrt(3.0) //00 **
     << -1.0/sqrt(4.0)     <<  -1.0/sqrt(4.0)     <<    1.0/sqrt(4.0)  <<  1.0/sqrt(4.0)  <<    0.0            <<  0.0           //20 **
     << -1.0/sqrt(12.0)    <<   1.0/sqrt(12.0)    <<   -1.0/sqrt(12.0) <<  1.0/sqrt(12.0) <<   -1.0/sqrt(3.0)  <<  1.0/sqrt(3.0) //20 **
     <<  1.0/sqrt(6.0)     <<  -1.0/sqrt(6.0)     <<    1.0/sqrt(6.0)  << -1.0/sqrt(6.0)  <<   -1.0/sqrt(6.0)  <<  1.0/sqrt(6.0) //20
     << -1.0/sqrt(6.0)     <<  -1.0/sqrt(6.0)     <<   -1.0/sqrt(6.0)  << -1.0/sqrt(6.0)  <<    1.0/sqrt(6.0)  <<  1.0/sqrt(6.0); //40 

  ColumnVector x(6), b(6);

  b = 0.0; x=0.0;
  b(1) = coeffs[0]; b(2) = coeffs[1];

  int ix = 2*indices[0], jx = 2*indices[1], kx = 2*indices[2], lx = 2*indices[3];

  if (order == CC_D_D)
  {
    //do nothing same as default
  }

  if (order == CD_D_C)
  {
    a.Row(1)  << -1.0/sqrt(4.0)     <<   0.0               <<    0.0            << -1.0/2.0        <<   -1.0/2.0        << -1.0/2.0 ;      //00 **
    a.Row(2) <<  1.0/sqrt(12.0)    <<  -1.0/sqrt(3.0)     <<   -1.0/sqrt(3.0)  <<  1.0/sqrt(12.0) <<    -1.0/sqrt(12.0) << -1.0/sqrt(12.0) ;//00 **
  }

  if (order == CD_CD)
  {
    a.Row(1)  << -1.0/sqrt(4.0)     <<   0.0               <<    0.0            << -1.0/2.0        <<   -1.0/2.0        << -1.0/2.0 ;      //00 **
    a.Row(2) << -1.0/sqrt(12.0)    <<   1.0/sqrt(3.0)     <<    1.0/sqrt(3.0)  << -1.0/sqrt(12.0) <<    1.0/sqrt(12.0) <<  1.0/sqrt(12.0) ;//00 **
  }

  if (order == D_CD_C)
  {
    a.Row(1)  << 0.0     <<    1.0/2.0         <<    1.0/2.0        <<   0.0          <<   1.0/2.0        <<  1.0/2.0 ;      //00 **
    a.Row(2) << -1.0/sqrt(3.0)    <<  1.0/sqrt(12.0)     <<   1.0/sqrt(12.0)  << -1.0/sqrt(3.0) <<  -1.0/sqrt(12.0) << -1.0/sqrt(12.0) ;//00 **
  }
  
  if (order == C_CD_D)
  {
    a.Row(1)  << 0.0     <<    1.0/2.0         <<    1.0/2.0        <<   0.0          <<   1.0/2.0        <<  1.0/2.0 ;      //00 **
    a.Row(2) << 1.0/sqrt(3.0)    <<  -1.0/sqrt(12.0)     <<   -1.0/sqrt(12.0)  <<  1.0/sqrt(3.0) <<   1.0/sqrt(12.0) <<  1.0/sqrt(12.0) ;//00 **
  }

  if (order == D_CC_D)
  {
    a.Row(1)  << 1.0/2.0    <<   -1.0/2.0         <<    -1.0/2.0        << 1.0/2.0        <<   0.0        <<  0.0 ;      //00 **
    a.Row(2) << -1.0/sqrt(12.0)    <<  -1.0/sqrt(12.0)     <<  -1.0/sqrt(12.0)  << -1.0/sqrt(12.0) <<-1.0/sqrt(3.0) << -1.0/sqrt(3.0) ;//00 **
  }

  xsolve_AxeqB(a, b, x);

  assign_antisymmetric(twopdm, ix  , jx+1, kx  , lx+1, x(1));
  assign_antisymmetric(twopdm, ix  , jx+1, kx+1, lx  , x(2));
  assign_antisymmetric(twopdm, ix+1, jx  , kx  , lx+1, x(3));
  assign_antisymmetric(twopdm, ix+1, jx  , kx+1, lx  , x(4));
  assign_antisymmetric(twopdm, ix+1, jx+1, kx+1, lx+1, x(5));
  assign_antisymmetric(twopdm, ix  , jx  , kx  , lx  , x(6));

  if (dotranspose)
  {
    int temp = ix; ix = lx; lx = temp; temp = kx; kx = jx; jx = temp;  
    assign_antisymmetric(twopdm, ix  , jx+1, kx  , lx+1, x(1));
    assign_antisymmetric(twopdm, ix  , jx+1, kx+1, lx  , x(2));
    assign_antisymmetric(twopdm, ix+1, jx  , kx  , lx+1, x(3));
    assign_antisymmetric(twopdm, ix+1, jx  , kx+1, lx  , x(4));
    assign_antisymmetric(twopdm, ix+1, jx+1, kx+1, lx+1, x(5));
    assign_antisymmetric(twopdm, ix  , jx  , kx  , lx  , x(6));
  }

}


void save_twopdm_text(const array_4d<double>& twopdm, const int &i, const int &j)
{
  if(!mpigetrank())
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/twopdm.", i, j, ".txt");
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

void save_spatial_twopdm_text(const array_4d<double>& twopdm, const int &i, const int &j)
{
  //the spatial has a factor of 1/2 in front of it 
  if(!mpigetrank())
  {
    std::vector<int> reorder;
    reorder.resize(twopdm.dim1()/2);
    for (int k=0; k<twopdm.dim2()/2; k++) {
      reorder.at(dmrginp.reorder_vector()[k]) = k;
    }

    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_twopdm.", i, j, ".txt");
    ofstream ofs(file);
    ofs << twopdm.dim1()/2 << endl;
    for(int k=0;k<twopdm.dim1()/2;++k)
      for(int l=0;l<twopdm.dim2()/2;++l)
        for(int m=0;m<twopdm.dim3()/2;++m)
          for(int n=0;n<twopdm.dim4()/2;++n) {
	    double pdm = 0.0;
	    for (int s=0; s<2; s++)
	      for (int t =0; t<2; t++)
		pdm += twopdm(2*reorder.at(k)+s, 2*reorder.at(l)+t, 2*reorder.at(m)+t, 2*reorder.at(n)+s)*0.5;
		
            ofs << boost::format("%d %d %d %d %20.14e\n") % k % l % m % n % pdm;
	  }
    ofs.close();
  }
}

void save_spatial_twopdm_binary(const array_4d<double>& twopdm, const int &i, const int &j)
{
  //the spatial has a factor of 1/2 in front of it 
  if(!mpigetrank())
  {
    std::vector<int> reorder;
    reorder.resize(twopdm.dim1()/2);
    for (int k=0; k<twopdm.dim2()/2; k++) {
      reorder.at(dmrginp.reorder_vector()[k]) = k;
    }

    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_binary_twopdm.", i, j, ".bin");
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
		pdm(k,l,m,n) += twopdm(2*reorder.at(k)+s, 2*reorder.at(l)+t, 2*reorder.at(m)+t, 2*reorder.at(n)+s)*0.5;
	  }
    int result = fwrite(&nrows,  1, sizeof(int), f);
    result = fwrite(&pdm(0,0,0,0), pdm.size(), sizeof(double), f);
    fclose(f);
  }
}

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


void assign_antisymmetric(array_4d<double>& twopdm, const int i, const int j, const int k, const int l, const double val)
{

  if ( twopdm(i, j, k, l) != 0.0 && (twopdm(i,j,k,l)-val) > 2e-4)
    {
      void *array[10];
      size_t size;
      size = backtrace(array, 10);
      pout << "Already calculated "<<i<<" "<<j<<" "<<k<<" "<<l<<endl;
      //backtrace_symbols_fd(array, size, 2);
      pout << "earlier value: "<<twopdm(i,j,k,l)<<endl<<"new value: "<<val<<endl;
      assert(1 == 0);
    }

  twopdm(i, j, k, l) = val;
  twopdm(i, j, l, k) = -val;
  twopdm(j, i, k, l) = -val;
  twopdm(j, i, l, k) = val;
}

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
}
