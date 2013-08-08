/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include <boost/format.hpp>
#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include "execinfo.h"

#include "twopdm.h"
#include "npdm_patterns.h"
#include "npdm_expectations.h"
#include "npdm_operator_wrappers.h"

namespace SpinAdapted{

// Forward declaration
boost::shared_ptr<NpdmSpinOps> select_op_wrapper( SpinBlock * spinBlock, std::vector<Npdm::CD> & cd_type );

//===========================================================================================================================================================

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
      	      for (int t =0; t<2; t++) {
//                 if ( (k==0) && (l==0) && (m==1) && (n==2) ) 
//                 if ( (k==0) && (l==0) && (m==2) && (n==1) ) 
//                     std::cout << 2*k+s<<","<< 2*l+t<<","<< 2*m+t<<","<< 2*n+s<<"\t\t"<<pdm<<"\t"<<twopdm(2*k+s, 2*l+t, 2*m+t, 2*n+s)*0.5 <<std::endl;
                 pdm += twopdm(2*k+s, 2*l+t, 2*m+t, 2*n+s)*0.5;
               }
		
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

void assign_twopdm_antisymmetric(array_4d<double>& twopdm, const int i, const int j, const int k, const int l, const double val)
{

//MAW
//if ( abs(val) > 1e-8 ) pout << "so-twopdm val: i,j,k,l = " << i << "," << j << "," << k << "," << l << "\t\t" << val << endl;
//pout << "so-twopdm val: i,j,k,l = " << i << "," << j << "," << k << "," << l << "\t\t" << val << endl;

  // Test for duplicates
  if ( twopdm(i, j, k, l) != 0.0 && abs(twopdm(i,j,k,l)-val) > 1e-6)
    {
      void *array[10];
      size_t size;
      size = backtrace(array, 10);
      cout << "WARNING: Already calculated "<<i<<" "<<j<<" "<<k<<" "<<l<<endl;
      //backtrace_symbols_fd(array, size, 2);
      cout << "earlier value: "<<twopdm(i,j,k,l)<<endl<< "new value:     "<<val<<endl;
      assert( false );
      return;
    }

  twopdm(i, j, k, l) = val;
  twopdm(i, j, l, k) = -val;
  twopdm(j, i, k, l) = -val;
  twopdm(j, i, l, k) = val;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void assign_twopdm_elements( std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements, array_4d<double> & twopdm )
{
  for (int i=0; i < new_spin_orbital_elements.size(); ++i) {
    int ix = new_spin_orbital_elements[i].first[0];
    int jx = new_spin_orbital_elements[i].first[1];
    int kx = new_spin_orbital_elements[i].first[2];
    int lx = new_spin_orbital_elements[i].first[3];
    double x = new_spin_orbital_elements[i].second;
    assign_twopdm_antisymmetric(twopdm, ix, jx, kx, lx, x);
  }

//FIXME is the transpose always needed?
  for (int i=0; i < new_spin_orbital_elements.size(); ++i) {
    int ix = new_spin_orbital_elements[i].first[3];
    int jx = new_spin_orbital_elements[i].first[2];
    int kx = new_spin_orbital_elements[i].first[1];
    int lx = new_spin_orbital_elements[i].first[0];
    double x = new_spin_orbital_elements[i].second;
    assign_twopdm_antisymmetric(twopdm, ix, jx, kx, lx, x);
  }
}

//===========================================================================================================================================================

void twopdm_loop_over_block_operators( Wavefunction & wavefunction, 
                                       const SpinBlock & big, 
                                       std::vector<Npdm::CD> & lhs_cd_type,
                                       std::vector<Npdm::CD> & dot_cd_type,
                                       std::vector<Npdm::CD> & rhs_cd_type,
                                       array_4d<double> & twopdm )
{
  SpinBlock* rhsBlock = big.get_rightBlock();
  SpinBlock* lhsdotBlock = big.get_leftBlock();
  SpinBlock* lhsBlock = lhsdotBlock->get_leftBlock();
  SpinBlock* dotBlock = lhsdotBlock->get_rightBlock();

  boost::shared_ptr<NpdmSpinOps> lhsOps = select_op_wrapper( lhsBlock, lhs_cd_type );
  boost::shared_ptr<NpdmSpinOps> rhsOps = select_op_wrapper( rhsBlock, rhs_cd_type );
  boost::shared_ptr<NpdmSpinOps> dotOps = select_op_wrapper( dotBlock, dot_cd_type );

  Npdm::Npdm_expectations npdm_expectations( wavefunction, big, *lhsOps, *dotOps, *rhsOps );

  // Only one spatial combination on the dot block
  assert( dotOps->size() == 1 );
  bool skip = dotOps->set_local_ops( 0 );
  if (skip) return;
  if ( lhsOps->opReps_.size() > 0 ) assert( dotOps->mults_.size() == dotOps->opReps_.size() );

  // Many spatial combinations on left block
  for ( int ilhs = 0; ilhs < lhsOps->size(); ++ilhs ) {
    skip = lhsOps->set_local_ops( ilhs );
    if (skip) continue;
    if ( lhsOps->opReps_.size() > 0 ) assert( lhsOps->mults_.size() == lhsOps->opReps_.size() );

    // Many spatial combinations on right block
    for ( int irhs = 0; irhs < rhsOps->size(); ++irhs ) {
      skip = rhsOps->set_local_ops( irhs );
      if (skip) continue;
      if ( rhsOps->opReps_.size() > 0 ) assert( rhsOps->mults_.size() == rhsOps->opReps_.size() );

      // Get non-spin-adapated 2PDM elements after building spin-adapted elements
      // FIXME magic number 6
      std::vector< std::pair< std::vector<int>, double > > new_spin_orbital_elements = npdm_expectations.get_nonspin_adapted_expectations( 6 );

      // Assign twopdm elements
      assign_twopdm_elements( new_spin_orbital_elements, twopdm );
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void compute_twopdm_sweep(std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int state, int sweepPos, int endPos)
{

pout << "Sweep position = "<< sweepPos << " (" << endPos << ")\n";

  // 2pdm array built so far
  array_4d<double> twopdm(big.size()*2, big.size()*2, big.size()*2, big.size()*2);
  load_twopdm_binary(twopdm, state, state);
  
  // Loop over NPDM operator patterns (here we initialize for 2PDM)
  Npdm::Npdm_patterns npdm_patterns( 2, sweepPos, endPos );

  for (auto pattern = npdm_patterns.ldr_cd_begin(); pattern != npdm_patterns.ldr_cd_end(); ++pattern) {
    std::vector<Npdm::CD> lhs_cd_type = pattern->at('l');
    std::vector<Npdm::CD> dot_cd_type = pattern->at('d');
    std::vector<Npdm::CD> rhs_cd_type = pattern->at('r');
    // Compute all irreducible PDM elements generated by this block operator pattern at this sweep position
    twopdm_loop_over_block_operators( wavefunctions.at(0), big, lhs_cd_type, dot_cd_type, rhs_cd_type, twopdm );
  }
  
  // Combine NPDM elements from this sweep point with others
  accumulate_twopdm(twopdm);
  save_twopdm_binary(twopdm, state, state);

}

//===========================================================================================================================================================

}

