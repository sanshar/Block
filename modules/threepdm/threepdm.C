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

#include "threepdm.h"
#include "npdm_patterns.h"
#include "npdm_expectations.h"
#include "npdm_operator_wrappers.h"
#include "npdm_epermute.h"

namespace SpinAdapted{

// Forward declaration
boost::shared_ptr<NpdmSpinOps> select_op_wrapper( SpinBlock * spinBlock, std::vector<Npdm::CD> & cd_type );

//===========================================================================================================================================================

void save_threepdm_text(const array_6d<double>& threepdm, const int &i, const int &j)
{
  if(!mpigetrank())
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/threepdm.", i, j);
    ofstream ofs(file);
    ofs << threepdm.dim1() << endl;

    for(int i=0; i<threepdm.dim1(); ++i)
      for(int j=0; j<threepdm.dim2(); ++j)
        for(int k=0; k<threepdm.dim3(); ++k)
          for(int l=0; l<threepdm.dim4(); ++l)
            for(int m=0; m<threepdm.dim5(); ++m)
              for(int n=0; n<threepdm.dim6(); ++n)
                if ( abs(threepdm(i,j,k,l,m,n)) > 1e-14 )
                  ofs << boost::format("%d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % threepdm(i,j,k,l,m,n);
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void save_spatial_threepdm_text(const array_6d<double>& threepdm, const int &i, const int &j)
{
//FIXME  2pdm spatial has a factor of 1/2 in front of it  ???
  double factor = 1.0;
  if(!mpigetrank())
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/spatial_threepdm.", i, j);
    ofstream ofs(file);
    ofs << threepdm.dim1()/2 << endl;

    for(int i=0; i<threepdm.dim1()/2; ++i)
      for(int j=0; j<threepdm.dim2()/2; ++j)
        for(int k=0; k<threepdm.dim3()/2; ++k)
          for(int l=0; l<threepdm.dim4()/2; ++l)
            for(int m=0; m<threepdm.dim5()/2; ++m)
              for(int n=0; n<threepdm.dim6()/2; ++n) {

                double pdm = 0.0;
                for (int s=0; s<2; s++)
                  for (int t=0; t<2; t++)
                    for (int u=0; u<2; u++)
                      pdm += factor * threepdm(2*i+s, 2*j+t, 2*k+u, 2*l+u, 2*m+t, 2*n+s);
    
                if ( abs(pdm) > 1e-14 ) ofs << boost::format("%d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % pdm;
              }
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void save_spatial_threepdm_binary(const array_6d<double>& threepdm, const int &i, const int &j)
{
//FIXME  2pdm spatial has a factor of 1/2 in front of it  ???
  double factor = 1.0;
  if(!mpigetrank())
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/spatial_binary_threepdm.", i, j);
    FILE* f = fopen(file, "wb");

    int nrows = threepdm.dim1()/2;
    array_6d<double> pdm; 
    pdm.resize(nrows, nrows, nrows, nrows, nrows, nrows);

    for(int i=0; i<threepdm.dim1()/2; ++i)
      for(int j=0; j<threepdm.dim2()/2; ++j)
        for(int k=0; k<threepdm.dim3()/2; ++k)
          for(int l=0; l<threepdm.dim4()/2; ++l)
            for(int m=0; m<threepdm.dim5()/2; ++m)
              for(int n=0; n<threepdm.dim6()/2; ++n) {

                pdm(i, j, k, l, m, n) = 0.0;
                for (int s=0; s<2; s++)
                  for (int t =0; t<2; t++)
                    for (int u =0; u<2; u++)
                      pdm(i, j, k, l, m, n) += factor * threepdm(2*i+s, 2*j+t, 2*k+u, 2*l+u, 2*m+t, 2*n+s);
              }

    int result = fwrite(&nrows,  1, sizeof(int), f);
    result = fwrite(&pdm(0,0,0,0,0,0), pdm.size(), sizeof(double), f);
    fclose(f);
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void save_threepdm_binary(const array_6d<double>& threepdm, const int &i, const int &j)
{
  if(!mpigetrank())
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/threepdm.", i, j);
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << threepdm;
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void load_threepdm_binary(array_6d<double>& threepdm, const int &i, const int &j)
{
  if(!mpigetrank())
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/threepdm.", i, j);
    std::ifstream ifs(file, std::ios::binary);
    boost::archive::binary_iarchive load(ifs);
    load >> threepdm;
    ifs.close();
  }
#ifndef SERIAL
  mpi::communicator world;
  mpi::broadcast(world,threepdm,0);
  if(mpigetrank())
    threepdm.Clear();
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void save_averaged_threepdm(const int &nroots)
{
///  NYI
  assert(false);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void accumulate_threepdm(array_6d<double>& threepdm)
{
#ifndef SERIAL
  array_6d<double> tmp_recv;
  mpi::communicator world;
  if (!mpigetrank())
  {
    for(int p=1; p<world.size(); ++p) {
      world.recv(p, p, tmp_recv);
      for(int i=0; i<threepdm.dim1()/2; ++i)
        for(int j=0; j<threepdm.dim2()/2; ++j)
          for(int k=0; k<threepdm.dim3()/2; ++k)
            for(int l=0; l<threepdm.dim4()/2; ++l)
              for(int m=0; m<threepdm.dim5()/2; ++m)
                for(int n=0; n<threepdm.dim6()/2; ++n) 
                  if(tmp_recv(i,j,k,l,m,n) != 0.0) threepdm(i,j,k,l,m,n) = tmp_recv(i,j,k,l,m,n);
    }
  }
  else 
  {
    world.send(0, mpigetrank(), threepdm);
  }
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// The number of possible combinations are (3!)**2

void assign_threepdm_antisymmetric(array_6d<double>& threepdm, 
                                   const int i, const int j, const int k, const int l, const int m, const int n, 
                                   const double val)
{
//if ( abs(val) > 1e-8 ) {
//  pout << "so-threepdm val: i,j,k,l,m,n = " 
//       << i << "," << j << "," << k << "," << l << "," << m << "," << n
//       << "\t\t" << val << endl;
//}

  // Test for duplicates
  if ( threepdm(i,j,k,l,m,n) != 0.0 && abs(threepdm(i,j,k,l,m,n)-val) > 1e-6) {
    void *array[10];
    size_t size;
    size = backtrace(array, 10);
    cout << "WARNING: Already calculated "<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<m<<" "<<n<<endl;
    //backtrace_symbols_fd(array, size, 2);
    cout << "earlier value: " << threepdm(i,j,k,l,m,n) << endl << "new value:     " <<val<<endl;
    assert( false );
    return;
  }

  if ( abs(val) < 1e-14 ) return;

  // If indices are not all unique, then all elements should be zero (and next_even_permutation fails)
  std::vector<int> v = {i,j,k};
  std::sort( v.begin(), v.end() );
  if ( (v[0]==v[1]) || (v[1]==v[2]) ) { if (abs(val) > 1e-15) assert(false); return; }
  std::vector<int> w = {l,m,n};
  std::sort( w.begin(), w.end() );
  if ( (w[0]==w[1]) || (w[1]==w[2]) ) { if (abs(val) > 1e-15) assert(false); return; }


  // The number of possible combinations are (3!)**2
  threepdm(i, j, k, l, m, n) =  val;
  threepdm(i, j, k, l, n, m) = -val;
  threepdm(i, j, k, n, m, l) = -val;
  threepdm(i, j, k, n, l, m) =  val;
  threepdm(i, j, k, m, l, n) = -val;
  threepdm(i, j, k, m, n, l) =  val;

  threepdm(i, k, j, l, m, n) = -val;
  threepdm(i, k, j, l, n, m) =  val;
  threepdm(i, k, j, n, m, l) =  val;
  threepdm(i, k, j, n, l, m) = -val;
  threepdm(i, k, j, m, l, n) =  val;
  threepdm(i, k, j, m, n, l) = -val;

  threepdm(j, i, k, l, m, n) = -val;
  threepdm(j, i, k, l, n, m) =  val;
  threepdm(j, i, k, n, m, l) =  val;
  threepdm(j, i, k, n, l, m) = -val;
  threepdm(j, i, k, m, l, n) =  val;
  threepdm(j, i, k, m, n, l) = -val;

  threepdm(j, k, i, l, m, n) =  val;
  threepdm(j, k, i, l, n, m) = -val;
  threepdm(j, k, i, n, m, l) = -val;
  threepdm(j, k, i, n, l, m) =  val;
  threepdm(j, k, i, m, l, n) = -val;
  threepdm(j, k, i, m, n, l) =  val;

  threepdm(k, j, i, l, m, n) = -val;
  threepdm(k, j, i, l, n, m) =  val;
  threepdm(k, j, i, n, m, l) =  val;
  threepdm(k, j, i, n, l, m) = -val;
  threepdm(k, j, i, m, l, n) =  val;
  threepdm(k, j, i, m, n, l) = -val;

  threepdm(k, i, j, l, m, n) =  val;
  threepdm(k, i, j, l, n, m) = -val;
  threepdm(k, i, j, n, m, l) = -val;
  threepdm(k, i, j, n, l, m) =  val;
  threepdm(k, i, j, m, l, n) = -val;
  threepdm(k, i, j, m, n, l) =  val;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void assign_threepdm_elements( std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements, array_6d<double> & threepdm )
{
  for (int i=0; i < new_spin_orbital_elements.size(); ++i) {
    assert( new_spin_orbital_elements[i].first.size() == 6 );
    int ix = new_spin_orbital_elements[i].first[0];
    int jx = new_spin_orbital_elements[i].first[1];
    int kx = new_spin_orbital_elements[i].first[2];
    int lx = new_spin_orbital_elements[i].first[3];
    int mx = new_spin_orbital_elements[i].first[4];
    int nx = new_spin_orbital_elements[i].first[5];
    double x = new_spin_orbital_elements[i].second;
    assign_threepdm_antisymmetric(threepdm, ix, jx, kx, lx, mx, nx, x);
  }

  // Assign transposed elements
  for (int i=0; i < new_spin_orbital_elements.size(); ++i) {
    assert( new_spin_orbital_elements[i].first.size() == 6 );
    int ix = new_spin_orbital_elements[i].first[5];
    int jx = new_spin_orbital_elements[i].first[4];
    int kx = new_spin_orbital_elements[i].first[3];
    int lx = new_spin_orbital_elements[i].first[2];
    int mx = new_spin_orbital_elements[i].first[1];
    int nx = new_spin_orbital_elements[i].first[0];
    double x = new_spin_orbital_elements[i].second;
    assign_threepdm_antisymmetric(threepdm, ix, jx, kx, lx, mx, nx, x);
  }
}

//===========================================================================================================================================================

void threepdm_loop_over_block_operators( Wavefunction & wavefunction, 
                                         const SpinBlock & big, 
                                         std::vector<Npdm::CD> & lhs_cd_type,
                                         std::vector<Npdm::CD> & dot_cd_type,
                                         std::vector<Npdm::CD> & rhs_cd_type,
                                         array_6d<double> & threepdm )
{
  SpinBlock* rhsBlock = big.get_rightBlock();
  SpinBlock* lhsdotBlock = big.get_leftBlock();
  SpinBlock* lhsBlock = lhsdotBlock->get_leftBlock();
  SpinBlock* dotBlock = lhsdotBlock->get_rightBlock();

  boost::shared_ptr<NpdmSpinOps> lhsOps = select_op_wrapper( lhsBlock, lhs_cd_type );
  boost::shared_ptr<NpdmSpinOps> rhsOps = select_op_wrapper( rhsBlock, rhs_cd_type );
  boost::shared_ptr<NpdmSpinOps> dotOps = select_op_wrapper( dotBlock, dot_cd_type );

  Npdm::Npdm_expectations npdm_expectations( wavefunction, big, *lhsOps, *dotOps, *rhsOps );
//pout << "-------------------------------------------------------------------------------------------\n";
//pout << "lhsOps->size()" << lhsOps->size() << std::endl;
//pout << "dotOps->size()" << dotOps->size() << std::endl;
//pout << "rhsOps->size()" << rhsOps->size() << std::endl;
//pout << "------\n";

  // Only one spatial combination on the dot block
  assert( dotOps->size() == 1 );
  bool skip = dotOps->set_local_ops( 0 );
  if (skip) return;
  if ( lhsOps->opReps_.size() > 0 ) assert( dotOps->mults_.size() == dotOps->opReps_.size() );

  // Many spatial combinations on left block
  for ( int ilhs = 0; ilhs < lhsOps->size(); ++ilhs ) {
    skip = lhsOps->set_local_ops( ilhs );
    if (skip) continue;
//pout << "lhsOps->indices_.size() " << lhsOps->indices_.size() << std::endl;
//pout << "lhsOps->opReps_.size() " << lhsOps->opReps_.size() << std::endl;
//pout << "dotOps->indices_.size() " << dotOps->indices_.size() << std::endl;
//pout << "dotOps->opReps_.size() " << dotOps->opReps_.size() << std::endl;
    if ( lhsOps->opReps_.size() > 0 ) assert( lhsOps->mults_.size() == lhsOps->opReps_.size() );

    // Many spatial combinations on right block
    for ( int irhs = 0; irhs < rhsOps->size(); ++irhs ) {
      skip = rhsOps->set_local_ops( irhs );
      if (skip) continue;
//pout << "rhsOps->indices_.size() " << rhsOps->indices_.size() << std::endl;
//pout << "rhsOps->opReps_.size() " << rhsOps->opReps_.size() << std::endl;
//pout << "-------------------------------------------------------------------------------------------\n";
//pout << "spatial: ilhs, irhs = " << ilhs << ", " << irhs << std::endl;
      if ( rhsOps->opReps_.size() > 0 ) assert( rhsOps->mults_.size() == rhsOps->opReps_.size() );

      // Get non-spin-adapated 3PDM elements after building spin-adapted elements
      std::vector< std::pair< std::vector<int>, double > > new_spin_orbital_elements = npdm_expectations.get_nonspin_adapted_expectations( 20 );

      // Assign threepdm elements
      assign_threepdm_elements( new_spin_orbital_elements, threepdm );
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void compute_threepdm_sweep(std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int state, int sweepPos, int endPos)
{

  pout << "===========================================================================================\n";
  pout << "3PDM sweep position = "<< sweepPos << " (" << endPos << ")\n";

  // 3pdm array built so far
  int dim = 2*big.size();
  array_6d<double> threepdm(dim,dim,dim,dim,dim,dim);
  load_threepdm_binary(threepdm, state, state);
  
  // Loop over NPDM operator patterns (here we initialize for 3PDM)
  Npdm::Npdm_patterns npdm_patterns( 3, sweepPos, endPos );

  for (auto pattern = npdm_patterns.ldr_cd_begin(); pattern != npdm_patterns.ldr_cd_end(); ++pattern) {

    pout << "===========================================================================================\n";
    std::cout << "Doing pattern:\n";
    npdm_patterns.print_cd_string( pattern->at('l') );
    npdm_patterns.print_cd_string( pattern->at('d') );
    npdm_patterns.print_cd_string( pattern->at('r') );
    std::cout << std::endl;

    std::vector<Npdm::CD> lhs_cd_type = pattern->at('l');
    std::vector<Npdm::CD> dot_cd_type = pattern->at('d');
    std::vector<Npdm::CD> rhs_cd_type = pattern->at('r');
//FIXME
//std::vector<Npdm::CD> foo = {Npdm::DESTRUCTION,Npdm::CREATION,Npdm::DESTRUCTION};
//if ( (lhs_cd_type == foo) || (rhs_cd_type == foo) || (dot_cd_type == foo) ) continue;
//if ( lhs_cd_type.size() == 2  &&
//     dot_cd_type.size() == 2  &&
//     rhs_cd_type.size() == 2 )
    // Compute all irreducible PDM elements generated by this block operator pattern at this sweep position
    threepdm_loop_over_block_operators( wavefunctions.at(0), big, lhs_cd_type, dot_cd_type, rhs_cd_type, threepdm );
  }
  
  // Combine NPDM elements from this sweep point with others
  accumulate_threepdm(threepdm);
  save_threepdm_binary(threepdm, state, state);

}

//===========================================================================================================================================================

}

