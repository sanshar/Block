///*                                                                           
//Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
//Copyright (c) 2012, Garnet K.-L. Chan                                        
//                                                                             
//This program is integrated in Molpro with the permission of 
//Sandeep Sharma and Garnet K.-L. Chan
//*/
//
//#include <algorithm>
//#include <boost/format.hpp>
//#ifndef SERIAL
//#include <boost/mpi/environment.hpp>
//#include <boost/mpi/communicator.hpp>
//#include <boost/mpi.hpp>
//#endif
//#include "execinfo.h"
//
//#include "fourpdm.h"
//#include "npdm_patterns.h"
//#include "npdm_expectations.h"
//#include "npdm_operator_wrappers.h"
//#include "npdm_epermute.h"
//
//namespace SpinAdapted{
//
//// Forward declaration
//boost::shared_ptr<NpdmSpinOps> select_op_wrapper( SpinBlock * spinBlock, std::vector<Npdm::CD> & cd_type );
//
////===========================================================================================================================================================
//
//void save_fourpdm_text(const array_8d<double>& fourpdm, const int &i, const int &j)
//{
//  if(!mpigetrank())
//  {
//    char file[5000];
//    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/fourpdm.", i, j);
//    ofstream ofs(file);
//    ofs << fourpdm.dim1() << endl;
//
//    for(int i=0; i<fourpdm.dim1(); ++i)
//      for(int j=0; j<fourpdm.dim2(); ++j)
//        for(int k=0; k<fourpdm.dim3(); ++k)
//          for(int l=0; l<fourpdm.dim4(); ++l)
//            for(int m=0; m<fourpdm.dim5(); ++m)
//              for(int n=0; n<fourpdm.dim6(); ++n)
//                for(int p=0; p<fourpdm.dim7(); ++p)
//                  for(int q=0; q<fourpdm.dim8(); ++q)
//                    if ( abs(fourpdm(i,j,k,l,m,n,p,q)) > 1e-14 )
//                      ofs << boost::format("%d %d %d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % p % q % fourpdm(i,j,k,l,m,n,p,q);
//    ofs.close();
//  }
//}
//
////-----------------------------------------------------------------------------------------------------------------------------------------------------------
//
//void save_spatial_fourpdm_text(const array_8d<double>& fourpdm, const int &i, const int &j)
//{
//std::cout << "Building spatial 4pdm\n";
//  double factor = 1.0;
//  if(!mpigetrank())
//  {
//    char file[5000];
//    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/spatial_fourpdm.", i, j);
//    ofstream ofs(file);
//    ofs << fourpdm.dim1()/2 << endl;
//
//    for(int i=0; i<fourpdm.dim1()/2; ++i)
//      for(int j=0; j<fourpdm.dim2()/2; ++j)
//        for(int k=0; k<fourpdm.dim3()/2; ++k)
//          for(int l=0; l<fourpdm.dim4()/2; ++l)
//            for(int m=0; m<fourpdm.dim5()/2; ++m)
//              for(int n=0; n<fourpdm.dim6()/2; ++n)
//                for(int p=0; p<fourpdm.dim7()/2; ++p)
//                  for(int q=0; q<fourpdm.dim8()/2; ++q) {
//
//                    double pdm = 0.0;
//                    for (int s=0; s<2; s++)
//                      for (int t=0; t<2; t++)
//                        for (int u=0; u<2; u++)
//                          for (int v=0; v<2; v++) {
//                            pdm += factor * fourpdm(2*i+s, 2*j+t, 2*k+u, 2*l+v, 2*m+v, 2*n+u, 2*p+t, 2*q+s);
//                              //if ( (i==1)&& (j==1)&& (k==0)&& (l==0)&& (m==2)&& (n==2) &&(p==1) &&(q==1) ) {
//                              //if ( (i==1)&& (j==1)&& (k==0)&& (l==0)&& (m==1)&& (n==1) &&(p==2) &&(q==2) ) {
//                              //  std::cout <<2*i+s<<","<<2*j+t<<","<<2*k+u<<","<<2*l+v<<","<<2*m+v<<","<<2*n+u<<","<<2*p+t<<","<<2*q+s<<"\t\t";
//                              //  std::cout << pdm << "\t\t" <<  fourpdm(2*i+s, 2*j+t, 2*k+u, 2*l+v, 2*m+v, 2*n+u, 2*p+t, 2*q+s) << std::endl;
//                              //}
//                          }
//        
//                    if ( abs(pdm) > 1e-14 ) ofs << boost::format("%d %d %d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % p % q % pdm;
//                  }
//    ofs.close();
//  }
//}
//
////-----------------------------------------------------------------------------------------------------------------------------------------------------------
//
//void save_spatial_fourpdm_binary(const array_8d<double>& fourpdm, const int &i, const int &j)
//{
//  assert(false);
////  double factor = 1.0;
////  if(!mpigetrank())
////  {
////    char file[5000];
////    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/spatial_binary_fourpdm.", i, j);
////    FILE* f = fopen(file, "wb");
////
////    int nrows = fourpdm.dim1()/2;
////    array_8d<double> pdm; 
////    pdm.resize(nrows, nrows, nrows, nrows, nrows, nrows);
////
////    for(int i=0; i<fourpdm.dim1()/2; ++i)
////      for(int j=0; j<fourpdm.dim2()/2; ++j)
////        for(int k=0; k<fourpdm.dim3()/2; ++k)
////          for(int l=0; l<fourpdm.dim4()/2; ++l)
////            for(int m=0; m<fourpdm.dim5()/2; ++m)
////              for(int n=0; n<fourpdm.dim6()/2; ++n) {
////
////                pdm(i, j, k, l, m, n) = 0.0;
////                for (int s=0; s<2; s++)
////                  for (int t =0; t<2; t++)
////                    for (int u =0; u<2; u++)
////                      pdm(i, j, k, l, m, n) += factor * fourpdm(2*i+s, 2*j+t, 2*k+u, 2*l+u, 2*m+t, 2*n+s) * factor;
////              }
////
////    int result = fwrite(&nrows,  1, sizeof(int), f);
////    result = fwrite(&pdm(0,0,0,0,0,0), pdm.size(), sizeof(double), f);
////    fclose(f);
////  }
//}
//
////-----------------------------------------------------------------------------------------------------------------------------------------------------------
//
//void save_fourpdm_binary(const array_8d<double>& fourpdm, const int &i, const int &j)
//{
//  if(!mpigetrank())
//  {
//    char file[5000];
//    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/fourpdm.", i, j);
//    std::ofstream ofs(file, std::ios::binary);
//    boost::archive::binary_oarchive save(ofs);
//    save << fourpdm;
//    ofs.close();
//  }
//}
//
////-----------------------------------------------------------------------------------------------------------------------------------------------------------
//
//void load_fourpdm_binary(array_8d<double>& fourpdm, const int &i, const int &j)
//{
//  if(!mpigetrank())
//  {
//    char file[5000];
//    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/fourpdm.", i, j);
//    std::ifstream ifs(file, std::ios::binary);
//    boost::archive::binary_iarchive load(ifs);
//    load >> fourpdm;
//    ifs.close();
//  }
//#ifndef SERIAL
//  mpi::communicator world;
//  mpi::broadcast(world,fourpdm,0);
//  if(mpigetrank())
//    fourpdm.Clear();
//#endif
//}
//
////-----------------------------------------------------------------------------------------------------------------------------------------------------------
//
//void save_averaged_fourpdm(const int &nroots)
//{
/////  NYI
//  assert(false);
//}
//
////-----------------------------------------------------------------------------------------------------------------------------------------------------------
//
//void accumulate_fourpdm(array_8d<double>& fourpdm)
//{
//#ifndef SERIAL
//  array_8d<double> tmp_recv;
//  mpi::communicator world;
//  if (!mpigetrank())
//  {
//    for(int u=1; u<world.size(); ++u) {
//      world.recv(u, u, tmp_recv);
//      for(int i=0; i<fourpdm.dim1()/2; ++i)
//        for(int j=0; j<fourpdm.dim2()/2; ++j)
//          for(int k=0; k<fourpdm.dim3()/2; ++k)
//            for(int l=0; l<fourpdm.dim4()/2; ++l)
//              for(int m=0; m<fourpdm.dim5()/2; ++m)
//                for(int n=0; n<fourpdm.dim6()/2; ++n) 
//                  for(int p=0; p<fourpdm.dim7()/2; ++p) 
//                    for(int q=0; q<fourpdm.dim8()/2; ++q) 
//                      if ( tmp_recv(i,j,k,l,m,n,p,q) != 0.0 ) fourpdm(i,j,k,l,m,n,p,q) = tmp_recv(i,j,k,l,m,n,p,q);
//    }
//  }
//  else 
//  {
//    world.send(0, mpigetrank(), fourpdm);
//  }
//#endif
//}
//
////-----------------------------------------------------------------------------------------------------------------------------------------------------------
//
//void get_even_and_odd_perms( const std::vector<int> mnpq, std::vector< std::vector<int> > & even_perms, std::vector< std::vector<int> > & odd_perms )
//{
//
//  // Get all even and odd mnpq permutations
//  bool even = false;
//
//  // Must sort them to get all possible permutations
//  std::vector<int> foo = mnpq;
//  std::sort( foo.begin(), foo.end() );
//
//  // Get first set
//  std::vector< std::vector<int> > perms1;
//  do { 
//    perms1.push_back( foo ); 
//    if (foo == mnpq) even = true; 
//  } while ( next_even_permutation(foo.begin(), foo.end()) );
//
//  // Re-sort and swap LAST TWO elements to ensure we get all the remaining permutations
//  std::sort( foo.begin(), foo.end() );
//  assert( foo.size() == 4 );
//  std::swap( foo[2], foo[3] );
//
//  // Get second set
//  std::vector< std::vector<int> > perms2;
//  do { 
//    perms2.push_back( foo ); 
//  } while ( next_even_permutation(foo.begin(), foo.end()) );
//
//  // Assign as even or odd permutations
//  even_perms = perms1;
//  odd_perms = perms2;
//  if (!even) std::swap( even_perms, odd_perms );
//
//}
//
////-----------------------------------------------------------------------------------------------------------------------------------------------------------
//// The number of possible combinations are (4!)**2
//
//void assign_fourpdm_antisymmetric(array_8d<double>& fourpdm, 
//                                  const int i, const int j, const int k, const int l, const int m, const int n, const int p, const int q,
//                                  const double val)
//{
////if ( abs(val) > 1e-8 ) {
////  pout << "so-fourpdm val: i,j,k,l,m,n,p,q = " 
////       << i << "," << j << "," << k << "," << l << "," << m << "," << n << "," << p << "," << q
////       << "\t\t" << val << endl;
////}
//
//  // Test for duplicates
//  if ( fourpdm(i,j,k,l,m,n,p,q) != 0.0 && abs(fourpdm(i,j,k,l,m,n,p,q)-val) > 1e-6) {
//    void *array[10];
//    size_t size;
//    size = backtrace(array, 10);
//    cout << "WARNING: Already calculated "<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<m<<" "<<n<<" "<<p<<" "<<q<<endl;
//    //backtrace_symbols_fd(array, size, 2);
//    cout << "earlier value: " << fourpdm(i,j,k,l,m,n,p,q) << endl << "new value:     " <<val<<endl;
//    assert( false );
//    return;
//  }
//
//  if ( abs(val) < 1e-14 ) return;
//
//  // If indices are not all unique, then all elements should be zero (and next_even_permutation fails)
//  std::vector<int> v = {i,j,k,l};
//  std::sort( v.begin(), v.end() );
//  if ( (v[0]==v[1]) || (v[1]==v[2]) || (v[2]==v[3]) ) { if (abs(val) > 1e-15) assert(false); return; }
//  std::vector<int> w = {m,n,p,q};
//  std::sort( w.begin(), w.end() );
//  if ( (w[0]==w[1]) || (w[1]==w[2]) || (w[2]==w[3]) ) { if (abs(val) > 1e-15) assert(false); return; }
//
//  // Get all even and odd permutations
//  const std::vector<int> ijkl = {i,j,k,l};
//  std::vector< std::vector<int> > ijkl_even, ijkl_odd;
//  get_even_and_odd_perms( ijkl, ijkl_even, ijkl_odd );
//  assert ( ijkl_even.size() + ijkl_odd.size() == 24 );
//
//  const std::vector<int> mnpq = {m,n,p,q};
//  std::vector< std::vector<int> > mnpq_even, mnpq_odd;
//  get_even_and_odd_perms( mnpq, mnpq_even, mnpq_odd );
//  assert ( mnpq_even.size() + mnpq_odd.size() == 24 );
//
//  // Even-Even terms
//  for ( auto u = ijkl_even.begin(); u != ijkl_even.end(); ++u )
//    for ( auto v = mnpq_even.begin(); v != mnpq_even.end(); ++v )
//      fourpdm( (*u)[0],(*u)[1],(*u)[2],(*u)[3], (*v)[0],(*v)[1],(*v)[2],(*v)[3] ) = val;
//  // Even-Odd terms
//  for ( auto u = ijkl_even.begin(); u != ijkl_even.end(); ++u )
//    for ( auto v = mnpq_odd.begin(); v != mnpq_odd.end(); ++v )
//      fourpdm( (*u)[0],(*u)[1],(*u)[2],(*u)[3], (*v)[0],(*v)[1],(*v)[2],(*v)[3] ) = -val;
//  // Odd-Even terms
//  for ( auto u = ijkl_odd.begin(); u != ijkl_odd.end(); ++u )
//    for ( auto v = mnpq_even.begin(); v != mnpq_even.end(); ++v )
//      fourpdm( (*u)[0],(*u)[1],(*u)[2],(*u)[3], (*v)[0],(*v)[1],(*v)[2],(*v)[3] ) = -val;
//  // Odd-Odd terms
//  for ( auto u = ijkl_odd.begin(); u != ijkl_odd.end(); ++u )
//    for ( auto v = mnpq_odd.begin(); v != mnpq_odd.end(); ++v )
//      fourpdm( (*u)[0],(*u)[1],(*u)[2],(*u)[3], (*v)[0],(*v)[1],(*v)[2],(*v)[3] ) = val;
//
//}
//
////-----------------------------------------------------------------------------------------------------------------------------------------------------------
//
//void assign_fourpdm_elements( std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements, array_8d<double> & fourpdm )
//{
//  for (int i=0; i < new_spin_orbital_elements.size(); ++i) {
//    assert( new_spin_orbital_elements[i].first.size() == 8 );
//    int ix = new_spin_orbital_elements[i].first[0];
//    int jx = new_spin_orbital_elements[i].first[1];
//    int kx = new_spin_orbital_elements[i].first[2];
//    int lx = new_spin_orbital_elements[i].first[3];
//    int mx = new_spin_orbital_elements[i].first[4];
//    int nx = new_spin_orbital_elements[i].first[5];
//    int px = new_spin_orbital_elements[i].first[6];
//    int qx = new_spin_orbital_elements[i].first[7];
//    double x = new_spin_orbital_elements[i].second;
//    assign_fourpdm_antisymmetric(fourpdm, ix, jx, kx, lx, mx, nx, px, qx, x);
//  }
//
//  for (int i=0; i < new_spin_orbital_elements.size(); ++i) {
//    assert( new_spin_orbital_elements[i].first.size() == 8 );
//    int ix = new_spin_orbital_elements[i].first[7];
//    int jx = new_spin_orbital_elements[i].first[6];
//    int kx = new_spin_orbital_elements[i].first[5];
//    int lx = new_spin_orbital_elements[i].first[4];
//    int mx = new_spin_orbital_elements[i].first[3];
//    int nx = new_spin_orbital_elements[i].first[2];
//    int px = new_spin_orbital_elements[i].first[1];
//    int qx = new_spin_orbital_elements[i].first[0];
//    double x = new_spin_orbital_elements[i].second;
//    assign_fourpdm_antisymmetric(fourpdm, ix, jx, kx, lx, mx, nx, px, qx, x);
//  }
//}
//
////===========================================================================================================================================================
//
//void fourpdm_loop_over_block_operators( Wavefunction & wavefunction, 
//                                         const SpinBlock & big, 
//                                         std::vector<Npdm::CD> & lhs_cd_type,
//                                         std::vector<Npdm::CD> & dot_cd_type,
//                                         std::vector<Npdm::CD> & rhs_cd_type,
//                                         array_8d<double> & fourpdm )
//{
//  SpinBlock* rhsBlock = big.get_rightBlock();
//  SpinBlock* lhsdotBlock = big.get_leftBlock();
//  SpinBlock* lhsBlock = lhsdotBlock->get_leftBlock();
//  SpinBlock* dotBlock = lhsdotBlock->get_rightBlock();
//
//  boost::shared_ptr<NpdmSpinOps> lhsOps = select_op_wrapper( lhsBlock, lhs_cd_type );
//  boost::shared_ptr<NpdmSpinOps> rhsOps = select_op_wrapper( rhsBlock, rhs_cd_type );
//  boost::shared_ptr<NpdmSpinOps> dotOps = select_op_wrapper( dotBlock, dot_cd_type );
//
////  Npdm::Npdm_expectations npdm_expectations( 4, wavefunction, big, *lhsOps, *dotOps, *rhsOps );
//
//  // Only one spatial combination on the dot block
//  assert( dotOps->size() == 1 );
//  bool skip = dotOps->set_local_ops( 0 );
//  if (skip) return;
//  if ( lhsOps->opReps_.size() > 0 ) assert( dotOps->mults_.size() == dotOps->opReps_.size() );
//
//  // Many spatial combinations on left block
//  for ( int ilhs = 0; ilhs < lhsOps->size(); ++ilhs ) {
//    skip = lhsOps->set_local_ops( ilhs );
//    if (skip) continue;
//    if ( lhsOps->opReps_.size() > 0 ) assert( lhsOps->mults_.size() == lhsOps->opReps_.size() );
//
//    // Many spatial combinations on right block
//    for ( int irhs = 0; irhs < rhsOps->size(); ++irhs ) {
//      skip = rhsOps->set_local_ops( irhs );
//      if (skip) continue;
//      if ( rhsOps->opReps_.size() > 0 ) assert( rhsOps->mults_.size() == rhsOps->opReps_.size() );
//
//      // Get non-spin-adapated 4PDM elements after building spin-adapted elements
////      std::vector< std::pair< std::vector<int>, double > > new_spin_orbital_elements = npdm_expectations.get_nonspin_adapted_expectations();
//
//      // Assign fourpdm elements
////      assign_fourpdm_elements( new_spin_orbital_elements, fourpdm );
//    }
//  }
//}
//
////-----------------------------------------------------------------------------------------------------------------------------------------------------------
//
//void compute_fourpdm_sweep(std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int state, int sweepPos, int endPos)
//{
//
//  pout << "===========================================================================================\n";
//  pout << "4PDM sweep position = "<< sweepPos << " (" << endPos << ")\n";
//
//  // 4pdm array built so far
//  int dim = 2*big.size();
//  array_8d<double> fourpdm(dim,dim,dim,dim,dim,dim,dim,dim);
//  load_fourpdm_binary(fourpdm, state, state);
//  
//  // Loop over NPDM operator patterns (here we initialize for 4PDM)
//  Npdm::Npdm_patterns npdm_patterns( 4, sweepPos, endPos );
//
//  for (auto pattern = npdm_patterns.ldr_cd_begin(); pattern != npdm_patterns.ldr_cd_end(); ++pattern) {
//
//    pout << "===========================================================================================\n";
//    std::cout << "Doing pattern:\n";
//    npdm_patterns.print_cd_string( pattern->at('l') );
//if (pattern->at('l').size() == 4) continue;
//    npdm_patterns.print_cd_string( pattern->at('d') );
//    npdm_patterns.print_cd_string( pattern->at('r') );
//    std::cout << std::endl;
//
//    std::vector<Npdm::CD> lhs_cd_type = pattern->at('l');
//    std::vector<Npdm::CD> dot_cd_type = pattern->at('d');
//    std::vector<Npdm::CD> rhs_cd_type = pattern->at('r');
//    // Compute all irreducible PDM elements generated by this block operator pattern at this sweep position
//if ( lhs_cd_type.size() == 2  &&
//     dot_cd_type.size() == 4  &&
//     rhs_cd_type.size() == 2 )
//    fourpdm_loop_over_block_operators( wavefunctions.at(0), big, lhs_cd_type, dot_cd_type, rhs_cd_type, fourpdm );
//  }
//  
//  // Combine NPDM elements from this sweep point with others
//  accumulate_fourpdm(fourpdm);
//  save_fourpdm_binary(fourpdm, state, state);
//  pout << "End of 4PDM at this sweep position\n";
//  pout << "===========================================================================================\n";
//
//}
//
////===========================================================================================================================================================
//
//}
//
