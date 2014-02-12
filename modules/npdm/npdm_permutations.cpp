/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "npdm_permutations.h"

namespace SpinAdapted{

//===========================================================================================================================================================

std::map< std::tuple<int,int,int,int>, int > Twopdm_permutations::get_spin_permutations( const std::vector<int>& indices )
{

  std::map< std::tuple<int,int,int,int>, int > perms;
  int i = indices[0];
  int j = indices[1];
  int k = indices[2];
  int l = indices[3];

  // 8 permutations
  //--------------------------
  perms[std::make_tuple(i, j, k, l)] = 1;
  perms[std::make_tuple(i, j, l, k)] = -1;
  perms[std::make_tuple(j, i, k, l)] = -1;
  perms[std::make_tuple(j, i, l, k)] = 1;

  perms[std::make_tuple(l, k, j, i)] = 1;
  perms[std::make_tuple(k, l, j, i)] = -1;
  perms[std::make_tuple(l, k, i, j)] = -1;
  perms[std::make_tuple(k, l, i, j)] = 1;

  return perms;
}

//===========================================================================================================================================================

std::map< std::vector<int>, int > Threepdm_permutations::get_spin_permutations( const std::vector<int>& indices )
{
  assert( indices.size() == 6 );
  std::map< std::vector<int>, int > perms;
  std::vector<int> idx;
  int i = indices[0];
  int j = indices[1];
  int k = indices[2];
  int l = indices[3];
  int m = indices[4];
  int n = indices[5];

  // The number of possible combinations is (3!)**2
  //------------------------------------------------

  idx = { i, j, k, l, m, n }; perms[ idx ] =  1;
  idx = { i, j, k, l, n, m }; perms[ idx ] = -1;
  idx = { i, j, k, n, m, l }; perms[ idx ] = -1;
  idx = { i, j, k, n, l, m }; perms[ idx ] =  1;
  idx = { i, j, k, m, l, n }; perms[ idx ] = -1;
  idx = { i, j, k, m, n, l }; perms[ idx ] =  1;

  idx = { i, k, j, l, m, n }; perms[ idx ] = -1;
  idx = { i, k, j, l, n, m }; perms[ idx ] =  1;
  idx = { i, k, j, n, m, l }; perms[ idx ] =  1;
  idx = { i, k, j, n, l, m }; perms[ idx ] = -1;
  idx = { i, k, j, m, l, n }; perms[ idx ] =  1;
  idx = { i, k, j, m, n, l }; perms[ idx ] = -1;

  idx = { j, i, k, l, m, n }; perms[ idx ] = -1;
  idx = { j, i, k, l, n, m }; perms[ idx ] =  1;
  idx = { j, i, k, n, m, l }; perms[ idx ] =  1;
  idx = { j, i, k, n, l, m }; perms[ idx ] = -1;
  idx = { j, i, k, m, l, n }; perms[ idx ] =  1;
  idx = { j, i, k, m, n, l }; perms[ idx ] = -1;

  idx = { j, k, i, l, m, n }; perms[ idx ] =  1;
  idx = { j, k, i, l, n, m }; perms[ idx ] = -1;
  idx = { j, k, i, n, m, l }; perms[ idx ] = -1;
  idx = { j, k, i, n, l, m }; perms[ idx ] =  1;
  idx = { j, k, i, m, l, n }; perms[ idx ] = -1;
  idx = { j, k, i, m, n, l }; perms[ idx ] =  1;

  idx = { k, j, i, l, m, n }; perms[ idx ] = -1;
  idx = { k, j, i, l, n, m }; perms[ idx ] =  1;
  idx = { k, j, i, n, m, l }; perms[ idx ] =  1;
  idx = { k, j, i, n, l, m }; perms[ idx ] = -1;
  idx = { k, j, i, m, l, n }; perms[ idx ] =  1;
  idx = { k, j, i, m, n, l }; perms[ idx ] = -1;

  idx = { k, i, j, l, m, n }; perms[ idx ] =  1;
  idx = { k, i, j, l, n, m }; perms[ idx ] = -1;
  idx = { k, i, j, n, m, l }; perms[ idx ] = -1;
  idx = { k, i, j, n, l, m }; perms[ idx ] =  1;
  idx = { k, i, j, m, l, n }; perms[ idx ] = -1;
  idx = { k, i, j, m, n, l }; perms[ idx ] =  1;

  // Get transpose elements with same parity factors
  std::map< std::vector<int>, int > trans_perms;
  for (auto it = perms.begin(); it != perms.end(); ++it) {
    std::vector<int> indices = it->first;
    std::reverse( indices.begin(), indices.end() );
    trans_perms[ indices ] = it->second;
  }

  // Now bundle them togther
  for (auto it = trans_perms.begin(); it != trans_perms.end(); ++it) {
    perms[ it->first ] = it->second;
  }

  return perms;
}

//===========================================================================================================================================================

void Fourpdm_permutations::get_even_and_odd_perms( const std::vector<int> mnpq, 
                                                   std::vector< std::vector<int> > & even_perms, 
                                                   std::vector< std::vector<int> > & odd_perms )
{
  // Get all even and odd mnpq permutations
  bool even = false;

  // Must sort them to get all possible permutations
  std::vector<int> foo = mnpq;
  std::sort( foo.begin(), foo.end() );

  // Get first set
  std::vector< std::vector<int> > perms1;
  do { 
    perms1.push_back( foo ); 
    if (foo == mnpq) even = true; 
  } while ( next_even_permutation(foo.begin(), foo.end()) );

  // Re-sort and swap LAST TWO elements to ensure we get all the remaining permutations
  std::sort( foo.begin(), foo.end() );
  assert( foo.size() == 4 );
  std::swap( foo[2], foo[3] );

  // Get second set
  std::vector< std::vector<int> > perms2;
  do { 
    perms2.push_back( foo ); 
  } while ( next_even_permutation(foo.begin(), foo.end()) );

  // Assign as even or odd permutations
  even_perms = perms1;
  odd_perms = perms2;
  if (!even) std::swap( even_perms, odd_perms );

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

std::map< std::vector<int>, int > Fourpdm_permutations::get_spin_permutations( const std::vector<int>& indices )
{
  assert( indices.size() == 8 );
  std::map< std::vector<int>, int > perms;
  std::vector<int> idx;
  int i = indices[0];
  int j = indices[1];
  int k = indices[2];
  int l = indices[3];
  int m = indices[4];
  int n = indices[5];
  int p = indices[6];
  int q = indices[7];

  // If indices are not all unique, then all elements should be zero (and next_even_permutation fails)
  std::vector<int> v = {i,j,k,l};
  std::sort( v.begin(), v.end() );
  if ( (v[0]==v[1]) || (v[1]==v[2]) || (v[2]==v[3]) ) return perms;
  std::vector<int> w = {m,n,p,q};
  std::sort( w.begin(), w.end() );
  if ( (w[0]==w[1]) || (w[1]==w[2]) || (w[2]==w[3]) ) return perms;

  // The number of possible combinations is (4!)**2 
  //------------------------------------------------

  // Get all even and odd permutations
  const std::vector<int> ijkl = {i,j,k,l};
  std::vector< std::vector<int> > ijkl_even, ijkl_odd;
  get_even_and_odd_perms( ijkl, ijkl_even, ijkl_odd );
  assert ( ijkl_even.size() + ijkl_odd.size() == 24 );

  const std::vector<int> mnpq = {m,n,p,q};
  std::vector< std::vector<int> > mnpq_even, mnpq_odd;
  get_even_and_odd_perms( mnpq, mnpq_even, mnpq_odd );
  assert ( mnpq_even.size() + mnpq_odd.size() == 24 );

  // Even-Even terms
  for ( auto u = ijkl_even.begin(); u != ijkl_even.end(); ++u ) {
    for ( auto v = mnpq_even.begin(); v != mnpq_even.end(); ++v ) {
      idx = { (*u)[0],(*u)[1],(*u)[2],(*u)[3], (*v)[0],(*v)[1],(*v)[2],(*v)[3] }; perms[ idx ] = 1;
    }
  }
  // Even-Odd terms
  for ( auto u = ijkl_even.begin(); u != ijkl_even.end(); ++u ) {
    for ( auto v = mnpq_odd.begin(); v != mnpq_odd.end(); ++v ) {
      idx = { (*u)[0],(*u)[1],(*u)[2],(*u)[3], (*v)[0],(*v)[1],(*v)[2],(*v)[3] }; perms[ idx ] = -1;
    }
  }
  // Odd-Even terms
  for ( auto u = ijkl_odd.begin(); u != ijkl_odd.end(); ++u ) {
    for ( auto v = mnpq_even.begin(); v != mnpq_even.end(); ++v ) {
      idx = { (*u)[0],(*u)[1],(*u)[2],(*u)[3], (*v)[0],(*v)[1],(*v)[2],(*v)[3] }; perms[ idx ] = -1;
    }
  }
  // Odd-Odd terms
  for ( auto u = ijkl_odd.begin(); u != ijkl_odd.end(); ++u ) {
    for ( auto v = mnpq_odd.begin(); v != mnpq_odd.end(); ++v ) {
      idx = { (*u)[0],(*u)[1],(*u)[2],(*u)[3], (*v)[0],(*v)[1],(*v)[2],(*v)[3] }; perms[ idx ] = 1;
    }
  }

  // Get transpose elements with same parity factors
  std::map< std::vector<int>, int > trans_perms;
  for (auto it = perms.begin(); it != perms.end(); ++it) {
    std::vector<int> indices = it->first;
    std::reverse( indices.begin(), indices.end() );
    trans_perms[ indices ] = it->second;
  }

  // Now bundle them togther
  for (auto it = trans_perms.begin(); it != trans_perms.end(); ++it) {
    perms[ it->first ] = it->second;
  }

  return perms;
}

//===========================================================================================================================================================

}

