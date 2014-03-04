/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012 Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "npdm_permutations.h"

namespace SpinAdapted{

//===========================================================================================================================================================

void Npdm_permutations::process_new_elements( const std::vector< std::pair< std::vector<int>, double > >& in, 
                                              std::vector< std::pair< std::vector<int>, double > >& nonredundant_elements,
                                              std::vector< std::pair< std::vector<int>, double > >& spin_perms )
{
  spin_perms.clear();
  int count = 0;
   
  // Loop over all input spin indices
  for (int i=0; i<in.size(); i++) {
    // Get all permutations of each set of spin indices
    std::vector< std::pair< std::vector<int>, double > > tmp;
    get_spin_permutations( tmp, in[i].first, in[i].second );
    // If permutations do not generate any of the following input indices, save the non-redundant original
    bool keep = true;
    for ( int j=0; j<tmp.size(); j++) {
      // Loop over remaining original indices
      for (int k=i+1; k<in.size(); k++) {
        if ( tmp[j].first == in[k].first ) {
          keep = false;
          break;
        }
      }
      if (!keep) break;
    }
    if (keep) {
       count++;
       nonredundant_elements.push_back( in[i] );
       spin_perms.insert( spin_perms.end(), tmp.begin(), tmp.end() );
    }
  }

  assert( count <= in.size() );
//cout << "nonredundant elements = " << count << endl;

}

//===========================================================================================================================================================

void Onepdm_permutations::get_spin_permutations( std::vector<std::pair<std::vector<int>,double> >& spin_batch, 
                                                 const std::vector<int>& indices, const double& val )
{
  assert( indices.size() == 2 );
  std::vector<int> idx;
  int i = indices[0];
  int j = indices[1];

  idx = { i, j };
  spin_batch.push_back( std::make_pair( idx, val ) );
  // Transpose is same 
  if ( i != j ) {
    idx = { j, i };
    spin_batch.push_back( std::make_pair( idx, val ) );
  }

}

//===========================================================================================================================================================

void Twopdm_permutations::get_spin_permutations( std::vector<std::pair<std::vector<int>,double> >& spin_batch, 
                                                 const std::vector<int>& indices, const double& val )
{
  assert( indices.size() == 4 );
  std::vector<int> idx;
  int i = indices[0];
  int j = indices[1];
  int k = indices[2];
  int l = indices[3];

  // If indices are not all unique, then all elements should be zero (and next_even_permutation fails)
  std::vector<int> v = {i,j};
  std::sort( v.begin(), v.end() );
  if ( (v[0]==v[1]) ) return;
  std::vector<int> w = {k,l};
  std::sort( w.begin(), w.end() );
  if ( (w[0]==w[1]) ) return;
  bool skip_transpose = ( v == w );

  // 8 permutations
  //--------------------------
  idx = { i, j, k, l }; spin_batch.push_back( std::make_pair( idx, val ) );
  idx = { i, j, l, k }; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx = { j, i, k, l }; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx = { j, i, l, k }; spin_batch.push_back( std::make_pair( idx, val ) );
                    
  if ( !skip_transpose ) {
    idx = { l, k, j, i }; spin_batch.push_back( std::make_pair( idx, val ) );
    idx = { k, l, j, i }; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx = { l, k, i, j }; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx = { k, l, i, j }; spin_batch.push_back( std::make_pair( idx, val ) );
  }

}

//===========================================================================================================================================================

void Threepdm_permutations::get_spin_permutations( std::vector<std::pair<std::vector<int>,double> >& spin_batch, 
                                                   const std::vector<int>& indices, const double& val )
{
  assert( indices.size() == 6 );
  std::vector<int> idx;
  int i = indices[0];
  int j = indices[1];
  int k = indices[2];
  int l = indices[3];
  int m = indices[4];
  int n = indices[5];

  // If indices are not all unique, then all elements should be zero (and next_even_permutation fails)
  std::vector<int> v = {i,j,k};
  std::sort( v.begin(), v.end() );
  if ( (v[0]==v[1]) || (v[1]==v[2]) ) return;
  std::vector<int> w = {l,m,n};
  std::sort( w.begin(), w.end() );
  if ( (w[0]==w[1]) || (w[1]==w[2]) ) return;
  bool skip_transpose = ( v == w );

  // The number of possible permutations is (3!)**2 twice
  idx = { i, j, k, l, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx = { i, j, k, l, n, m }; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx = { i, j, k, n, m, l }; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx = { i, j, k, n, l, m }; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx = { i, j, k, m, l, n }; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx = { i, j, k, m, n, l }; spin_batch.push_back( std::make_pair( idx,  val ) );

  idx = { i, k, j, l, m, n }; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx = { i, k, j, l, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx = { i, k, j, n, m, l }; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx = { i, k, j, n, l, m }; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx = { i, k, j, m, l, n }; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx = { i, k, j, m, n, l }; spin_batch.push_back( std::make_pair( idx, -val ) );

  idx = { j, i, k, l, m, n }; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx = { j, i, k, l, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx = { j, i, k, n, m, l }; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx = { j, i, k, n, l, m }; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx = { j, i, k, m, l, n }; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx = { j, i, k, m, n, l }; spin_batch.push_back( std::make_pair( idx, -val ) );

  idx = { j, k, i, l, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx = { j, k, i, l, n, m }; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx = { j, k, i, n, m, l }; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx = { j, k, i, n, l, m }; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx = { j, k, i, m, l, n }; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx = { j, k, i, m, n, l }; spin_batch.push_back( std::make_pair( idx,  val ) );

  idx = { k, j, i, l, m, n }; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx = { k, j, i, l, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx = { k, j, i, n, m, l }; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx = { k, j, i, n, l, m }; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx = { k, j, i, m, l, n }; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx = { k, j, i, m, n, l }; spin_batch.push_back( std::make_pair( idx, -val ) );

  idx = { k, i, j, l, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx = { k, i, j, l, n, m }; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx = { k, i, j, n, m, l }; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx = { k, i, j, n, l, m }; spin_batch.push_back( std::make_pair( idx,  val ) );
  idx = { k, i, j, m, l, n }; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx = { k, i, j, m, n, l }; spin_batch.push_back( std::make_pair( idx,  val ) );

  // Get transpose elements with same parity factors, hardcoded for speed
  if ( !skip_transpose ) {
    idx = { n, m, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx = { n, m, l, k, i, j }; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx = { n, m, l, i, j, k }; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx = { n, m, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx = { n, m, l, j, k, i }; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx = { n, m, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) );
  
    idx = { n, l, m, k, j, i }; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx = { n, l, m, k, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx = { n, l, m, i, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx = { n, l, m, i, k, j }; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx = { n, l, m, j, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx = { n, l, m, j, i, k }; spin_batch.push_back( std::make_pair( idx, -val ) );
  
    idx = { m, n, l, k, j, i }; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx = { m, n, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx = { m, n, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx = { m, n, l, i, k, j }; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx = { m, n, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx = { m, n, l, j, i, k }; spin_batch.push_back( std::make_pair( idx, -val ) );
  
    idx = { m, l, n, k, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx = { m, l, n, k, i, j }; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx = { m, l, n, i, j, k }; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx = { m, l, n, i, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx = { m, l, n, j, k, i }; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx = { m, l, n, j, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) );
  
    idx = { l, m, n, k, j, i }; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx = { l, m, n, k, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx = { l, m, n, i, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx = { l, m, n, i, k, j }; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx = { l, m, n, j, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx = { l, m, n, j, i, k }; spin_batch.push_back( std::make_pair( idx, -val ) );
  
    idx = { l, n, m, k, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx = { l, n, m, k, i, j }; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx = { l, n, m, i, j, k }; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx = { l, n, m, i, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) );
    idx = { l, n, m, j, k, i }; spin_batch.push_back( std::make_pair( idx, -val ) );
    idx = { l, n, m, j, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) );
  }

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

void Fourpdm_permutations::get_spin_permutations( std::vector<std::pair<std::vector<int>,double> >& spin_batch, 
                                                  const std::vector<int>& indices, const double& val )
{
  assert( indices.size() == 8 );
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
  if ( (v[0]==v[1]) || (v[1]==v[2]) || (v[2]==v[3]) ) return;
  std::vector<int> w = {m,n,p,q};
  std::sort( w.begin(), w.end() );
  if ( (w[0]==w[1]) || (w[1]==w[2]) || (w[2]==w[3]) ) return;
  bool skip_transpose = ( v == w );

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
      idx = { (*u)[0],(*u)[1],(*u)[2],(*u)[3], (*v)[0],(*v)[1],(*v)[2],(*v)[3] };
      spin_batch.push_back( std::make_pair( idx, val ) );
      if ( !skip_transpose ) {
        // Include transpose
        std::reverse( idx.begin(), idx.end() );
        spin_batch.push_back( std::make_pair( idx, val ) );
      }
    }
  }
  // Even-Odd terms
  for ( auto u = ijkl_even.begin(); u != ijkl_even.end(); ++u ) {
    for ( auto v = mnpq_odd.begin(); v != mnpq_odd.end(); ++v ) {
      idx = { (*u)[0],(*u)[1],(*u)[2],(*u)[3], (*v)[0],(*v)[1],(*v)[2],(*v)[3] };
      spin_batch.push_back( std::make_pair( idx, -val ) );
      if ( !skip_transpose ) {
        // Include transpose
        std::reverse( idx.begin(), idx.end() );
        spin_batch.push_back( std::make_pair( idx, -val ) );
      }
    }
  }
  // Odd-Even terms
  for ( auto u = ijkl_odd.begin(); u != ijkl_odd.end(); ++u ) {
    for ( auto v = mnpq_even.begin(); v != mnpq_even.end(); ++v ) {
      idx = { (*u)[0],(*u)[1],(*u)[2],(*u)[3], (*v)[0],(*v)[1],(*v)[2],(*v)[3] }; 
      spin_batch.push_back( std::make_pair( idx, -val ) );
      if ( !skip_transpose ) {
        // Include transpose
        std::reverse( idx.begin(), idx.end() );
        spin_batch.push_back( std::make_pair( idx, -val ) );
      }
    }
  }
  // Odd-Odd terms
  for ( auto u = ijkl_odd.begin(); u != ijkl_odd.end(); ++u ) {
    for ( auto v = mnpq_odd.begin(); v != mnpq_odd.end(); ++v ) {
      idx = { (*u)[0],(*u)[1],(*u)[2],(*u)[3], (*v)[0],(*v)[1],(*v)[2],(*v)[3] };
      spin_batch.push_back( std::make_pair( idx, val ) );
      if ( !skip_transpose ) {
        // Include transpose
        std::reverse( idx.begin(), idx.end() );
        spin_batch.push_back( std::make_pair( idx, val ) );
      }
    }
  }

}

//===========================================================================================================================================================

}

