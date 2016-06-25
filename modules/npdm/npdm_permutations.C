/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012 Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include <cassert>
#include <algorithm>

#include "npdm_permutations.h"
#include "spinblock.h"
#include "pario.h"

namespace SpinAdapted{

//===========================================================================================================================================================

void Npdm_permutations::get_permute(const std::pair<std::vector<int>,int>& origin, int start, int n, std::vector<std::pair<std::vector<int>,int> >& reorders)
{
  if (n<2 ||start+1 >= n){
    reorders.push_back(origin);
    return;
  }
  get_permute(origin,start+1,n,reorders);
  for(int i=start+1; i<n;i++)
  {
    std::pair<std::vector<int>,int> neworder = origin;
    int tmp = neworder.first[i];
    neworder.first[i]= neworder.first[start];
    neworder.first[start] = tmp;
    neworder.second *= -1;
    get_permute(neworder,start+1,n,reorders);
  }


}


void Onepdm_permutations::get_spatial_batch( const std::vector< std::pair< std::vector<int>, double > >& in, 
                                              std::vector< std::pair< std::vector<int>, double > >& spatial_perms )
{
  spatial_perms.clear();
  std::vector<int> spatial(2);

  if (in.size()==0)
    return ;
  spatial[0] = in[0].first[0]/2;
  spatial[1] = in[0].first[1]/2;
  double value = in[0].second + in[1].second;
  spatial_perms.push_back(std::make_pair(spatial,value));
  if (spatial[0] !=spatial[1] && dmrginp.doimplicitTranspose())
    {
      int tmp = spatial[0];
      spatial[0] = spatial[1];
      spatial[1] = tmp;
      spatial_perms.push_back(std::make_pair(spatial,value));
    }
}

//===========================================================================================================================================================

void Twopdm_permutations::get_spatial_batch( const std::vector< std::pair< std::vector<int>, double > >& in, 
                                              std::vector< std::pair< std::vector<int>, double > >& spatial_perms )
{
  spatial_perms.clear();
  std::vector<int> indices;
  std::vector<std::vector<int>> newindices;
  std::vector<int> tmp;

  if (in.size()==0)
    return ;

  for(int i=0;i<in[0].first.size();i++)
    indices.push_back(in[0].first.at(i)/2);

  for(int i=0;i<reorders.size();i++)
  {
    
    std::vector<int> spatial_indices = indices;
    for(int j=0; j< reorders[i].first.size();j++)
    {
      spatial_indices[j] = indices[reorders[i].first[j]];
    }
    bool skip = false;
    for(int j=0; j< spatial_perms.size();j++)
    {
      if (spatial_perms[j].first == spatial_indices)
      {
        skip = true;
        break;
      }
    }
    if(!skip)
    {
      double value = 0;
      for(int j=0; j<in.size();j++)
      {
        if (in[j].first.at(reorders[i].first[0])%2 == in[j].first.at(3)%2 &&
            in[j].first.at(reorders[i].first[1])%2 == in[j].first.at(2)%2 )
          value += in[j].second*reorders[i].second;
      }
      if(abs(value)< NUMERICAL_ZERO)
        continue;
      std::set<std::vector<int> > spatial_batch;
      std::vector<int> tmp(4);
      for(int k=0;k<reorders.size();k++)
      {
        tmp[0]= spatial_indices[reorders[k].first[0]];
        tmp[1]= spatial_indices[reorders[k].first[1]];
        tmp[2]= spatial_indices[3-reorders[k].first[1]];
        tmp[3]= spatial_indices[3-reorders[k].first[0]];
        spatial_batch.insert(tmp);
        if(dmrginp.doimplicitTranspose())
        {
          tmp[0]= spatial_indices[3-reorders[k].first[0]];
          tmp[1]= spatial_indices[3-reorders[k].first[1]];
          tmp[2]= spatial_indices[reorders[k].first[1]];
          tmp[3]= spatial_indices[reorders[k].first[0]];
          spatial_batch.insert(tmp);
        }
      }
      for(auto x: spatial_batch)
        spatial_perms.push_back(std::make_pair(x,value));
    }
  }
}

//===========================================================================================================================================================

void Threepdm_permutations::get_spatial_batch( const std::vector< std::pair< std::vector<int>, double > >& in, 
                                              std::vector< std::pair< std::vector<int>, double > >& spatial_perms )
{
  spatial_perms.clear();
  std::vector<int> indices;
  std::vector<std::vector<int>> newindices;
  std::vector<int> tmp;

  if (in.size()==0)
    return ;

  for(int i=0;i<in[0].first.size();i++)
    indices.push_back(in[0].first.at(i)/2);

  for(int i=0;i<reorders.size();i++)
  {
    
    std::vector<int> spatial_indices = indices;
    for(int j=0; j< reorders[i].first.size();j++)
    {
      spatial_indices[j] = indices[reorders[i].first[j]];
    }
    bool skip = false;
    for(int j=0; j< spatial_perms.size();j++)
    {
      if (spatial_perms[j].first == spatial_indices)
      {
        skip = true;
        break;
      }
    }
    if(!skip)
    {
      double value = 0;
      for(int j=0; j<in.size();j++)
      {
        if (in[j].first.at(reorders[i].first[0])%2 == in[j].first.at(5)%2 &&
            in[j].first.at(reorders[i].first[1])%2 == in[j].first.at(4)%2 &&
            in[j].first.at(reorders[i].first[2])%2 == in[j].first.at(3)%2 )
          value += in[j].second*reorders[i].second;
      }
      if(abs(value)< NUMERICAL_ZERO)
        continue;
      std::set<std::vector<int> > spatial_batch;
      std::vector<int> tmp(6);
      for(int k=0;k<reorders.size();k++)
      {
        tmp[0]= spatial_indices[reorders[k].first[0]];
        tmp[1]= spatial_indices[reorders[k].first[1]];
        tmp[2]= spatial_indices[reorders[k].first[2]];
        tmp[3]= spatial_indices[5-reorders[k].first[2]];
        tmp[4]= spatial_indices[5-reorders[k].first[1]];
        tmp[5]= spatial_indices[5-reorders[k].first[0]];
        spatial_batch.insert(tmp);
        if(dmrginp.doimplicitTranspose())
        {
          tmp[0]= spatial_indices[5-reorders[k].first[0]];
          tmp[1]= spatial_indices[5-reorders[k].first[1]];
          tmp[2]= spatial_indices[5-reorders[k].first[2]];
          tmp[3]= spatial_indices[reorders[k].first[2]];
          tmp[4]= spatial_indices[reorders[k].first[1]];
          tmp[5]= spatial_indices[reorders[k].first[0]];
          spatial_batch.insert(tmp);
        }
      }
   //   int k = spatial_indices[0];
   //   int l = spatial_indices[1];
   //   int m = spatial_indices[2];
   //   int n = spatial_indices[3];
   //   int p = spatial_indices[4];
   //   int q = spatial_indices[5];
   //   tmp={k,l,m,n,p,q}; spatial_batch.insert(tmp);
   //   tmp={k,m,l,p,n,q}; spatial_batch.insert(tmp);
   //   tmp={l,k,m,n,q,p}; spatial_batch.insert(tmp);
   //   tmp={l,m,k,q,n,p}; spatial_batch.insert(tmp);
   //   tmp={m,k,l,p,q,n}; spatial_batch.insert(tmp);
   //   tmp={m,l,k,q,p,n}; spatial_batch.insert(tmp);
   //   if(dmrginp.doimplicitTranspose())
   //   {
   //     tmp={q,p,n,m,l,k}; spatial_batch.insert(tmp);
   //     tmp={q,n,p,l,m,k}; spatial_batch.insert(tmp);
   //     tmp={p,q,n,m,k,l}; spatial_batch.insert(tmp);
   //     tmp={p,n,q,k,m,l}; spatial_batch.insert(tmp);
   //     tmp={n,q,p,l,k,m}; spatial_batch.insert(tmp);
   //     tmp={n,p,q,k,l,m}; spatial_batch.insert(tmp);
   //   }
      for(auto x: spatial_batch)
        spatial_perms.push_back(std::make_pair(x,value));
    }
  }
}

//===========================================================================================================================================================

void Fourpdm_permutations::get_spatial_batch( const std::vector< std::pair< std::vector<int>, double > >& in, 
                                              std::vector< std::pair< std::vector<int>, double > >& spatial_perms )
{
  spatial_perms.clear();
  std::vector<int> indices;
  std::vector<std::vector<int>> newindices;
  std::vector<int> tmp;

  if (in.size()==0)
    return ;

  for(int i=0;i<in[0].first.size();i++)
    indices.push_back(in[0].first.at(i)/2);

  for(int i=0;i<reorders.size();i++)
  {
    
    std::vector<int> spatial_indices = indices;
    for(int j=0; j< reorders[i].first.size();j++)
    {
      spatial_indices[j] = indices[reorders[i].first[j]];
    }
    bool skip = false;
    for(int j=0; j< spatial_perms.size();j++)
    {
      if (spatial_perms[j].first == spatial_indices)
      {
        skip = true;
        break;
      }
    }
    if(!skip)
    {
      double value = 0;
      for(int j=0; j<in.size();j++)
      {
        if (in[j].first.at(reorders[i].first[0])%2 == in[j].first.at(7)%2 &&
            in[j].first.at(reorders[i].first[1])%2 == in[j].first.at(6)%2 &&
            in[j].first.at(reorders[i].first[2])%2 == in[j].first.at(5)%2 &&
            in[j].first.at(reorders[i].first[3])%2 == in[j].first.at(4)%2 )
          value += in[j].second*reorders[i].second;
      }
      if(abs(value)< NUMERICAL_ZERO)
        continue;
      std::vector<int> tmp(8);
      std::set<std::vector<int> > spatial_batch;
      for(int k=0;k<reorders.size();k++)
      {
        tmp[0]= spatial_indices[reorders[k].first[0]];
        tmp[1]= spatial_indices[reorders[k].first[1]];
        tmp[2]= spatial_indices[reorders[k].first[2]];
        tmp[3]= spatial_indices[reorders[k].first[3]];
        tmp[4]= spatial_indices[7-reorders[k].first[3]];
        tmp[5]= spatial_indices[7-reorders[k].first[2]];
        tmp[6]= spatial_indices[7-reorders[k].first[1]];
        tmp[7]= spatial_indices[7-reorders[k].first[0]];
        spatial_batch.insert(tmp);
        if(dmrginp.doimplicitTranspose())
        {
          tmp[0]= spatial_indices[7-reorders[k].first[0]];
          tmp[1]= spatial_indices[7-reorders[k].first[1]];
          tmp[2]= spatial_indices[7-reorders[k].first[2]];
          tmp[3]= spatial_indices[7-reorders[k].first[3]];
          tmp[4]= spatial_indices[reorders[k].first[3]];
          tmp[5]= spatial_indices[reorders[k].first[2]];
          tmp[6]= spatial_indices[reorders[k].first[1]];
          tmp[7]= spatial_indices[reorders[k].first[0]];
          spatial_batch.insert(tmp);
        }
      }
      for(auto x: spatial_batch)
        spatial_perms.push_back(std::make_pair(x,value));
    }
  }
}

//===========================================================================================================================================================

void Npdm_permutations::process_new_elements( const std::vector< std::pair< std::vector<int>, double > >& in, 
                                              std::vector< std::pair< std::vector<int>, double > >& nonredundant_elements,
                                              std::vector< std::pair< std::vector<int>, double > >& spin_perms )
{
  spin_perms.clear();
  int count = 0;
   
  std::vector<int> spatial_indices;
  if (in.size() != 0)
  {
    for(int i=0;i<in[0].first.size();i++)
      spatial_indices.push_back(in[0].first.at(i)/2);

  }
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
//pout << "nonredundant elements = " << count << endl;

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
  if(dmrginp.doimplicitTranspose())
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
  if ( v[0]==v[1] ) return;
  std::vector<int> w = {k,l};
  std::sort( w.begin(), w.end() );
  if ( w[0]==w[1] ) return;
  bool skip_transpose = ( v == w );

  // 8 permutations
  //--------------------------
  idx = { i, j, k, l }; spin_batch.push_back( std::make_pair( idx, val ) );
  idx = { i, j, l, k }; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx = { j, i, k, l }; spin_batch.push_back( std::make_pair( idx, -val ) );
  idx = { j, i, l, k }; spin_batch.push_back( std::make_pair( idx, val ) );
                    
  if ( !skip_transpose && dmrginp.doimplicitTranspose()) {
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
// This is a general routine for arbitrary order permutations.
/*
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
// This is a general routine for arbitrary order permutations.

void Fourpdm_permutations::get_spin_permutations_general( std::vector<std::pair<std::vector<int>,double> >& spin_batch, 
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
*/

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//FIXME can speed up by replacing push_backs with [] ? Or .reserve() ?
// Note we can use the arbitary order routines above instead, but probably slower.

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
  idx = { i, j, k, l, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, j, k, l, m, p, q, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, j, k, l, m, q, n, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, j, k, l, n, m, q, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, j, k, l, n, p, m, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, j, k, l, n, q, p, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, j, k, l, p, m, n, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, j, k, l, p, n, q, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, j, k, l, p, q, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, j, k, l, q, m, p, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, j, k, l, q, n, m, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, j, k, l, q, p, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, k, l, j, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, k, l, j, m, p, q, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, k, l, j, m, q, n, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, k, l, j, n, m, q, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, k, l, j, n, p, m, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, k, l, j, n, q, p, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, k, l, j, p, m, n, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, k, l, j, p, n, q, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, k, l, j, p, q, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, k, l, j, q, m, p, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, k, l, j, q, n, m, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, k, l, j, q, p, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, l, j, k, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, l, j, k, m, p, q, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, l, j, k, m, q, n, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, l, j, k, n, m, q, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, l, j, k, n, p, m, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, l, j, k, n, q, p, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, l, j, k, p, m, n, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, l, j, k, p, n, q, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, l, j, k, p, q, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, l, j, k, q, m, p, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, l, j, k, q, n, m, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, l, j, k, q, p, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, i, l, k, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, i, l, k, m, p, q, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, i, l, k, m, q, n, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, i, l, k, n, m, q, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, i, l, k, n, p, m, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, i, l, k, n, q, p, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, i, l, k, p, m, n, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, i, l, k, p, n, q, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, i, l, k, p, q, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, i, l, k, q, m, p, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, i, l, k, q, n, m, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, i, l, k, q, p, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, k, i, l, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, k, i, l, m, p, q, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, k, i, l, m, q, n, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, k, i, l, n, m, q, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, k, i, l, n, p, m, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, k, i, l, n, q, p, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, k, i, l, p, m, n, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, k, i, l, p, n, q, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, k, i, l, p, q, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, k, i, l, q, m, p, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, k, i, l, q, n, m, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, k, i, l, q, p, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, l, k, i, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, l, k, i, m, p, q, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, l, k, i, m, q, n, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, l, k, i, n, m, q, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, l, k, i, n, p, m, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, l, k, i, n, q, p, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, l, k, i, p, m, n, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, l, k, i, p, n, q, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, l, k, i, p, q, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, l, k, i, q, m, p, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, l, k, i, q, n, m, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, l, k, i, q, p, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, i, j, l, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, i, j, l, m, p, q, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, i, j, l, m, q, n, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, i, j, l, n, m, q, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, i, j, l, n, p, m, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, i, j, l, n, q, p, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, i, j, l, p, m, n, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, i, j, l, p, n, q, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, i, j, l, p, q, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, i, j, l, q, m, p, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, i, j, l, q, n, m, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, i, j, l, q, p, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, j, l, i, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, j, l, i, m, p, q, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, j, l, i, m, q, n, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, j, l, i, n, m, q, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, j, l, i, n, p, m, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, j, l, i, n, q, p, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, j, l, i, p, m, n, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, j, l, i, p, n, q, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, j, l, i, p, q, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, j, l, i, q, m, p, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, j, l, i, q, n, m, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, j, l, i, q, p, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, l, i, j, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, l, i, j, m, p, q, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, l, i, j, m, q, n, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, l, i, j, n, m, q, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, l, i, j, n, p, m, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, l, i, j, n, q, p, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, l, i, j, p, m, n, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, l, i, j, p, n, q, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, l, i, j, p, q, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, l, i, j, q, m, p, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, l, i, j, q, n, m, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, l, i, j, q, p, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, i, k, j, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, i, k, j, m, p, q, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, i, k, j, m, q, n, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, i, k, j, n, m, q, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, i, k, j, n, p, m, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, i, k, j, n, q, p, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, i, k, j, p, m, n, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, i, k, j, p, n, q, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, i, k, j, p, q, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, i, k, j, q, m, p, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, i, k, j, q, n, m, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, i, k, j, q, p, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, j, i, k, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, j, i, k, m, p, q, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, j, i, k, m, q, n, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, j, i, k, n, m, q, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, j, i, k, n, p, m, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, j, i, k, n, q, p, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, j, i, k, p, m, n, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, j, i, k, p, n, q, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, j, i, k, p, q, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, j, i, k, q, m, p, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, j, i, k, q, n, m, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, j, i, k, q, p, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, k, j, i, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, k, j, i, m, p, q, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, k, j, i, m, q, n, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, k, j, i, n, m, q, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, k, j, i, n, p, m, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, k, j, i, n, q, p, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, k, j, i, p, m, n, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, k, j, i, p, n, q, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, k, j, i, p, q, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, k, j, i, q, m, p, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, k, j, i, q, n, m, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, k, j, i, q, p, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, j, k, l, m, n, q, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, j, k, l, m, p, n, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, j, k, l, m, q, p, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, j, k, l, n, m, p, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, j, k, l, n, p, q, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, j, k, l, n, q, m, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, j, k, l, p, m, q, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, j, k, l, p, n, m, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, j, k, l, p, q, n, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, j, k, l, q, m, n, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, j, k, l, q, n, p, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, j, k, l, q, p, m, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, k, l, j, m, n, q, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, k, l, j, m, p, n, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, k, l, j, m, q, p, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, k, l, j, n, m, p, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, k, l, j, n, p, q, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, k, l, j, n, q, m, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, k, l, j, p, m, q, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, k, l, j, p, n, m, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, k, l, j, p, q, n, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, k, l, j, q, m, n, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, k, l, j, q, n, p, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, k, l, j, q, p, m, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, l, j, k, m, n, q, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, l, j, k, m, p, n, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, l, j, k, m, q, p, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, l, j, k, n, m, p, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, l, j, k, n, p, q, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, l, j, k, n, q, m, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, l, j, k, p, m, q, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, l, j, k, p, n, m, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, l, j, k, p, q, n, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, l, j, k, q, m, n, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, l, j, k, q, n, p, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, l, j, k, q, p, m, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, i, l, k, m, n, q, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, i, l, k, m, p, n, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, i, l, k, m, q, p, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, i, l, k, n, m, p, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, i, l, k, n, p, q, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, i, l, k, n, q, m, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, i, l, k, p, m, q, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, i, l, k, p, n, m, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, i, l, k, p, q, n, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, i, l, k, q, m, n, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, i, l, k, q, n, p, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, i, l, k, q, p, m, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, k, i, l, m, n, q, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, k, i, l, m, p, n, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, k, i, l, m, q, p, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, k, i, l, n, m, p, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, k, i, l, n, p, q, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, k, i, l, n, q, m, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, k, i, l, p, m, q, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, k, i, l, p, n, m, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, k, i, l, p, q, n, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, k, i, l, q, m, n, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, k, i, l, q, n, p, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, k, i, l, q, p, m, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, l, k, i, m, n, q, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, l, k, i, m, p, n, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, l, k, i, m, q, p, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, l, k, i, n, m, p, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, l, k, i, n, p, q, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, l, k, i, n, q, m, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, l, k, i, p, m, q, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, l, k, i, p, n, m, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, l, k, i, p, q, n, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, l, k, i, q, m, n, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, l, k, i, q, n, p, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, l, k, i, q, p, m, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, i, j, l, m, n, q, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, i, j, l, m, p, n, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, i, j, l, m, q, p, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, i, j, l, n, m, p, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, i, j, l, n, p, q, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, i, j, l, n, q, m, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, i, j, l, p, m, q, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, i, j, l, p, n, m, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, i, j, l, p, q, n, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, i, j, l, q, m, n, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, i, j, l, q, n, p, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, i, j, l, q, p, m, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, j, l, i, m, n, q, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, j, l, i, m, p, n, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, j, l, i, m, q, p, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, j, l, i, n, m, p, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, j, l, i, n, p, q, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, j, l, i, n, q, m, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, j, l, i, p, m, q, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, j, l, i, p, n, m, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, j, l, i, p, q, n, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, j, l, i, q, m, n, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, j, l, i, q, n, p, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, j, l, i, q, p, m, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, l, i, j, m, n, q, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, l, i, j, m, p, n, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, l, i, j, m, q, p, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, l, i, j, n, m, p, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, l, i, j, n, p, q, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, l, i, j, n, q, m, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, l, i, j, p, m, q, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, l, i, j, p, n, m, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, l, i, j, p, q, n, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, l, i, j, q, m, n, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, l, i, j, q, n, p, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, l, i, j, q, p, m, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, i, k, j, m, n, q, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, i, k, j, m, p, n, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, i, k, j, m, q, p, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, i, k, j, n, m, p, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, i, k, j, n, p, q, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, i, k, j, n, q, m, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, i, k, j, p, m, q, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, i, k, j, p, n, m, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, i, k, j, p, q, n, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, i, k, j, q, m, n, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, i, k, j, q, n, p, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, i, k, j, q, p, m, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, j, i, k, m, n, q, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, j, i, k, m, p, n, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, j, i, k, m, q, p, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, j, i, k, n, m, p, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, j, i, k, n, p, q, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, j, i, k, n, q, m, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, j, i, k, p, m, q, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, j, i, k, p, n, m, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, j, i, k, p, q, n, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, j, i, k, q, m, n, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, j, i, k, q, n, p, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, j, i, k, q, p, m, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, k, j, i, m, n, q, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, k, j, i, m, p, n, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, k, j, i, m, q, p, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, k, j, i, n, m, p, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, k, j, i, n, p, q, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, k, j, i, n, q, m, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, k, j, i, p, m, q, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, k, j, i, p, n, m, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, k, j, i, p, q, n, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, k, j, i, q, m, n, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, k, j, i, q, n, p, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, k, j, i, q, p, m, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, j, l, k, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, j, l, k, m, p, q, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, j, l, k, m, q, n, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, j, l, k, n, m, q, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, j, l, k, n, p, m, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, j, l, k, n, q, p, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, j, l, k, p, m, n, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, j, l, k, p, n, q, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, j, l, k, p, q, m, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, j, l, k, q, m, p, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, j, l, k, q, n, m, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, j, l, k, q, p, n, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, k, j, l, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, k, j, l, m, p, q, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, k, j, l, m, q, n, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, k, j, l, n, m, q, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, k, j, l, n, p, m, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, k, j, l, n, q, p, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, k, j, l, p, m, n, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, k, j, l, p, n, q, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, k, j, l, p, q, m, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, k, j, l, q, m, p, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, k, j, l, q, n, m, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, k, j, l, q, p, n, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, l, k, j, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, l, k, j, m, p, q, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, l, k, j, m, q, n, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, l, k, j, n, m, q, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, l, k, j, n, p, m, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, l, k, j, n, q, p, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, l, k, j, p, m, n, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, l, k, j, p, n, q, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, l, k, j, p, q, m, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, l, k, j, q, m, p, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, l, k, j, q, n, m, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, l, k, j, q, p, n, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, i, k, l, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, i, k, l, m, p, q, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, i, k, l, m, q, n, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, i, k, l, n, m, q, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, i, k, l, n, p, m, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, i, k, l, n, q, p, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, i, k, l, p, m, n, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, i, k, l, p, n, q, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, i, k, l, p, q, m, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, i, k, l, q, m, p, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, i, k, l, q, n, m, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, i, k, l, q, p, n, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, k, l, i, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, k, l, i, m, p, q, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, k, l, i, m, q, n, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, k, l, i, n, m, q, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, k, l, i, n, p, m, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, k, l, i, n, q, p, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, k, l, i, p, m, n, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, k, l, i, p, n, q, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, k, l, i, p, q, m, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, k, l, i, q, m, p, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, k, l, i, q, n, m, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, k, l, i, q, p, n, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, l, i, k, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, l, i, k, m, p, q, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, l, i, k, m, q, n, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, l, i, k, n, m, q, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, l, i, k, n, p, m, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, l, i, k, n, q, p, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, l, i, k, p, m, n, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, l, i, k, p, n, q, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, l, i, k, p, q, m, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, l, i, k, q, m, p, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, l, i, k, q, n, m, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { j, l, i, k, q, p, n, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, i, l, j, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, i, l, j, m, p, q, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, i, l, j, m, q, n, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, i, l, j, n, m, q, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, i, l, j, n, p, m, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, i, l, j, n, q, p, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, i, l, j, p, m, n, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, i, l, j, p, n, q, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, i, l, j, p, q, m, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, i, l, j, q, m, p, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, i, l, j, q, n, m, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, i, l, j, q, p, n, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, j, i, l, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, j, i, l, m, p, q, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, j, i, l, m, q, n, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, j, i, l, n, m, q, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, j, i, l, n, p, m, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); idx = { k, j, i, l, n, q, p, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); idx = { k, j, i, l, p, m, n, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); idx = { k, j, i, l, p, n, q, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, j, i, l, p, q, m, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, j, i, l, q, m, p, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, j, i, l, q, n, m, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, j, i, l, q, p, n, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, l, j, i, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, l, j, i, m, p, q, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, l, j, i, m, q, n, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, l, j, i, n, m, q, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, l, j, i, n, p, m, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, l, j, i, n, q, p, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, l, j, i, p, m, n, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, l, j, i, p, n, q, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, l, j, i, p, q, m, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, l, j, i, q, m, p, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, l, j, i, q, n, m, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { k, l, j, i, q, p, n, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, i, j, k, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, i, j, k, m, p, q, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, i, j, k, m, q, n, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, i, j, k, n, m, q, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, i, j, k, n, p, m, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, i, j, k, n, q, p, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, i, j, k, p, m, n, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, i, j, k, p, n, q, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, i, j, k, p, q, m, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, i, j, k, q, m, p, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, i, j, k, q, n, m, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, i, j, k, q, p, n, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, j, k, i, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, j, k, i, m, p, q, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, j, k, i, m, q, n, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, j, k, i, n, m, q, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, j, k, i, n, p, m, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, j, k, i, n, q, p, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, j, k, i, p, m, n, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, j, k, i, p, n, q, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, j, k, i, p, q, m, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, j, k, i, q, m, p, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, j, k, i, q, n, m, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, j, k, i, q, p, n, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, k, i, j, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, k, i, j, m, p, q, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, k, i, j, m, q, n, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, k, i, j, n, m, q, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, k, i, j, n, p, m, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, k, i, j, n, q, p, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, k, i, j, p, m, n, q }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, k, i, j, p, n, q, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, k, i, j, p, q, m, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, k, i, j, q, m, p, n }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, k, i, j, q, n, m, p }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { l, k, i, j, q, p, n, m }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
  idx = { i, j, l, k, m, n, q, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, j, l, k, m, p, n, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, j, l, k, m, q, p, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, j, l, k, n, m, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, j, l, k, n, p, q, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, j, l, k, n, q, m, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, j, l, k, p, m, q, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, j, l, k, p, n, m, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, j, l, k, p, q, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, j, l, k, q, m, n, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, j, l, k, q, n, p, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, j, l, k, q, p, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, k, j, l, m, n, q, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, k, j, l, m, p, n, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, k, j, l, m, q, p, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, k, j, l, n, m, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, k, j, l, n, p, q, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, k, j, l, n, q, m, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, k, j, l, p, m, q, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, k, j, l, p, n, m, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, k, j, l, p, q, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, k, j, l, q, m, n, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, k, j, l, q, n, p, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, k, j, l, q, p, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, l, k, j, m, n, q, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, l, k, j, m, p, n, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, l, k, j, m, q, p, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, l, k, j, n, m, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, l, k, j, n, p, q, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, l, k, j, n, q, m, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, l, k, j, p, m, q, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, l, k, j, p, n, m, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, l, k, j, p, q, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, l, k, j, q, m, n, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, l, k, j, q, n, p, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { i, l, k, j, q, p, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, i, k, l, m, n, q, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, i, k, l, m, p, n, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, i, k, l, m, q, p, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, i, k, l, n, m, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, i, k, l, n, p, q, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, i, k, l, n, q, m, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, i, k, l, p, m, q, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, i, k, l, p, n, m, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, i, k, l, p, q, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, i, k, l, q, m, n, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, i, k, l, q, n, p, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, i, k, l, q, p, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, k, l, i, m, n, q, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, k, l, i, m, p, n, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, k, l, i, m, q, p, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, k, l, i, n, m, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, k, l, i, n, p, q, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, k, l, i, n, q, m, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, k, l, i, p, m, q, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, k, l, i, p, n, m, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, k, l, i, p, q, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, k, l, i, q, m, n, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, k, l, i, q, n, p, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, k, l, i, q, p, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, l, i, k, m, n, q, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, l, i, k, m, p, n, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, l, i, k, m, q, p, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, l, i, k, n, m, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, l, i, k, n, p, q, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, l, i, k, n, q, m, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, l, i, k, p, m, q, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, l, i, k, p, n, m, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, l, i, k, p, q, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, l, i, k, q, m, n, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, l, i, k, q, n, p, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { j, l, i, k, q, p, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, i, l, j, m, n, q, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, i, l, j, m, p, n, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, i, l, j, m, q, p, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, i, l, j, n, m, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, i, l, j, n, p, q, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, i, l, j, n, q, m, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, i, l, j, p, m, q, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, i, l, j, p, n, m, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, i, l, j, p, q, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, i, l, j, q, m, n, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, i, l, j, q, n, p, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, i, l, j, q, p, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, j, i, l, m, n, q, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, j, i, l, m, p, n, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, j, i, l, m, q, p, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, j, i, l, n, m, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, j, i, l, n, p, q, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, j, i, l, n, q, m, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, j, i, l, p, m, q, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, j, i, l, p, n, m, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, j, i, l, p, q, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, j, i, l, q, m, n, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, j, i, l, q, n, p, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, j, i, l, q, p, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, l, j, i, m, n, q, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, l, j, i, m, p, n, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, l, j, i, m, q, p, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, l, j, i, n, m, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, l, j, i, n, p, q, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, l, j, i, n, q, m, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, l, j, i, p, m, q, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, l, j, i, p, n, m, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, l, j, i, p, q, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, l, j, i, q, m, n, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, l, j, i, q, n, p, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { k, l, j, i, q, p, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, i, j, k, m, n, q, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, i, j, k, m, p, n, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, i, j, k, m, q, p, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, i, j, k, n, m, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, i, j, k, n, p, q, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, i, j, k, n, q, m, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, i, j, k, p, m, q, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, i, j, k, p, n, m, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, i, j, k, p, q, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, i, j, k, q, m, n, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, i, j, k, q, n, p, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, i, j, k, q, p, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, j, k, i, m, n, q, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, j, k, i, m, p, n, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, j, k, i, m, q, p, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, j, k, i, n, m, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, j, k, i, n, p, q, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, j, k, i, n, q, m, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, j, k, i, p, m, q, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, j, k, i, p, n, m, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, j, k, i, p, q, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, j, k, i, q, m, n, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, j, k, i, q, n, p, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, j, k, i, q, p, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, k, i, j, m, n, q, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, k, i, j, m, p, n, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, k, i, j, m, q, p, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, k, i, j, n, m, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, k, i, j, n, p, q, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, k, i, j, n, q, m, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, k, i, j, p, m, q, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, k, i, j, p, n, m, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, k, i, j, p, q, n, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, k, i, j, q, m, n, p }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, k, i, j, q, n, p, m }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  idx = { l, k, i, j, q, p, m, n }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  // Get transpose elements with same parity factors, hardcoded for speed
  if ( !skip_transpose ) {
    idx = { q, p, n, m, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, q, p, m, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, n, q, m, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, q, m, n, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, m, p, n, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, p, q, n, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, n, m, p, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, q, n, p, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, m, q, p, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, p, m, q, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, m, n, q, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, n, p, q, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, p, n, m, j, l, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, q, p, m, j, l, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, n, q, m, j, l, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, q, m, n, j, l, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, m, p, n, j, l, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, p, q, n, j, l, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, n, m, p, j, l, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, q, n, p, j, l, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, m, q, p, j, l, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, p, m, q, j, l, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, m, n, q, j, l, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, n, p, q, j, l, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, p, n, m, k, j, l, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, q, p, m, k, j, l, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, n, q, m, k, j, l, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, q, m, n, k, j, l, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, m, p, n, k, j, l, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, p, q, n, k, j, l, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, n, m, p, k, j, l, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, q, n, p, k, j, l, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, m, q, p, k, j, l, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, p, m, q, k, j, l, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, m, n, q, k, j, l, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, n, p, q, k, j, l, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, p, n, m, k, l, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, q, p, m, k, l, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, n, q, m, k, l, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, q, m, n, k, l, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, m, p, n, k, l, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, p, q, n, k, l, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, n, m, p, k, l, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, q, n, p, k, l, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, m, q, p, k, l, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, p, m, q, k, l, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, m, n, q, k, l, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, n, p, q, k, l, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, p, n, m, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, q, p, m, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, n, q, m, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, q, m, n, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, m, p, n, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, p, q, n, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, n, m, p, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, q, n, p, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, m, q, p, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, p, m, q, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, m, n, q, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, n, p, q, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, p, n, m, i, k, l, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, q, p, m, i, k, l, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, n, q, m, i, k, l, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, q, m, n, i, k, l, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, m, p, n, i, k, l, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, p, q, n, i, k, l, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, n, m, p, i, k, l, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, q, n, p, i, k, l, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, m, q, p, i, k, l, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, p, m, q, i, k, l, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, m, n, q, i, k, l, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, n, p, q, i, k, l, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, p, n, m, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, q, p, m, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, n, q, m, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, q, m, n, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, m, p, n, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, p, q, n, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, n, m, p, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, q, n, p, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, m, q, p, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, p, m, q, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, m, n, q, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, n, p, q, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, p, n, m, i, l, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, q, p, m, i, l, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, n, q, m, i, l, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, q, m, n, i, l, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, m, p, n, i, l, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, p, q, n, i, l, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, n, m, p, i, l, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, q, n, p, i, l, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, m, q, p, i, l, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, p, m, q, i, l, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, m, n, q, i, l, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, n, p, q, i, l, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, p, n, m, j, i, l, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, q, p, m, j, i, l, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, n, q, m, j, i, l, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, q, m, n, j, i, l, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, m, p, n, j, i, l, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, p, q, n, j, i, l, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, n, m, p, j, i, l, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, q, n, p, j, i, l, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, m, q, p, j, i, l, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, p, m, q, j, i, l, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, m, n, q, j, i, l, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, n, p, q, j, i, l, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, p, n, m, j, k, i, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, q, p, m, j, k, i, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, n, q, m, j, k, i, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, q, m, n, j, k, i, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, m, p, n, j, k, i, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, p, q, n, j, k, i, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, n, m, p, j, k, i, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, q, n, p, j, k, i, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, m, q, p, j, k, i, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, p, m, q, j, k, i, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, m, n, q, j, k, i, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, n, p, q, j, k, i, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, p, n, m, k, i, j, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, q, p, m, k, i, j, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, n, q, m, k, i, j, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, q, m, n, k, i, j, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, m, p, n, k, i, j, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, p, q, n, k, i, j, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, n, m, p, k, i, j, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, q, n, p, k, i, j, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, m, q, p, k, i, j, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, p, m, q, k, i, j, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, m, n, q, k, i, j, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, n, p, q, k, i, j, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, p, n, m, i, j, k, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, q, p, m, i, j, k, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, n, q, m, i, j, k, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, q, m, n, i, j, k, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, m, p, n, i, j, k, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, p, q, n, i, j, k, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, n, m, p, i, j, k, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, q, n, p, i, j, k, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, m, q, p, i, j, k, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, p, m, q, i, j, k, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, m, n, q, i, j, k, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, n, p, q, i, j, k, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, q, n, m, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, n, p, m, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, p, q, m, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, p, m, n, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, q, p, n, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, m, q, n, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, q, m, p, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, m, n, p, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, n, q, p, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, n, m, q, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, p, n, q, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, m, p, q, l, k, j, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, q, n, m, j, l, k, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, n, p, m, j, l, k, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, p, q, m, j, l, k, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, p, m, n, j, l, k, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, q, p, n, j, l, k, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, m, q, n, j, l, k, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, q, m, p, j, l, k, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, m, n, p, j, l, k, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, n, q, p, j, l, k, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, n, m, q, j, l, k, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, p, n, q, j, l, k, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, m, p, q, j, l, k, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, q, n, m, k, j, l, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, n, p, m, k, j, l, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, p, q, m, k, j, l, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, p, m, n, k, j, l, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, q, p, n, k, j, l, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, m, q, n, k, j, l, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, q, m, p, k, j, l, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, m, n, p, k, j, l, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, n, q, p, k, j, l, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, n, m, q, k, j, l, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, p, n, q, k, j, l, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, m, p, q, k, j, l, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, q, n, m, k, l, i, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, n, p, m, k, l, i, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, p, q, m, k, l, i, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, p, m, n, k, l, i, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, q, p, n, k, l, i, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, m, q, n, k, l, i, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, q, m, p, k, l, i, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, m, n, p, k, l, i, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, n, q, p, k, l, i, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, n, m, q, k, l, i, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, p, n, q, k, l, i, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, m, p, q, k, l, i, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, q, n, m, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, n, p, m, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, p, q, m, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, p, m, n, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, q, p, n, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, m, q, n, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, q, m, p, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, m, n, p, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, n, q, p, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, n, m, q, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, p, n, q, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, m, p, q, l, i, k, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, q, n, m, i, k, l, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, n, p, m, i, k, l, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, p, q, m, i, k, l, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, p, m, n, i, k, l, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, q, p, n, i, k, l, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, m, q, n, i, k, l, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, q, m, p, i, k, l, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, m, n, p, i, k, l, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, n, q, p, i, k, l, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, n, m, q, i, k, l, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, p, n, q, i, k, l, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, m, p, q, i, k, l, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, q, n, m, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, n, p, m, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, p, q, m, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, p, m, n, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, q, p, n, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, m, q, n, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, q, m, p, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, m, n, p, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, n, q, p, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, n, m, q, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, p, n, q, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, m, p, q, l, j, i, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, q, n, m, i, l, j, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, n, p, m, i, l, j, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, p, q, m, i, l, j, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, p, m, n, i, l, j, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, q, p, n, i, l, j, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, m, q, n, i, l, j, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, q, m, p, i, l, j, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, m, n, p, i, l, j, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, n, q, p, i, l, j, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, n, m, q, i, l, j, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, p, n, q, i, l, j, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, m, p, q, i, l, j, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, q, n, m, j, i, l, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, n, p, m, j, i, l, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, p, q, m, j, i, l, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, p, m, n, j, i, l, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, q, p, n, j, i, l, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, m, q, n, j, i, l, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, q, m, p, j, i, l, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, m, n, p, j, i, l, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, n, q, p, j, i, l, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, n, m, q, j, i, l, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, p, n, q, j, i, l, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, m, p, q, j, i, l, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, q, n, m, j, k, i, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, n, p, m, j, k, i, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, p, q, m, j, k, i, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, p, m, n, j, k, i, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, q, p, n, j, k, i, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, m, q, n, j, k, i, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, q, m, p, j, k, i, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, m, n, p, j, k, i, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, n, q, p, j, k, i, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, n, m, q, j, k, i, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, p, n, q, j, k, i, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, m, p, q, j, k, i, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, q, n, m, k, i, j, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, n, p, m, k, i, j, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, p, q, m, k, i, j, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, p, m, n, k, i, j, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, q, p, n, k, i, j, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, m, q, n, k, i, j, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, q, m, p, k, i, j, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, m, n, p, k, i, j, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, n, q, p, k, i, j, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, n, m, q, k, i, j, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, p, n, q, k, i, j, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, m, p, q, k, i, j, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, q, n, m, i, j, k, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, n, p, m, i, j, k, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, p, q, m, i, j, k, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, p, m, n, i, j, k, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, q, p, n, i, j, k, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, m, q, n, i, j, k, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, q, m, p, i, j, k, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, m, n, p, i, j, k, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, n, q, p, i, j, k, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, n, m, q, i, j, k, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, p, n, q, i, j, k, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, m, p, q, i, j, k, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, p, n, m, k, l, j, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, q, p, m, k, l, j, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, n, q, m, k, l, j, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, q, m, n, k, l, j, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, m, p, n, k, l, j, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, p, q, n, k, l, j, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, n, m, p, k, l, j, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, q, n, p, k, l, j, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, m, q, p, k, l, j, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, p, m, q, k, l, j, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, m, n, q, k, l, j, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, n, p, q, k, l, j, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, p, n, m, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, q, p, m, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, n, q, m, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, q, m, n, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, m, p, n, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, p, q, n, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, n, m, p, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, q, n, p, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, m, q, p, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, p, m, q, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, m, n, q, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, n, p, q, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, p, n, m, j, k, l, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, q, p, m, j, k, l, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, n, q, m, j, k, l, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, q, m, n, j, k, l, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, m, p, n, j, k, l, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, p, q, n, j, k, l, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, n, m, p, j, k, l, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, q, n, p, j, k, l, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, m, q, p, j, k, l, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, p, m, q, j, k, l, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, m, n, q, j, k, l, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, n, p, q, j, k, l, i }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, p, n, m, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, q, p, m, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, n, q, m, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, q, m, n, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, m, p, n, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, p, q, n, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, n, m, p, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, q, n, p, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, m, q, p, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, p, m, q, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, m, n, q, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, n, p, q, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, p, n, m, i, l, k, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, q, p, m, i, l, k, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, n, q, m, i, l, k, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, q, m, n, i, l, k, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, m, p, n, i, l, k, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, p, q, n, i, l, k, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, n, m, p, i, l, k, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, q, n, p, i, l, k, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, m, q, p, i, l, k, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, p, m, q, i, l, k, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, m, n, q, i, l, k, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, n, p, q, i, l, k, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, p, n, m, k, i, l, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, q, p, m, k, i, l, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, n, q, m, k, i, l, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, q, m, n, k, i, l, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, m, p, n, k, i, l, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, p, q, n, k, i, l, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, n, m, p, k, i, l, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, q, n, p, k, i, l, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, m, q, p, k, i, l, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, p, m, q, k, i, l, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, m, n, q, k, i, l, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, n, p, q, k, i, l, j }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, p, n, m, j, l, i, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, q, p, m, j, l, i, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, n, q, m, j, l, i, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, q, m, n, j, l, i, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, m, p, n, j, l, i, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, p, q, n, j, l, i, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, n, m, p, j, l, i, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, q, n, p, j, l, i, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, m, q, p, j, l, i, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, p, m, q, j, l, i, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, m, n, q, j, l, i, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, n, p, q, j, l, i, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, p, n, m, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, q, p, m, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, n, q, m, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, q, m, n, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, m, p, n, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, p, q, n, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, n, m, p, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, q, n, p, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, m, q, p, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, p, m, q, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, m, n, q, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, n, p, q, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, p, n, m, i, j, l, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, q, p, m, i, j, l, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, n, q, m, i, j, l, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, q, m, n, i, j, l, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, m, p, n, i, j, l, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, p, q, n, i, j, l, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, n, m, p, i, j, l, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, q, n, p, i, j, l, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, m, q, p, i, j, l, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, p, m, q, i, j, l, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, m, n, q, i, j, l, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, n, p, q, i, j, l, k }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, p, n, m, k, j, i, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, q, p, m, k, j, i, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, n, q, m, k, j, i, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, q, m, n, k, j, i, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, m, p, n, k, j, i, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, p, q, n, k, j, i, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, n, m, p, k, j, i, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, q, n, p, k, j, i, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, m, q, p, k, j, i, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, p, m, q, k, j, i, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, m, n, q, k, j, i, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, n, p, q, k, j, i, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, p, n, m, i, k, j, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, q, p, m, i, k, j, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, n, q, m, i, k, j, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, q, m, n, i, k, j, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, m, p, n, i, k, j, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, p, q, n, i, k, j, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, n, m, p, i, k, j, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, q, n, p, i, k, j, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, m, q, p, i, k, j, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, p, m, q, i, k, j, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, m, n, q, i, k, j, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, n, p, q, i, k, j, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, p, n, m, j, i, k, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, q, p, m, j, i, k, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, n, q, m, j, i, k, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, q, m, n, j, i, k, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, m, p, n, j, i, k, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, p, q, n, j, i, k, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { q, n, m, p, j, i, k, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, q, n, p, j, i, k, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, m, q, p, j, i, k, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { n, p, m, q, j, i, k, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, m, n, q, j, i, k, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { m, n, p, q, j, i, k, l }; spin_batch.push_back( std::make_pair( idx,  -val ) ); 
    idx = { p, q, n, m, k, l, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, n, p, m, k, l, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, p, q, m, k, l, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, p, m, n, k, l, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, q, p, n, k, l, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, m, q, n, k, l, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, q, m, p, k, l, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, m, n, p, k, l, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, n, q, p, k, l, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, n, m, q, k, l, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, p, n, q, k, l, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, m, p, q, k, l, j, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, q, n, m, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, n, p, m, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, p, q, m, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, p, m, n, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, q, p, n, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, m, q, n, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, q, m, p, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, m, n, p, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, n, q, p, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, n, m, q, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, p, n, q, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, m, p, q, l, j, k, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, q, n, m, j, k, l, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, n, p, m, j, k, l, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, p, q, m, j, k, l, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, p, m, n, j, k, l, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, q, p, n, j, k, l, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, m, q, n, j, k, l, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, q, m, p, j, k, l, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, m, n, p, j, k, l, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, n, q, p, j, k, l, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, n, m, q, j, k, l, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, p, n, q, j, k, l, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, m, p, q, j, k, l, i }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, q, n, m, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, n, p, m, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, p, q, m, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, p, m, n, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, q, p, n, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, m, q, n, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, q, m, p, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, m, n, p, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, n, q, p, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, n, m, q, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, p, n, q, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, m, p, q, l, k, i, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, q, n, m, i, l, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, n, p, m, i, l, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, p, q, m, i, l, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, p, m, n, i, l, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, q, p, n, i, l, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, m, q, n, i, l, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, q, m, p, i, l, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, m, n, p, i, l, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, n, q, p, i, l, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, n, m, q, i, l, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, p, n, q, i, l, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, m, p, q, i, l, k, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, q, n, m, k, i, l, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, n, p, m, k, i, l, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, p, q, m, k, i, l, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, p, m, n, k, i, l, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, q, p, n, k, i, l, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, m, q, n, k, i, l, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, q, m, p, k, i, l, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, m, n, p, k, i, l, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, n, q, p, k, i, l, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, n, m, q, k, i, l, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, p, n, q, k, i, l, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, m, p, q, k, i, l, j }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, q, n, m, j, l, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, n, p, m, j, l, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, p, q, m, j, l, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, p, m, n, j, l, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, q, p, n, j, l, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, m, q, n, j, l, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, q, m, p, j, l, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, m, n, p, j, l, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, n, q, p, j, l, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, n, m, q, j, l, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, p, n, q, j, l, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, m, p, q, j, l, i, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, q, n, m, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, n, p, m, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, p, q, m, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, p, m, n, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, q, p, n, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, m, q, n, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, q, m, p, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, m, n, p, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, n, q, p, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, n, m, q, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, p, n, q, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, m, p, q, l, i, j, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, q, n, m, i, j, l, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, n, p, m, i, j, l, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, p, q, m, i, j, l, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, p, m, n, i, j, l, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, q, p, n, i, j, l, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, m, q, n, i, j, l, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, q, m, p, i, j, l, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, m, n, p, i, j, l, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, n, q, p, i, j, l, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, n, m, q, i, j, l, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, p, n, q, i, j, l, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, m, p, q, i, j, l, k }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, q, n, m, k, j, i, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, n, p, m, k, j, i, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, p, q, m, k, j, i, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, p, m, n, k, j, i, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, q, p, n, k, j, i, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, m, q, n, k, j, i, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, q, m, p, k, j, i, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, m, n, p, k, j, i, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, n, q, p, k, j, i, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, n, m, q, k, j, i, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, p, n, q, k, j, i, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, m, p, q, k, j, i, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, q, n, m, i, k, j, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, n, p, m, i, k, j, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, p, q, m, i, k, j, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, p, m, n, i, k, j, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, q, p, n, i, k, j, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, m, q, n, i, k, j, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, q, m, p, i, k, j, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, m, n, p, i, k, j, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, n, q, p, i, k, j, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, n, m, q, i, k, j, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, p, n, q, i, k, j, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, m, p, q, i, k, j, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, q, n, m, j, i, k, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, n, p, m, j, i, k, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, p, q, m, j, i, k, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, p, m, n, j, i, k, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, q, p, n, j, i, k, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, m, q, n, j, i, k, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, q, m, p, j, i, k, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { q, m, n, p, j, i, k, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, n, q, p, j, i, k, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { p, n, m, q, j, i, k, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { m, p, n, q, j, i, k, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
    idx = { n, m, p, q, j, i, k, l }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
  }
}

//===========================================================================================================================================================
//
//void Fourpdm_permutations::get_spatial_permutations( std::vector<std::pair<std::vector<int>,double> >& spin_batch, 
//                                                  const std::vector<int>& indices, const double& val )
//{
//  assert( indices.size() == 8 );
//  std::vector<int> idx;
//  int i = indices[0];
//  int j = indices[1];
//  int k = indices[2];
//  int l = indices[3];
//  int m = indices[4];
//  int n = indices[5];
//  int p = indices[6];
//  int q = indices[7];
//
//  // If indices are not all unique, then all elements should be zero (and next_even_permutation fails)
//  std::vector<int> v = {i,j,k,l};
//  std::sort( v.begin(), v.end() );
//  if ( (v[0]==v[1]) || (v[1]==v[2]) || (v[2]==v[3]) ) return;
//  std::vector<int> w = {m,n,p,q};
//  std::sort( w.begin(), w.end() );
//  if ( (w[0]==w[1]) || (w[1]==w[2]) || (w[2]==w[3]) ) return;
//  bool skip_transpose = ( v == w );
//
//  // The number of possible combinations is (4!)**2 
//  idx = { i, j, k, l, m, n, p, q }; spin_batch.push_back( std::make_pair( idx,  val ) ); 
//
//}
//
//===========================================================================================================================================================

void Pairpdm_permutations::get_spin_permutations( std::vector<std::pair<std::vector<int>,double> >& spin_batch, 
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
    spin_batch.push_back( std::make_pair( idx, -val ) );
    // <m|a_ia_j|n>=-<m|a_ja_i|n>
  }

}
//===========================================================================================================================================================

}
