/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "nevpt2_3pdm_container.h"

namespace SpinAdapted{

//===========================================================================================================================================================

void Nevpt2_3pdm_container::store_a16_contribution( std::map< std::vector<int>, double >& spatial_batch )
{
  const TwoElectronArray& twoInt = v_2;
  int dim = a16_matrix_.dim1();
  int i,j,k,l,m,n,p,q;

  for (auto it = spatial_batch.begin(); it != spatial_batch.end(); ++it) {
    // delta_jk terms
    i = (it->first)[0]; m = (it->first)[1]; p = (it->first)[2]; q = (it->first)[3]; n = (it->first)[4]; l = (it->first)[5];
    for ( int jk=0; jk<dim; ++jk ) {
      for ( int a=0; a<dim; ++a )
        a16_matrix_(jk,jk,i,a,l,q) += twoInt(2*m,2*p,2*n,2*a) * it->second;
      for ( int c=0; c<dim; ++c )
        a16_matrix_(jk,jk,i,p,l,c) -= twoInt(2*m,2*c,2*n,2*q) * it->second;
      for ( int b=0; b<dim; ++b )
        a16_matrix_(jk,jk,i,p,b,q) -= twoInt(2*m,2*b,2*n,2*l) * it->second;
    }
    // delta_jm terms
    i = (it->first)[0]; k = (it->first)[1]; p = (it->first)[2]; q = (it->first)[3]; l = (it->first)[4]; n = (it->first)[5];
    for ( int jm=0; jm<dim; ++jm ) {
      for ( int a=0; a<dim; ++a )
        a16_matrix_(jm,k,i,a,l,q) += twoInt(2*jm,2*p,2*n,2*a) * it->second;
      for ( int c=0; c<dim; ++c )
        a16_matrix_(jm,k,i,p,l,c) -= twoInt(2*jm,2*c,2*n,2*q) * it->second;
      for ( int b=0; b<dim; ++b )
        a16_matrix_(jm,k,i,p,b,q) -= twoInt(2*jm,2*b,2*n,2*l) * it->second;
    }
    // delta_lm terms
    i = (it->first)[0]; k = (it->first)[1]; p = (it->first)[2]; q = (it->first)[3]; n = (it->first)[4]; j = (it->first)[5];
    for ( int lm=0; lm<dim; ++lm ) {
      for ( int a=0; a<dim; ++a )
        a16_matrix_(j,k,i,a,lm,q) += twoInt(2*lm,2*p,2*n,2*a) * it->second;
      for ( int c=0; c<dim; ++c )
        a16_matrix_(j,k,i,p,lm,c) -= twoInt(2*lm,2*c,2*n,2*q) * it->second;
      for ( int b=0; b<dim; ++b )
        a16_matrix_(j,k,i,p,b,q) -= twoInt(2*lm,2*b,2*n,2*lm) * it->second;
    }
    // delta_np terms
    i = (it->first)[0]; k = (it->first)[1]; m = (it->first)[2]; q = (it->first)[3]; l = (it->first)[4]; j = (it->first)[5];
    for ( int np=0; np<dim; ++np ) {
      for ( int a=0; a<dim; ++a )
        a16_matrix_(j,k,i,a,l,q) += twoInt(2*m,2*np,2*np,2*a) * it->second;
      for ( int c=0; c<dim; ++c )
        a16_matrix_(j,k,i,np,l,c) -= twoInt(2*m,2*c,2*np,2*q) * it->second;
      for ( int b=0; b<dim; ++b )
        a16_matrix_(j,k,i,np,b,q) -= twoInt(2*m,2*b,2*np,2*l) * it->second;
    }
    // delta_lp terms
    i = (it->first)[0]; k = (it->first)[1]; m = (it->first)[2]; n = (it->first)[3]; q = (it->first)[4]; j = (it->first)[5];
    for ( int lp=0; lp<dim; ++lp ) {
      for ( int a=0; a<dim; ++a )
        a16_matrix_(j,k,i,a,lp,q) += twoInt(2*m,2*lp,2*n,2*a) * it->second;
      for ( int c=0; c<dim; ++c )
        a16_matrix_(j,k,i,lp,lp,c) -= twoInt(2*m,2*c,2*n,2*q) * it->second;
      for ( int b=0; b<dim; ++b )
        a16_matrix_(j,k,i,lp,b,q) -= twoInt(2*m,2*b,2*n,2*lp) * it->second;
    }
    // delta_jp terms
    i = (it->first)[0]; k = (it->first)[1]; m = (it->first)[2]; n = (it->first)[3]; l = (it->first)[4]; q = (it->first)[5];
    for ( int jp=0; jp<dim; ++jp ) {
      for ( int a=0; a<dim; ++a )
        a16_matrix_(jp,k,i,a,l,q) += twoInt(2*m,2*jp,2*n,2*a) * it->second;
      for ( int c=0; c<dim; ++c )
        a16_matrix_(jp,k,i,jp,l,c) -= twoInt(2*m,2*c,2*n,2*q) * it->second;
      for ( int b=0; b<dim; ++b )
        a16_matrix_(jp,k,i,jp,b,q) -= twoInt(2*m,2*b,2*n,2*l) * it->second;
    }
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_3pdm_container::build_spatial_elements( std::map< std::vector<int>, double >& spin_batch, 
                                                    std::map< std::vector<int>, double >& spatial_batch )
{
  double factor = 1.0;

  for (auto it = spatial_batch.begin(); it != spatial_batch.end(); ++it) {
    assert( (it->first).size() == 6 );
    int i = (it->first)[0];
    int j = (it->first)[1];
    int k = (it->first)[2];
    int l = (it->first)[3];
    int m = (it->first)[4];
    int n = (it->first)[5];
    // Sum over spin indices
    double val = 0.0;
    for (int s=0; s<2; s++) {
      for (int t=0; t<2; t++) {
        for (int u=0; u<2; u++) {
          std::vector<int> idx = { 2*i+s, 2*j+t, 2*k+u, 2*l+u, 2*m+t, 2*n+s };
          val += spin_batch[ idx ];
        }
      }
    }
    // Save spatial element
    if ( abs(val) > 1e-14 ) spatial_batch[ it->first ] = val;
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_3pdm_container::store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements)
{
  assert( new_spin_orbital_elements.size() == 20 );
  // Temporary batches of npdm elements
  std::map< std::vector<int>, double > spin_batch;
  std::map< std::vector<int>, double > spatial_batch;

  for (int idx=0; idx < new_spin_orbital_elements.size(); ++idx) {
    // Get all spin-index permutations
    Threepdm_permutations p;
    std::map< std::vector<int>, int > spin_indices = p.get_spin_permutations( new_spin_orbital_elements[idx].first );
    double val = new_spin_orbital_elements[idx].second;
    for (auto it = spin_indices.begin(); it != spin_indices.end(); ++it) {
      // Initialize spatial indices
      std::vector<int> vec;
      for (int i=0; i < (it->first).size(); ++i)
        vec.push_back( (it->first)[i]/2 );
      spatial_batch[ vec ] = 0.0;
      // Assign temporary batch of spin-orbital elements
      spin_batch[ it->first ] = it->second * val;
    }
  }

  // Build and store new spatial elements
  build_spatial_elements( spin_batch, spatial_batch );
  // Build and store new NEVPT2 A_matrix contributions
  store_a16_contribution( spatial_batch );
}

//===========================================================================================================================================================

}

