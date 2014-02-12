/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "nevpt2_4pdm_container.h"

namespace SpinAdapted{

//===========================================================================================================================================================

void Nevpt2_4pdm_container::store_a16_contribution( std::map< std::vector<int>, double >& spatial_batch )
{
  const TwoElectronArray& twoInt = v_2;
  int dim = a16_matrix_.dim1();

  for (auto it = spatial_batch.begin(); it != spatial_batch.end(); ++it) {
    // Note i,j,k... ordering
    int i = (it->first)[0];
    int k = (it->first)[1];
    int m = (it->first)[2];
    int p = (it->first)[3];
    int j = (it->first)[4];
    int l = (it->first)[5];
    int n = (it->first)[6];
    int q = (it->first)[7];
    for ( int a=0; a<dim; ++a )
      a16_matrix_(j,k,i,a,l,q) += twoInt(2*m,2*p,2*n,2*a) * it->second;
    for ( int c=0; c<dim; ++c )
      a16_matrix_(j,k,i,p,l,c) -= twoInt(2*m,2*c,2*n,2*q) * it->second;
    for ( int b=0; b<dim; ++b )
      a16_matrix_(j,k,i,p,b,q) -= twoInt(2*m,2*b,2*n,2*l) * it->second;
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_4pdm_container::build_spatial_elements( std::map< std::vector<int>, double >& spin_batch, 
                                                    std::map< std::vector<int>, double >& spatial_batch )
{

  for (auto it = spatial_batch.begin(); it != spatial_batch.end(); ++it) {
    assert( (it->first).size() == 8 );
    int i = (it->first)[0];
    int j = (it->first)[1];
    int k = (it->first)[2];
    int l = (it->first)[3];
    int m = (it->first)[4];
    int n = (it->first)[5];
    int p = (it->first)[6];
    int q = (it->first)[7];
    // Sum over spin indices
    double val = 0.0;
    for (int s=0; s<2; s++) {
      for (int t=0; t<2; t++) {
        for (int u=0; u<2; u++) {
          for (int v=0; v<2; v++) {
            std::vector<int> idx = { 2*i+s, 2*j+t, 2*k+u, 2*l+v, 2*m+v, 2*n+u, 2*p+t, 2*q+s };
            val += spin_batch[ idx ];
          }
        }
      }
    }
    // Save spatial element
    if ( abs(val) > 1e-14 ) spatial_batch[ it->first ] = val;
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_4pdm_container::store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements)
{
  assert( new_spin_orbital_elements.size() == 70 );
  // Temporary batches of npdm elements
  std::map< std::vector<int>, double > spin_batch;
  std::map< std::vector<int>, double > spatial_batch;

  for (int idx=0; idx < new_spin_orbital_elements.size(); ++idx) {
    // Get all spin-index permutations
    Fourpdm_permutations p;
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
