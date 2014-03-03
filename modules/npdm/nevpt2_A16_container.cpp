/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "nevpt2_A16_container.h"

namespace SpinAdapted{

//===========================================================================================================================================================

Nevpt2_A16_matrix::Nevpt2_A16_matrix( int sites )
{
  a16_matrix_.resize(sites,sites,sites,sites,sites,sites);
  a16_matrix_.Clear();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_A16_matrix::save_npdms(const int& i, const int& j)
{

//FIXME parallel    
//boost::mpi::communicator world;
//world.barrier();
//    Timer timer;
//    // accumulate_npdm();
//    if ( store_full_spin_array_ ) {
      save_A16_matrix_text();
//      save_npdm_binary(i, j);
//    }
//    if ( store_full_spatial_array_ ) {
//      save_spatial_npdm_text(i, j);
//      save_spatial_npdm_binary(i, j);
//    }
//
//world.barrier();
//    pout << "NEVPT2 A16_matrix save full array time " << timer.elapsedwalltime() << " " << timer.elapsedcputime() << endl;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_A16_matrix::save_A16_matrix_text()
{

 if( mpigetrank() == 0)
  {
    int dim = a16_matrix_.dim1();
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/A16_matrix.", 0, 0,".txt");
    ofstream ofs(file);
    ofs << dim << endl;

    double trace = 0.0;
    for(int i=0; i<dim; ++i)
      for(int j=0; j<dim; ++j)
        for(int k=0; k<dim; ++k)
          for(int l=0; l<dim; ++l)
            for(int m=0; m<dim; ++m)
              for(int n=0; n<dim; ++n) {
                if ( abs(a16_matrix_(i,j,k,l,m,n)) > 1e-14 ) {
                  ofs << boost::format("%d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % a16_matrix_(i,j,k,l,m,n);
                }
              }
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_A16_matrix::store_A16_4pdm_contribution( std::map< std::vector<int>, double >& spatial_batch )
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

void Nevpt2_A16_matrix::store_A16_3pdm_contribution( std::map< std::vector<int>, double >& spatial_batch )
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

void Nevpt2_A16_matrix::store_A16_2pdm_contribution( std::map< std::vector<int>, double >& spatial_batch )
{
  const TwoElectronArray& twoInt = v_2;
  int dim = a16_matrix_.dim1();
  int i,j,k,l,m,n,p,q;

assert(false);
  for (auto it = spatial_batch.begin(); it != spatial_batch.end(); ++it) {
    // delta_jk * delta_lm terms
    i = (it->first)[0]; p = (it->first)[1]; q = (it->first)[2]; n = (it->first)[3];
    for ( int jk=0; jk<dim; ++jk ) {
      for ( int lm=0; lm<dim; ++lm ) {
        for ( int a=0; a<dim; ++a )
          a16_matrix_(jk,jk,i,a,lm,q) += twoInt(2*lm,2*p,2*n,2*a) * it->second;
        for ( int c=0; c<dim; ++c )
          a16_matrix_(jk,jk,i,p,lm,c) -= twoInt(2*lm,2*c,2*n,2*q) * it->second;
        for ( int b=0; b<dim; ++b )
          a16_matrix_(jk,jk,i,p,b,q) -= twoInt(2*lm,2*b,2*n,2*lm) * it->second;
      }
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_A16_matrix::build_spatial_4pdm_elements( std::map< std::vector<int>, double >& spin_batch, 
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

void Nevpt2_A16_matrix::build_spatial_3pdm_elements( std::map< std::vector<int>, double >& spin_batch, 
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

void Nevpt2_A16_matrix::build_spatial_2pdm_elements( std::map< std::vector<int>, double >& spin_batch, 
                                                     std::map< std::vector<int>, double >& spatial_batch )
{
  // Note we multiply the spatial 2PDM by a factor of 1/2 to be consistent with the old BLOCK code
  double factor = 1.0;

  for (auto it = spatial_batch.begin(); it != spatial_batch.end(); ++it) {
    int i = (it->first)[0];
    int j = (it->first)[1];
    int k = (it->first)[2];
    int l = (it->first)[3];
    // Sum over spin indices 
    double val = 0.0;
    for (int s=0; s<2; s++) {
      for (int t =0; t<2; t++) {
        std::vector<int> idx = { 2*i+s, 2*j+t, 2*k+t, 2*l+s };
        val += spin_batch[ idx ];
      }
    }
    // Save spatial element
    if ( abs(val) > 1e-14 ) spatial_batch[ it->first ] = factor * val;
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_A16_matrix::build_npdm_batch( Npdm_permutations& p, const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements,
                                          std::map< std::vector<int>, double >& spatial_batch )
{
//  std::map< std::vector<int>, double > spin_batch;
//
//  for (int idx=0; idx < new_spin_orbital_elements.size(); ++idx) {
//    // Get all spin-index permutations
//    std::map< std::vector<int>, int > spin_indices = p.get_spin_permutations( new_spin_orbital_elements[idx].first );
//    double val = new_spin_orbital_elements[idx].second;
//    for (auto it = spin_indices.begin(); it != spin_indices.end(); ++it) {
//      // Initialize spatial indices
//      std::vector<int> vec;
//      for (int i=0; i < (it->first).size(); ++i)
//        vec.push_back( (it->first)[i]/2 );
//      spatial_batch[ vec ] = 0.0;
//      // Assign temporary batch of spin-orbital elements
//      spin_batch[ it->first ] = it->second * val;
//    }
//  }
//
//  // Build spatial elements
//  if ( new_spin_orbital_elements[0].first.size() == 4 ) {
//    build_spatial_2pdm_elements( spin_batch, spatial_batch );
//  }
//  else if ( new_spin_orbital_elements[0].first.size() == 6 ) {
//    build_spatial_3pdm_elements( spin_batch, spatial_batch );
//  }
//  else if ( new_spin_orbital_elements[0].first.size() == 8 ) {
//    build_spatial_4pdm_elements( spin_batch, spatial_batch );
//  }
//  else assert(false);

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_A16_matrix::store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements)
{
//  // Temporary batch of spatial npdm elements
//  std::map< std::vector<int>, double > spatial_batch;
//  
//  // 2PDM elements
//  if ( new_spin_orbital_elements[0].first.size() == 4 ) {
//    assert( new_spin_orbital_elements.size() == 6 );
//    Twopdm_permutations p;
//    build_npdm_batch( p, new_spin_orbital_elements, spatial_batch );
//    store_A16_2pdm_contribution( spatial_batch );
//  }
//  // 3PDM elements
//  else if ( new_spin_orbital_elements[0].first.size() == 6 ) {
//    assert( new_spin_orbital_elements.size() == 20 );
//    Threepdm_permutations p;
//    build_npdm_batch( p, new_spin_orbital_elements, spatial_batch );
//    store_A16_3pdm_contribution( spatial_batch );
//  }
//  // 4PDM elements
//  else if ( new_spin_orbital_elements[0].first.size() == 8 ) {
//    assert( new_spin_orbital_elements.size() == 70 );
//    Fourpdm_permutations p;
//    build_npdm_batch( p, new_spin_orbital_elements, spatial_batch );
//    store_A16_4pdm_contribution( spatial_batch );
//  }
//  else assert(false);
}

//===========================================================================================================================================================

}

