/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "global.h"
#include "nevpt2_A16_container.h"
#include <boost/format.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/map.hpp>
#include "npdm_array_buffer.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

Nevpt2_A16_matrix::Nevpt2_A16_matrix( int sites )
#ifndef DEBUG_A16_FULL_MATRIX
 : a16_matrix_( Npdm_array_buffer(6,sites) )
#endif
{
#ifdef DEBUG_A16_FULL_MATRIX
  a16_matrix_.resize(sites,sites,sites,sites,sites,sites);
  a16_matrix_.Clear();
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_A16_matrix::save_npdms(const int& i, const int& j)
{
#ifdef DEBUG_A16_FULL_MATRIX
  accumulate_full_array();
  save_full_array_text( a16_matrix_ );
#else
  a16_matrix_.close_buffer();
  accumulate_files();
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

#ifdef DEBUG_A16_FULL_MATRIX
void Nevpt2_A16_matrix::accumulate_full_array()
{
#ifndef SERIAL
  int dim = a16_matrix_.dim1();
  array_6d<double> tmp_recv;
  mpi::communicator world;
  if( mpigetrank() == 0)
  {
    for(int p=1; p<world.size(); ++p) {
      world.recv(p, p, tmp_recv);
      for(int i=0; i<dim; ++i)
        for(int j=0; j<dim; ++j)
          for(int k=0; k<dim; ++k)
            for(int l=0; l<dim; ++l)
              for(int m=0; m<dim; ++m)
                for(int n=0; n<dim; ++n)
                  // This only works because we carefully ensured no duplicate npdm elements are built on different processors
                  if( abs(tmp_recv(i,j,k,l,m,n)) > NUMERICAL_ZERO )
                    a16_matrix_(i,j,k,l,m,n) += tmp_recv(i,j,k,l,m,n);
    }
  }
  else
  {
    world.send(0, mpigetrank(), a16_matrix_);
  }
#endif
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

#ifndef DEBUG_A16_FULL_MATRIX
void Nevpt2_A16_matrix::accumulate_files()
{
  if( mpigetrank() == 0) {

    // Initialize full array to accumulate files into
    int dim = a16_matrix_.dim1();
    array_6d<double> mat;
    mat.resize(dim,dim,dim,dim,dim,dim);
    mat.Clear();

    // Get file names
    namespace fs = boost::filesystem ;
    std::vector<std::string> names ;
    std::string dir = dmrginp.save_prefix();
    assert( fs::exists(dir) );
    fs::directory_iterator it(dir) ;
    fs::directory_iterator end ;
    while ( it != end ) {
      names.push_back(it->path().filename().string()) ;
      ++it ;
    }

    // Filter out required files
    std::string pattern = "partial_A16_matrix.0.0";
    auto pos = std::remove_if( std::begin(names), std::end(names), [&](std::string& s) { return s.find(pattern) == std::string::npos ; } ) ; 
    names.erase(pos, std::end(names)) ;

    // Read in and accumulate data from all files
    for ( auto& name : names ) {
      std::map< std::vector<int>, double > buff;
      std::ifstream ifs(name.c_str(), std::ios::binary);
      boost::archive::binary_iarchive load(ifs);
      load >> buff;
      ifs.close();
      for ( auto it = buff.begin(); it != buff.end(); ++it ) {
        int i = (it->first)[0]; int j = (it->first)[1]; int k = (it->first)[2]; 
        int l = (it->first)[3]; int m = (it->first)[4]; int n = (it->first)[5];
        mat(i,j,k,l,m,n) += it->second;
      }
    }
    save_full_array_text( mat );
  }

}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_A16_matrix::save_full_array_text( array_6d<double>& matrix )
{
 if( mpigetrank() == 0)
  {
    int dim = matrix.dim1();
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
                if ( abs(matrix(i,j,k,l,m,n)) > NUMERICAL_ZERO ) {
                  ofs << boost::format("%d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % matrix(i,j,k,l,m,n);
                }
              }
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_A16_matrix::update_4pdm_contribution( std::vector< std::pair< std::vector<int>, double > >& spin_batch )
{
  const TwoElectronArray& twoInt = v_2;
  int dim = a16_matrix_.dim1();
  int i,j,k,l,m,n,p,q;

  // Take into account orbital reordering
  const std::vector<int>& ro = dmrginp.reorder_vector();

  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    // Store significant elements only
    if ( abs(it->second) < NUMERICAL_ZERO ) continue;
    // Spin indices
    int i = (it->first)[0]; int k = (it->first)[1]; int m = (it->first)[2]; int p = (it->first)[3];
    int j = (it->first)[4]; int l = (it->first)[5]; int n = (it->first)[6]; int q = (it->first)[7];

    if ( i%2 != q%2 ) continue;
    if ( k%2 != n%2 ) continue;
    if ( m%2 != l%2 ) continue;
    if ( p%2 != j%2 ) continue;

    // Space indices
    i = ro.at(i/2); k = ro.at(k/2); m = ro.at(m/2); p = ro.at(p/2);
    j = ro.at(j/2); l = ro.at(l/2); n = ro.at(n/2); q = ro.at(q/2);

    for ( int a=0; a<dim; ++a )
      a16_matrix_(j,k,i,a,l,q) += twoInt(2*m,2*p,2*n,2*a) * it->second;
    for ( int c=0; c<dim; ++c )
      a16_matrix_(j,k,i,p,l,c) -= twoInt(2*m,2*c,2*n,2*q) * it->second;
    for ( int b=0; b<dim; ++b )
      a16_matrix_(j,k,i,p,b,q) -= twoInt(2*m,2*b,2*n,2*l) * it->second;
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_A16_matrix::update_3pdm_contribution( std::vector< std::pair< std::vector<int>, double > >& spin_batch )
{
  const TwoElectronArray& twoInt = v_2;
  int dim = a16_matrix_.dim1();
  int i,j,k,l,m,n,p,q;

  // Take into account orbital reordering
  const std::vector<int>& ro = dmrginp.reorder_vector();

  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    // Store significant elements only
    if ( abs(it->second) < NUMERICAL_ZERO ) continue;

    i = (it->first)[0]; j = (it->first)[1]; k = (it->first)[2]; l = (it->first)[3]; m = (it->first)[4]; n = (it->first)[5];
    if ( i%2 != n%2 ) continue;
    if ( j%2 != m%2 ) continue;
    if ( k%2 != l%2 ) continue;

    // delta_jk terms
    //--------------------
    // Spin indices
    i = (it->first)[0]; m = (it->first)[1]; p = (it->first)[2]; q = (it->first)[3]; n = (it->first)[4]; l = (it->first)[5];
    // Space indices
    i = ro.at(i/2); m = ro.at(m/2); p = ro.at(p/2); q = ro.at(q/2); n = ro.at(n/2); l = ro.at(l/2); 
    for ( int jk=0; jk<dim; ++jk ) {
      for ( int a=0; a<dim; ++a )
        a16_matrix_(jk,jk,i,a,l,q) += twoInt(2*m,2*p,2*n,2*a) * it->second;
      for ( int c=0; c<dim; ++c )
        a16_matrix_(jk,jk,i,p,l,c) -= twoInt(2*m,2*c,2*n,2*q) * it->second;
      for ( int b=0; b<dim; ++b )
        a16_matrix_(jk,jk,i,p,b,q) -= twoInt(2*m,2*b,2*n,2*l) * it->second;
    }
    // delta_jm terms
    //--------------------
    // Spin indices
    i = (it->first)[0]; k = (it->first)[1]; p = (it->first)[2]; q = (it->first)[3]; l = (it->first)[4]; n = (it->first)[5];
    // Space indices
    i = ro.at(i/2); k = ro.at(k/2); p = ro.at(p/2); q = ro.at(q/2); l = ro.at(l/2); n = ro.at(n/2); 
    for ( int jm=0; jm<dim; ++jm ) {
      for ( int a=0; a<dim; ++a )
        a16_matrix_(jm,k,i,a,l,q) += twoInt(2*jm,2*p,2*n,2*a) * it->second;
      for ( int c=0; c<dim; ++c )
        a16_matrix_(jm,k,i,p,l,c) -= twoInt(2*jm,2*c,2*n,2*q) * it->second;
      for ( int b=0; b<dim; ++b )
        a16_matrix_(jm,k,i,p,b,q) -= twoInt(2*jm,2*b,2*n,2*l) * it->second;
    }
    // delta_lm terms
    //--------------------
    // Spin indices
    i = (it->first)[0]; k = (it->first)[1]; p = (it->first)[2]; q = (it->first)[3]; n = (it->first)[4]; j = (it->first)[5];
    // Space indices
    i = ro.at(i/2); k = ro.at(k/2); p = ro.at(p/2); q = ro.at(q/2); n = ro.at(n/2); j = ro.at(j/2); 
    for ( int lm=0; lm<dim; ++lm ) {
      for ( int a=0; a<dim; ++a )
        a16_matrix_(j,k,i,a,lm,q) += twoInt(2*lm,2*p,2*n,2*a) * it->second;
      for ( int c=0; c<dim; ++c )
        a16_matrix_(j,k,i,p,lm,c) -= twoInt(2*lm,2*c,2*n,2*q) * it->second;
      for ( int b=0; b<dim; ++b )
        a16_matrix_(j,k,i,p,b,q) -= twoInt(2*lm,2*b,2*n,2*lm) * it->second;
    }
    // delta_np terms
    //--------------------
    // Spin indices
    i = (it->first)[0]; k = (it->first)[1]; m = (it->first)[2]; q = (it->first)[3]; l = (it->first)[4]; j = (it->first)[5];
    // Space indices
    i = ro.at(i/2); k = ro.at(k/2); m = ro.at(m/2); q = ro.at(q/2); l = ro.at(l/2); j = ro.at(j/2); 
    for ( int np=0; np<dim; ++np ) {
      for ( int a=0; a<dim; ++a )
        a16_matrix_(j,k,i,a,l,q) += twoInt(2*m,2*np,2*np,2*a) * it->second;
      for ( int c=0; c<dim; ++c )
        a16_matrix_(j,k,i,np,l,c) -= twoInt(2*m,2*c,2*np,2*q) * it->second;
      for ( int b=0; b<dim; ++b )
        a16_matrix_(j,k,i,np,b,q) -= twoInt(2*m,2*b,2*np,2*l) * it->second;
    }
    // delta_lp terms
    //--------------------
    // Spin indices
    i = (it->first)[0]; k = (it->first)[1]; m = (it->first)[2]; n = (it->first)[3]; q = (it->first)[4]; j = (it->first)[5];
    // Space indices
    i = ro.at(i/2); k = ro.at(k/2); m = ro.at(m/2); n = ro.at(n/2); q = ro.at(q/2); j = ro.at(j/2); 
    for ( int lp=0; lp<dim; ++lp ) {
      for ( int a=0; a<dim; ++a )
        a16_matrix_(j,k,i,a,lp,q) += twoInt(2*m,2*lp,2*n,2*a) * it->second;
      for ( int c=0; c<dim; ++c )
        a16_matrix_(j,k,i,lp,lp,c) -= twoInt(2*m,2*c,2*n,2*q) * it->second;
      for ( int b=0; b<dim; ++b )
        a16_matrix_(j,k,i,lp,b,q) -= twoInt(2*m,2*b,2*n,2*lp) * it->second;
    }
    // delta_jp terms
    //--------------------
    // Spin indices
    i = (it->first)[0]; k = (it->first)[1]; m = (it->first)[2]; n = (it->first)[3]; l = (it->first)[4]; q = (it->first)[5];
    // Space indices
    i = ro.at(i/2); k = ro.at(k/2); m = ro.at(m/2); n = ro.at(n/2); l = ro.at(l/2); q = ro.at(q/2); 
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

void Nevpt2_A16_matrix::update_2pdm_contribution( std::vector< std::pair< std::vector<int>, double > >& spin_batch )
{
  const TwoElectronArray& twoInt = v_2;
  int dim = a16_matrix_.dim1();
  int i,j,k,l,m,n,p,q;

  // Take into account orbital reordering
  const std::vector<int>& ro = dmrginp.reorder_vector();

  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    // Store significant elements only
    if ( abs(it->second) < NUMERICAL_ZERO ) continue;

    i = (it->first)[0]; j = (it->first)[1]; k = (it->first)[2]; l = (it->first)[3];
    if ( i%2 != l%2 ) continue;
    if ( j%2 != k%2 ) continue;

    // delta_jk * delta_lm
    //--------------------
    // Spin indices
    i = (it->first)[0]; p = (it->first)[1]; q = (it->first)[2]; n = (it->first)[3];
    // Space indices
    i = ro.at(i/2); p = ro.at(p/2); q = ro.at(q/2); n = ro.at(n/2); 
    for ( int jk=0; jk<dim; ++jk ) {
      for ( int lm=0; lm<dim; ++lm ) {
        // Three <EEEE> terms in eqn.(A16)
        for ( int a=0; a<dim; ++a )
          a16_matrix_(jk,jk,i,a,lm,q) += twoInt(2*lm,2*p,2*n,2*a) * it->second;
        for ( int c=0; c<dim; ++c )
          a16_matrix_(jk,jk,i,p,lm,c) -= twoInt(2*lm,2*c,2*n,2*q) * it->second;
        for ( int b=0; b<dim; ++b )
          a16_matrix_(jk,jk,i,p,b,q) -= twoInt(2*lm,2*b,2*n,2*lm) * it->second;
      }
    }
    // delta_jk * delta_np
    //--------------------
    // Spin indices
    i = (it->first)[0]; m = (it->first)[1]; q = (it->first)[2]; l = (it->first)[3];
    // Space indices
    i = ro.at(i/2); m = ro.at(m/2); q = ro.at(q/2); l = ro.at(l/2); 
    for ( int jk=0; jk<dim; ++jk ) {
      for ( int np=0; np<dim; ++np ) {
        // Three <EEEE> terms in eqn.(A16)
        for ( int a=0; a<dim; ++a )
          a16_matrix_(jk,jk,i,a,l,q) += twoInt(2*m,2*np,2*np,2*a) * it->second;
        for ( int c=0; c<dim; ++c )
          a16_matrix_(jk,jk,i,np,l,c) -= twoInt(2*m,2*c,2*np,2*q) * it->second;
        for ( int b=0; b<dim; ++b )
          a16_matrix_(jk,jk,i,np,b,q) -= twoInt(2*m,2*b,2*np,2*l) * it->second;
      }
    }
    // delta_jk * delta_lp
    //--------------------
    // Spin indices
    i = (it->first)[0]; m = (it->first)[1]; n = (it->first)[2]; q = (it->first)[3];
    // Space indices
    i = ro.at(i/2); m = ro.at(m/2); n = ro.at(n/2); q = ro.at(q/2); 
    for ( int jk=0; jk<dim; ++jk ) {
      for ( int lp=0; lp<dim; ++lp ) {
        // Three <EEEE> terms in eqn.(A16)
        for ( int a=0; a<dim; ++a )
          a16_matrix_(jk,jk,i,a,lp,q) += twoInt(2*m,2*lp,2*n,2*a) * it->second;
        for ( int c=0; c<dim; ++c )
          a16_matrix_(jk,jk,i,lp,lp,c) -= twoInt(2*m,2*c,2*n,2*q) * it->second;
        for ( int b=0; b<dim; ++b )
          a16_matrix_(jk,jk,i,lp,b,q) -= twoInt(2*m,2*b,2*n,2*lp) * it->second;
      }
    }
    // delta_jm * delta_np
    //--------------------
    // Spin indices
    i = (it->first)[0]; k = (it->first)[1]; l = (it->first)[2]; q = (it->first)[3];
    // Space indices
    i = ro.at(i/2); k = ro.at(k/2); l = ro.at(l/2); q = ro.at(q/2); 
    for ( int jm=0; jm<dim; ++jm ) {
      for ( int np=0; np<dim; ++np ) {
        // Three <EEEE> terms in eqn.(A16)
        for ( int a=0; a<dim; ++a )
          a16_matrix_(jm,k,i,a,l,q) += twoInt(2*jm,2*np,2*np,2*a) * it->second;
        for ( int c=0; c<dim; ++c )
          a16_matrix_(jm,k,i,np,l,c) -= twoInt(2*jm,2*c,2*np,2*q) * it->second;
        for ( int b=0; b<dim; ++b )
          a16_matrix_(jm,k,i,np,b,q) -= twoInt(2*jm,2*b,2*np,2*l) * it->second;
      }
    }
    // delta_jm * delta_lp
    //--------------------
    // Spin indices
    i = (it->first)[0]; k = (it->first)[1]; q = (it->first)[2]; n = (it->first)[3];
    // Space indices
    i = ro.at(i/2); k = ro.at(k/2); q = ro.at(q/2); n = ro.at(n/2); 
    for ( int jm=0; jm<dim; ++jm ) {
      for ( int lp=0; lp<dim; ++lp ) {
        // Three <EEEE> terms in eqn.(A16)
        for ( int a=0; a<dim; ++a )
          a16_matrix_(jm,k,i,a,lp,q) += twoInt(2*jm,2*lp,2*n,2*a) * it->second;
        for ( int c=0; c<dim; ++c )
          a16_matrix_(jm,k,i,lp,lp,c) -= twoInt(2*jm,2*c,2*n,2*q) * it->second;
        for ( int b=0; b<dim; ++b )
          a16_matrix_(jm,k,i,lp,b,q) -= twoInt(2*jm,2*b,2*n,2*lp) * it->second;
      }
    }
    // delta_lm * delta_np
    //--------------------
    // Spin indices
    i = (it->first)[0]; k = (it->first)[1]; q = (it->first)[2]; j = (it->first)[3];
    // Space indices
    i = ro.at(i/2); k = ro.at(k/2); q = ro.at(q/2); j = ro.at(j/2); 
    for ( int lm=0; lm<dim; ++lm ) {
      for ( int np=0; np<dim; ++np ) {
        // Three <EEEE> terms in eqn.(A16)
        for ( int a=0; a<dim; ++a )
          a16_matrix_(j,k,i,a,lm,q) += twoInt(2*lm,2*np,2*np,2*a) * it->second;
        for ( int c=0; c<dim; ++c )
          a16_matrix_(j,k,i,np,lm,c) -= twoInt(2*lm,2*c,2*np,2*q) * it->second;
        for ( int b=0; b<dim; ++b )
          a16_matrix_(j,k,i,np,b,q) -= twoInt(2*lm,2*b,2*np,2*lm) * it->second;
      }
    }
    // delta_lm * delta_jp
    //--------------------
    // Spin indices
    i = (it->first)[0]; k = (it->first)[1]; n = (it->first)[2]; q = (it->first)[3];
    // Space indices
    i = ro.at(i/2); k = ro.at(k/2); n = ro.at(n/2); q = ro.at(q/2); 
    for ( int lm=0; lm<dim; ++lm ) {
      for ( int jp=0; jp<dim; ++jp ) {
        // Three <EEEE> terms in eqn.(A16)
        for ( int a=0; a<dim; ++a )
          a16_matrix_(jp,k,i,a,lm,q) += twoInt(2*lm,2*jp,2*n,2*a) * it->second;
        for ( int c=0; c<dim; ++c )
          a16_matrix_(jp,k,i,jp,lm,c) -= twoInt(2*lm,2*c,2*n,2*q) * it->second;
        for ( int b=0; b<dim; ++b )
          a16_matrix_(jp,k,i,jp,b,q) -= twoInt(2*lm,2*b,2*n,2*lm) * it->second;
      }
    }

  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_A16_matrix::update_1pdm_contribution( std::vector< std::pair< std::vector<int>, double > >& spin_batch )
{

  const TwoElectronArray& twoInt = v_2;
  int dim = a16_matrix_.dim1();

  // Take into account orbital reordering
  const std::vector<int>& ro = dmrginp.reorder_vector();

  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    // Store significant elements only
    if ( abs(it->second) < NUMERICAL_ZERO ) continue;
    // Spin indices
    int i = (it->first)[0];
    int q = (it->first)[1];
    if ( i%2 != q%2 ) continue;
    // Space indices
    i = ro.at(i/2); 
    q = ro.at(q/2); 
    // delta_jk * delta_lm * delta_np
    for ( int jk=0; jk<dim; ++jk ) {
      for ( int lm=0; lm<dim; ++lm ) {
        for ( int np=0; np<dim; ++np ) {
          // Three <EEEE> terms in eqn.(A16)
          for ( int a=0; a<dim; ++a )
            a16_matrix_(jk,jk,i,a,lm,q) += twoInt(2*lm,2*np,2*np,2*a) * it->second;
          for ( int c=0; c<dim; ++c )
            a16_matrix_(jk,jk,i,np,lm,c) -= twoInt(2*lm,2*c,2*np,2*q) * it->second;
          for ( int b=0; b<dim; ++b )
            a16_matrix_(jk,jk,i,np,b,q) -= twoInt(2*lm,2*b,2*np,2*lm) * it->second;
        }
      }
    }
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_A16_matrix::store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements)
{

  std::vector< std::pair< std::vector<int>, double > > dummy;
  std::vector< std::pair< std::vector<int>, double > > spin_batch;

  if ( new_spin_orbital_elements[0].first.size() == 2 ) {
    // 1PDM elements
    Onepdm_permutations perm;
    perm.process_new_elements( new_spin_orbital_elements, dummy, spin_batch );
    update_1pdm_contribution( spin_batch );
  }
  else if ( new_spin_orbital_elements[0].first.size() == 4 ) {
    // 2PDM elements
    Twopdm_permutations perm;
    perm.process_new_elements( new_spin_orbital_elements, dummy, spin_batch );
    update_2pdm_contribution( spin_batch );
  }
  else if ( new_spin_orbital_elements[0].first.size() == 6 ) {
    // 3PDM elements
    Threepdm_permutations perm;
    perm.process_new_elements( new_spin_orbital_elements, dummy, spin_batch );
    update_3pdm_contribution( spin_batch );
  }
  else if ( new_spin_orbital_elements[0].first.size() == 8 ) {
    // 4PDM elements
    Fourpdm_permutations perm;
    perm.process_new_elements( new_spin_orbital_elements, dummy, spin_batch );
    update_4pdm_contribution( spin_batch );
  }
  else abort();

}

//===========================================================================================================================================================

}
}

