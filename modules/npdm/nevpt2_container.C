/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "global.h"
#include "nevpt2_container.h"
#include <boost/format.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/map.hpp>

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

Nevpt2_container::Nevpt2_container( int sites )
{
//  if ( dmrginp.spatpdm_disk_dump() ){
    a16_matrix_.resize(sites,sites,sites,sites,sites,sites);
    a16_matrix_.Clear();
    a22_matrix_.resize(sites,sites,sites,sites,sites,sites);
    a22_matrix_.Clear();
//  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_container::save_npdms(const int& i, const int& j)
{
//  if ( dmrginp.spatpdm_disk_dump() ){
  char file[5000];
  accumulate_full_array(a16_matrix_ );
  sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/A16_matrix.", i, j,".txt");
  save_full_array_text( a16_matrix_ ,file);

//
  accumulate_full_array(a22_matrix_ );
  sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/A22_matrix.", i, j,".txt");
  save_full_array_text( a22_matrix_ ,file);


//  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_container::save_full_array_text( array_6d<double>& matrix, char file[])
{
 if( mpigetrank() == 0)
  {
    const std::vector<int>& ro = dmrginp.reorder_vector();
    int dim = matrix.dim1();
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
                  ofs << boost::format("%d %d %d %d %d %d %20.14e\n") % ro[i] % ro[j] % ro[k] % ro[l] % ro[m] % ro[n] % matrix(i,j,k,l,m,n);
                }
            }
    ofs.close();
  }
}

void Nevpt2_container::save_full_array_text( array_4d<double>& matrix, char file[])
{
 if( mpigetrank() == 0)
  {
    const std::vector<int>& ro = dmrginp.reorder_vector();
    int dim = matrix.dim1();
    ofstream ofs(file);
    ofs << dim << endl;

    double trace = 0.0;
    for(int i=0; i<dim; ++i)
      for(int j=0; j<dim; ++j)
        for(int k=0; k<dim; ++k)
          for(int l=0; l<dim; ++l)
            if ( abs(matrix(i,j,k,l)) > NUMERICAL_ZERO ) {
              ofs << boost::format("%d %d %d %d %20.14e\n") % ro[i] % ro[j] % ro[k] % ro[l] %  matrix(i,j,k,l);
            }
    ofs.close();
  }
}

void Nevpt2_container::save_full_array_text( array_2d<double>& matrix, char file[])
{
 if( mpigetrank() == 0)
  {
    const std::vector<int>& ro = dmrginp.reorder_vector();
    int dim = matrix.dim1();
    ofstream ofs(file);
    ofs << dim << endl;

    double trace = 0.0;
    for(int i=0; i<dim; ++i)
      for(int j=0; j<dim; ++j)
        if ( abs(matrix(i,j)) > NUMERICAL_ZERO ) {
          ofs << boost::format("%d %d %20.14e\n") % ro[i] % ro[j] %  matrix(i,j);
        }
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//
void Nevpt2_container::save_full_array_text( std::vector<double>& matrix, char file[])
{
 if( mpigetrank() == 0)
  {
    ofstream ofs(file);
    ofs << matrix.size() << endl;

    double trace = 0.0;
    for(int i=0; i<matrix.size(); ++i)
        if ( abs(matrix[i]) > NUMERICAL_ZERO ) {
          ofs << boost::format("%d %20.14e\n") % i %  matrix[i];
        }
    ofs.close();
  }
}

void Nevpt2_container::accumulate_full_array(array_6d<double>& matrix)
{
#ifndef SERIAL
  mpi::communicator world;
#ifdef NDEBUG
  if( mpigetrank() == 0)
    MPI_Reduce(MPI_IN_PLACE, matrix.data(),  matrix.size(), MPI_DOUBLE, MPI_SUM, 0, world);
  else
    MPI_Reduce(matrix.data(), matrix.data(),  matrix.size(), MPI_DOUBLE, MPI_SUM, 0, world);
#else
  int dim = matrix.dim1();
  array_6d<double> tmp_recv;
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
                    matrix(i,j,k,l,m,n) += tmp_recv(i,j,k,l,m,n);
    }
  }
  else
  {
    world.send(0, mpigetrank(), matrix);
  }
#endif
#endif
}

void Nevpt2_container::accumulate_full_array(array_4d<double>& matrix)
{
#ifndef SERIAL
  mpi::communicator world;
#ifdef NDEBUG
  if( mpigetrank() == 0)
    MPI_Reduce(MPI_IN_PLACE, matrix.data(),  matrix.size(), MPI_DOUBLE, MPI_SUM, 0, world);
  else
    MPI_Reduce(matrix.data(), matrix.data(),  matrix.size(), MPI_DOUBLE, MPI_SUM, 0, world);
#else
  int dim = matrix.dim1();
  array_4d<double> tmp_recv;
  if( mpigetrank() == 0)
  {
    for(int p=1; p<world.size(); ++p) {
      world.recv(p, p, tmp_recv);
      for(int i=0; i<dim; ++i)
        for(int j=0; j<dim; ++j)
          for(int k=0; k<dim; ++k)
            for(int l=0; l<dim; ++l)
              // This only works because we carefully ensured no duplicate npdm elements are built on different processors
              if( abs(tmp_recv(i,j,k,l)) > NUMERICAL_ZERO )
                  matrix(i,j,k,l) += tmp_recv(i,j,k,l);
    }
  }
  else
  {
    world.send(0, mpigetrank(), matrix);
  }
#endif
#endif
}

void Nevpt2_container::accumulate_full_array(array_2d<double>& matrix)
{
#ifndef SERIAL
  mpi::communicator world;
#ifdef NDEBUG
  if( mpigetrank() == 0)
    MPI_Reduce(MPI_IN_PLACE, matrix.data(),  matrix.size(), MPI_DOUBLE, MPI_SUM, 0, world);
  else
    MPI_Reduce(matrix.data(), matrix.data(),  matrix.size(), MPI_DOUBLE, MPI_SUM, 0, world);
#else
  int dim = matrix.dim1();
  array_2d<double> tmp_recv;
  if( mpigetrank() == 0)
  {
    for(int p=1; p<world.size(); ++p) {
      world.recv(p, p, tmp_recv);
      for(int i=0; i<dim; ++i)
        for(int j=0; j<dim; ++j)
              // This only works because we carefully ensured no duplicate npdm elements are built on different processors
              if( abs(tmp_recv(i,j)) > NUMERICAL_ZERO )
                  matrix(i,j) += tmp_recv(i,j);
    }
  }
  else
  {
    world.send(0, mpigetrank(), matrix);
  }
#endif
#endif
}

void Nevpt2_container::accumulate_full_array(std::vector<double>& matrix)
{
#ifndef SERIAL
  mpi::communicator world;
#ifdef NDEBUG
  if( mpigetrank() == 0)
    MPI_Reduce(MPI_IN_PLACE, matrix.data(),  matrix.size(), MPI_DOUBLE, MPI_SUM, 0, world);
  else
    MPI_Reduce(matrix.data(), matrix.data(),  matrix.size(), MPI_DOUBLE, MPI_SUM, 0, world);
#else
  int dim = matrix.size();
  std::vector<double> tmp_recv;
  if( mpigetrank() == 0)
  {
    for(int p=1; p<world.size(); ++p) {
      world.recv(p, p, tmp_recv);
      for(int i=0; i<dim; ++i)
        // This only works because we carefully ensured no duplicate npdm elements are built on different processors
        if( abs(tmp_recv[i]) > NUMERICAL_ZERO )
            matrix[i] += tmp_recv[i];
    }
  }
  else
  {
    world.send(0, mpigetrank(), matrix);
  }
#endif
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_container::update_4pdm_contribution( std::vector< std::pair< std::vector<int>, double > >& spin_batch )
{
  //S_r^{(-1)}
  const TwoElectronArray& twoInt = v_2[0];
  int dim = a16_matrix_.dim1();
  int i,j,k,l,m,n,p,q;

  // Take into account orbital reordering
  const std::vector<int>& ro = dmrginp.reorder_vector();

  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    // Store significant elements only
    double value = it->second;
    if ( abs(value) < NUMERICAL_ZERO ) continue;
    // Spin indices
    int i = (it->first)[0]; int j = (it->first)[1]; int k = (it->first)[2]; int l = (it->first)[3];
    int m = (it->first)[4]; int n = (it->first)[5]; int p = (it->first)[6]; int q = (it->first)[7];

  


    {
    for ( int a=0; a<dim; ++a )
    {

      a16_matrix_(q,j,i,a,p,n) += twoInt(2*k,2*l,2*a,2*m) * value;

      a16_matrix_(p,k,j,i,n,a) -= twoInt(2*a,2*l,2*q,2*m) * value;

      a16_matrix_(p,k,j,i,a,q) -= twoInt(2*a,2*l,2*n,2*m) * value;



      a22_matrix_(p,n,j,i,a,q) -= twoInt(2*k,2*l,2*a,2*m) * value;

      a22_matrix_(p,q,j,a,i,n) -= twoInt(2*k,2*l,2*a,2*m) * value;

      a22_matrix_(n,q,k,j,i,a) += twoInt(2*a,2*l,2*p,2*m) * value;
    }
    }
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_container::update_3pdm_contribution( std::vector< std::pair< std::vector<int>, double > >& spin_batch )
{
  const TwoElectronArray& twoInt = v_2[0];
  const OneElectronArray& oneInt = v_1[0];
//  const std::map<TwoPerturbType,PerturbTwoElectronArray>& vp_2 = vpt_2;
//  const OneElectronArray& vp_1 = vpt_1;
  int dim = a16_matrix_.dim1();
  int i,j,k,l,m,n;

  // Take into account orbital reordering
  const std::vector<int>& ro = dmrginp.reorder_vector();

  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    double value = it->second;
    // Store significant elements only
    if ( abs(value) < NUMERICAL_ZERO ) continue;

    i = (it->first)[0]; j = (it->first)[1]; k = (it->first)[2]; l = (it->first)[3]; m = (it->first)[4]; n = (it->first)[5];

    {
    //S_r^{(-1)} subspace
    //A16 intermediate
    
      for(int a =0 ; a < dim; a++ )
      {
        a16_matrix_(n,j,i,a,m,l) += oneInt(2*k,2*a) *value;
        a16_matrix_(m,k,j,i,l,a) -= oneInt(2*a,2*n) *value;
        a16_matrix_(m,k,j,i,a,n) -= oneInt(2*a,2*l) *value;
      }
      for(int a =0 ; a < dim; a++ )
      for(int b =0 ; b < dim; b++ )
      {
        a16_matrix_(a,a,i,b,n,m) += twoInt(2*j,2*k,2*b,2*l) *value;
        a16_matrix_(a,j,i,b,m,n) += twoInt(2*a,2*k,2*b,2*l) *value;
        a16_matrix_(a,j,i,b,m,l) += twoInt(2*a,2*k,2*n,2*b) *value;
        a16_matrix_(a,a,j,i,b,n) -= twoInt(2*b,2*k,2*m,2*l) *value;
        a16_matrix_(n,j,i,a,a,b) -= twoInt(2*b,2*k,2*m,2*l) *value;
        a16_matrix_(a,a,j,i,m,b) -= twoInt(2*b,2*k,2*n,2*l) *value;
        a16_matrix_(a,j,i,a,m,b) -= twoInt(2*b,2*k,2*n,2*l) *value;
        a16_matrix_(a,j,i,a,b,n) -= twoInt(2*b,2*k,2*m,2*l) *value;
        a16_matrix_(m,k,j,i,a,b) -= twoInt(2*b,2*a,2*n,2*l) *value;
        a16_matrix_(a,k,j,i,l,b) -= twoInt(2*b,2*a,2*n,2*m) *value;
        a16_matrix_(a,k,j,i,b,n) -= twoInt(2*b,2*a,2*l,2*m) *value;
      }
    }


    //A22 in subspace S_i^{(1)}

    //one body interaction in H_v
    {
      for(int a =0 ; a < dim; a++ )
      {
      a22_matrix_(m,l,j,i,a,n) -= oneInt(2*k,2*a) *value;

      a22_matrix_(m,n,j,a,i,l) -= oneInt(2*k,2*a) *value;

      a22_matrix_(l,n,k,j,i,a) += oneInt(2*a,2*m) *value;
      }

    }
    // Two body interaction in H_v

    {

      for(int a =0 ; a < dim; a++ )
      for(int b =0 ; b < dim; b++ )
      {

      a22_matrix_(a, n, i, b, a, m) -=   twoInt(2*j,2*k,2*b,2*l) *value;
      a22_matrix_(n, l, i, a, b, m) -=   twoInt(2*k,2*j,2*b,2*a) *value;
      a22_matrix_(n, b, i, b, a, m) -=   twoInt(2*j,2*k,2*a,2*l) *value;
      a22_matrix_(m, b, j, a, i, n) -=   twoInt(2*b,2*k,2*a,2*l) *value;
      a22_matrix_(m, b, j, a, i, l) -=   twoInt(2*b,2*k,2*n,2*a) *value;
      a22_matrix_(a, n, j, b, i, m) -=   twoInt(2*a,2*k,2*b,2*l) *value;
      a22_matrix_(a, n, j, b, i, l) -=   twoInt(2*a,2*k,2*m,2*b) *value;
      a22_matrix_(a, m, i, a, b, n) -=   twoInt(2*j,2*k,2*b,2*l) *value;
      a22_matrix_(m, b, j, i, a, n) -=   twoInt(2*k,2*b,2*a,2*l) *value;
      a22_matrix_(m, b, j, i, a, n) += 2*twoInt(2*b,2*k,2*a,2*l) *value;
      a22_matrix_(n, b, i, a, b, m) += 2*twoInt(2*j,2*k,2*a,2*l) *value;
      a22_matrix_(a, m, j, i, b, n) -=   twoInt(2*a,2*k,2*b,2*l) *value;
      a22_matrix_(a, l, j, i, b, n) -=   twoInt(2*k,2*a,2*b,2*m) *value;
      a22_matrix_(m, b, j, i, b, a) -= 2*twoInt(2*a,2*k,2*n,2*l) *value;
      a22_matrix_(a, m, j, i, a, b) +=   twoInt(2*b,2*k,2*n,2*l) *value;
      a22_matrix_(m, b, j, b, i, a) +=   twoInt(2*a,2*k,2*n,2*l) *value;
      a22_matrix_(a, n, j, a, i, b) +=   twoInt(2*b,2*k,2*m,2*l) *value;
      a22_matrix_(l, b, k, j, i, a) +=   twoInt(2*a,2*b,2*m,2*n) *value;
      a22_matrix_(a, n, k, j, i, b) +=   twoInt(2*b,2*a,2*m,2*l) *value;
      }
    }
    

  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_container::update_2pdm_contribution( std::vector< std::pair< std::vector<int>, double > >& spin_batch )
{
  const TwoElectronArray& twoInt = v_2[0];
  const OneElectronArray& oneInt = v_1[0];
//  const std::map<TwoPerturbType,PerturbTwoElectronArray>& vp_2 = vpt_2;
//  const OneElectronArray& vp_1 = vpt_1;
  int dim = a16_matrix_.dim1();
  int i,j,k,l;

  // Take into account orbital reordering
  const std::vector<int>& ro = dmrginp.reorder_vector();

  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    double value = it->second;
    // Store significant elements only
    if ( abs(value) < NUMERICAL_ZERO ) continue;

    i = (it->first)[0]; j = (it->first)[1]; k = (it->first)[2]; l = (it->first)[3];

    //Subspace S_r^{(-1)}
    //a16 in subspace S_r^{(-1)}
    {
      for(int a=0; a< dim; a++)
      for(int b=0; b< dim; b++)
      {
        a16_matrix_(a,j,i,b,k,l) += oneInt(2*a,2*b) *value;
        a16_matrix_(a,a,i,b,l,k) += oneInt(2*j,2*b) *value;
        a16_matrix_(a,a,j,i,k,b) -= oneInt(2*b,2*l) *value;
        a16_matrix_(a,a,j,i,b,l) -= oneInt(2*b,2*k) *value;
        a16_matrix_(l,j,i,a,a,b) -= oneInt(2*b,2*k) *value;
        a16_matrix_(a,j,i,a,k,b) -= oneInt(2*b,2*l) *value;
        a16_matrix_(a,j,i,a,b,l) -= oneInt(2*b,2*k) *value;
      }
      for(int a=0; a< dim; a++)
      for(int b=0; b< dim; b++)
      for(int c=0; c< dim; c++)
      {
        a16_matrix_(a,a,i,b,b,c) -= twoInt(2*c,2*j,2*l,2*k) *value;
        a16_matrix_(a,a,j,i,b,c) -= twoInt(2*c,2*b,2*l,2*k) *value;
        a16_matrix_(a,j,i,b,b,c) -= twoInt(2*c,2*a,2*k,2*l) *value;
        a16_matrix_(a,j,i,a,b,c) -= twoInt(2*c,2*b,2*l,2*k) *value;
      }


    }

    //A22 in S_i^{(1)} subspace
    //One body interation in H_v
    {
      for( int a=0; a < dim; a++)
      for( int b=0; b < dim; b++)
      {

      a22_matrix_(k, b, j, i, a, l) += 2*oneInt(2*b,2*a) *value;
      a22_matrix_(l, b, i, a, b, k) += 2*oneInt(2*j,2*a) *value;
      a22_matrix_(a, k, j, i, b, l) -=   oneInt(2*a,2*b) *value;
      a22_matrix_(a, l, i, b, a, k) -=   oneInt(2*j,2*b) *value;
      a22_matrix_(k, b, j, a, i, l) -=   oneInt(2*b,2*a) *value;
      a22_matrix_(a, l, j, b, i, k) -=   oneInt(2*a,2*b) *value;
      a22_matrix_(l, b, i, b, a, k) -=   oneInt(2*j,2*a) *value;
      a22_matrix_(a, k, i, a, b, l) -=   oneInt(2*j,2*b) *value;
      a22_matrix_(k, b, j, i, b, a) -= 2*oneInt(2*a,2*l) *value;
      a22_matrix_(a, k, j, i, a, b) +=   oneInt(2*b,2*l) *value;
      a22_matrix_(k, b, j, b, i, a) +=   oneInt(2*a,2*l) *value;
      a22_matrix_(a, l, j, a, i, b) +=   oneInt(2*b,2*k) *value;
      }

    }
    //Two body interation in H_v
    {
      for( int a=0; a < dim; a++)
      for( int b=0; b < dim; b++)
      for( int c=0; c < dim; c++)
      {
      a22_matrix_(l, b, i, a, c, k) += 2*twoInt(2*b,2*j,2*c,2*a) *value;
      a22_matrix_(a, b, i, a, c, l) += 2*twoInt(2*b,2*j,2*c,2*k) *value;
      a22_matrix_(a, b, j, i, c, l) += 2*twoInt(2*b,2*a,2*c,2*k) *value;
      a22_matrix_(a, b, j, i, c, l) -=   twoInt(2*a,2*b,2*c,2*k) *value;
      a22_matrix_(a, l, i, b, c, k) -=   twoInt(2*a,2*j,2*c,2*b) *value;
      a22_matrix_(a, b, i, b, c, l) -=   twoInt(2*a,2*j,2*c,2*k) *value;
      a22_matrix_(l, b, i, a, c, k) -=   twoInt(2*j,2*b,2*c,2*a) *value;
      a22_matrix_(a, b, i, a, c, l) -=   twoInt(2*j,2*b,2*c,2*k) *value;
      a22_matrix_(a, k, i, b, c, l) -=   twoInt(2*j,2*a,2*c,2*b) *value;
      a22_matrix_(a, b, i, b, c, k) -=   twoInt(2*j,2*a,2*c,2*l) *value;
      a22_matrix_(a, b, i, c, b, l) += 2*twoInt(2*a,2*j,2*c,2*k) *value;
      a22_matrix_(a, b, i, c, b, k) += 2*twoInt(2*a,2*j,2*l,2*c) *value;
      a22_matrix_(a, b, i, c, a, l) -=   twoInt(2*b,2*j,2*c,2*k) *value;
      a22_matrix_(a, b, i, c, a, k) -=   twoInt(2*b,2*j,2*l,2*c) *value;
      a22_matrix_(a, b, j, c, i, l) -=   twoInt(2*b,2*a,2*c,2*k) *value;
      a22_matrix_(a, b, j, c, i, k) -=   twoInt(2*b,2*a,2*l,2*c) *value;
      a22_matrix_(a, b, i, a, b, c) -= 2*twoInt(2*c,2*j,2*l,2*k) *value;
      a22_matrix_(a, b, j, i, b, c) -= 2*twoInt(2*c,2*a,2*l,2*k) *value;
      a22_matrix_(a, b, i, b, a, c) +=   twoInt(2*c,2*j,2*l,2*k) *value;
      a22_matrix_(a, b, j, i, a, c) +=   twoInt(2*c,2*b,2*l,2*k) *value;
      a22_matrix_(a, b, j, b, i, c) +=   twoInt(2*c,2*a,2*l,2*k) *value;
      a22_matrix_(a, b, j, a, i, c) +=   twoInt(2*c,2*b,2*k,2*l) *value;

      }

    }

  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_container::update_1pdm_contribution( std::vector< std::pair< std::vector<int>, double > >& spin_batch )
{

  const TwoElectronArray& twoInt = v_2[0];
  const OneElectronArray& oneInt = v_1[0];
  int dim = a16_matrix_.dim1();
  int i,j;


  // Take into account orbital reordering
  const std::vector<int>& ro = dmrginp.reorder_vector();

  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    double value = it->second;
    // Store significant elements only
    if ( abs(value) < NUMERICAL_ZERO ) continue;
    int i = (it->first)[0];
    int j = (it->first)[1];
    //Subspace S_r^{(-1)}

    {

    for(int aa = 0; aa <dim; aa++)
    {

      for(int a =0; a <dim; a++)
      for(int c =0; c <dim; c++)
        a16_matrix_(aa,aa,i,a,a,c) -=oneInt(2*c,2*j) *value;

    }
    }
    //A22 for S_i^{(1)} subspace
    {

      for( int a=0; a < dim; a++)
      for( int b=0; b < dim; b++)
      for( int c=0; c < dim; c++)
      {
      a22_matrix_(a,b,i,a,c,j) += 2*oneInt(2*b,2*c) *value;
      a22_matrix_(a,b,i,b,c,j) -= oneInt(2*a,2*c) *value;
      a22_matrix_(a,b,i,c,b,j) += 2*oneInt(2*a,2*c) *value;
      a22_matrix_(a,b,i,c,a,j) -= oneInt(2*b,2*c) *value;
      a22_matrix_(a,b,i,a,b,c) -= 2*oneInt(2*c,2*j) *value;
      a22_matrix_(a,b,i,b,a,c) += oneInt(2*c,2*j) *value;
      }

      for( int a=0; a < dim; a++)
      for( int b=0; b < dim; b++)
      for( int c=0; c < dim; c++)
      for( int d=0; d < dim; d++)
      {
      a22_matrix_(a,b,i,c,d,j) += 2*twoInt(2*b,2*a,2*d,2*c) *value;
      a22_matrix_(a,b,i,c,d,j) -= twoInt(2*a,2*b,2*d,2*c) *value;
      }


    }
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_container::store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements)
{

  std::vector< std::pair< std::vector<int>, double > > spatial_batch;

  if ( new_spin_orbital_elements[0].first.size() == 2 ) {
    // 1PDM elements
    Onepdm_permutations perm;
    perm.get_spatial_batch( new_spin_orbital_elements, spatial_batch );
    update_1pdm_contribution( spatial_batch );
  }
  else if ( new_spin_orbital_elements[0].first.size() == 4 ) {
    // 2PDM elements
    Twopdm_permutations perm;
    perm.get_spatial_batch( new_spin_orbital_elements, spatial_batch );
    update_2pdm_contribution( spatial_batch );
  }
  else if ( new_spin_orbital_elements[0].first.size() == 6 ) {
    // 3PDM elements
    Threepdm_permutations perm;
    perm.get_spatial_batch( new_spin_orbital_elements, spatial_batch );
    update_3pdm_contribution( spatial_batch );
  }
  else if ( new_spin_orbital_elements[0].first.size() == 8 ) {
    // 4PDM elements
    Fourpdm_permutations perm;
    perm.get_spatial_batch( new_spin_orbital_elements, spatial_batch );
    update_4pdm_contribution( spatial_batch );
  }
  else abort();

}

//===========================================================================================================================================================

}
}

