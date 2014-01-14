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
#include "threepdm_driver.h"

namespace SpinAdapted{

//===========================================================================================================================================================

void Threepdm_driver::save_npdm_text(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/threepdm.", i, j,".txt");
    ofstream ofs(file);
    ofs << threepdm.dim1() << endl;

    double trace = 0.0;
    for(int i=0; i<threepdm.dim1(); ++i)
      for(int j=0; j<threepdm.dim2(); ++j)
        for(int k=0; k<threepdm.dim3(); ++k)
          for(int l=0; l<threepdm.dim4(); ++l)
            for(int m=0; m<threepdm.dim5(); ++m)
              for(int n=0; n<threepdm.dim6(); ++n) {
                if ( abs(threepdm(i,j,k,l,m,n)) > 1e-14 ) {
                  ofs << boost::format("%d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % threepdm(i,j,k,l,m,n);
                  if ( (i==n) && (j==m) && (k==l) ) trace += threepdm(i,j,k,l,m,n);
                }
              }
    ofs.close();
    std::cout << "Spin-orbital 3PDM trace = " << trace << "\n";
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_driver::save_spatial_npdm_text(const int &i, const int &j)
{
//FIXME  2pdm spatial has a factor of 1/2 in front of it  ???
  double factor = 1.0;
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.", i, j,".txt");
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

void Threepdm_driver::save_spatial_npdm_binary(const int &i, const int &j)
{
//FIXME  2pdm spatial has a factor of 1/2 in front of it  ???
  double factor = 1.0;
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.", i, j,".bin");
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

void Threepdm_driver::save_npdm_binary(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/threepdm.", i, j,".bin");
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << threepdm;
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_driver::load_npdm_binary(const int &i, const int &j)
{
pout << "load_threepdm_binary\n";
cout.flush();
assert(false); // <<<< CAN WE RETHINK USE OF DISK FOR NPDM?
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/threepdm.", i, j,".bin");
    std::ifstream ifs(file, std::ios::binary);
    boost::archive::binary_iarchive load(ifs);
cout << "loading... load_threepdm_binary\n";
cout.flush();
    load >> threepdm;
    ifs.close();
  }

//cout << threepdm(0,0,0,0,0,0) << endl;
//cout << threepdm(1,2,3,4,5,6) << endl;
//#ifndef SERIAL
//cout << "broadcasting1... load_threepdm_binary\n";
//cout.flush();
//
//assert(false); // MEMORY USE LARGE HERE!!!
//  mpi::communicator world;
//  mpi::broadcast(world,threepdm,0);
//cout << "broadcasting2... load_threepdm_binary\n";
//cout.flush();
//  if( mpigetrank() != 0)
//    threepdm.Clear();
//#endif

cout << "all done... load_threepdm_binary\n";
cout.flush();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//
//void Threepdm_driver::save_averaged_npdm(const int &nroots)
//{
/////  NYI
//  assert(false);
//}
//
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_driver::accumulate_npdm()
{
#ifndef SERIAL
  array_6d<double> tmp_recv;
  mpi::communicator world;
//cout << "threepdm size = " << threepdm.get_size() << " ; rank = " << mpigetrank() << endl;
//cout.flush();
//threepdm.resize(26,26,26,26,26,26);
  if( mpigetrank() == 0)
  {
    for(int p=1; p<world.size(); ++p) {
//assert(false); // MEMORY USE LARGE HERE!!!  Is this why it crashes sometimes?
      world.recv(p, p, tmp_recv);
      for(int i=0; i<threepdm.dim1(); ++i)
        for(int j=0; j<threepdm.dim2(); ++j)
          for(int k=0; k<threepdm.dim3(); ++k)
            for(int l=0; l<threepdm.dim4(); ++l)
              for(int m=0; m<threepdm.dim5(); ++m)
                for(int n=0; n<threepdm.dim6(); ++n) 
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

void Threepdm_driver::assign_npdm_antisymmetric(const int i, const int j, const int k, const int l, const int m, const int n, const double val)
{
if ( abs(val) > 1e-8 ) {
  cout << "so-threepdm val: i,j,k,l,m,n = " 
       << i << "," << j << "," << k << "," << l << "," << m << "," << n
       << "\t\t" << val << endl;
}

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

  // If indices are not all unique, then all elements should be zero
  std::vector<int> v = {i,j,k};
  std::sort( v.begin(), v.end() );
  if ( (v[0]==v[1]) || (v[1]==v[2]) ) { if (abs(val) > 1e-15) { std::cout << abs(val) << std::endl; assert(false); } }
  std::vector<int> w = {l,m,n};
  std::sort( w.begin(), w.end() );
  if ( (w[0]==w[1]) || (w[1]==w[2]) ) { if (abs(val) > 1e-15) { std::cout << abs(val) << std::endl; assert(false); } }


  // The number of possible combinations is (3!)**2  (For 4pdm and higher we use a general implementation)
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

void Threepdm_driver::assign_npdm_elements( std::vector< std::pair<std::vector<int>, double >>& new_spin_orbital_elements )
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
    assign_npdm_antisymmetric(ix, jx, kx, lx, mx, nx, x);
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
    assign_npdm_antisymmetric(ix, jx, kx, lx, mx, nx, x);
  }
}

//===========================================================================================================================================================

}

