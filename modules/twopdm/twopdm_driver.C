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

#include "twopdm_driver.h"
//#include "npdm_patterns.h"
//#include "npdm_expectations.h"
//#include "npdm_operator_wrappers.h"

namespace SpinAdapted{

//===========================================================================================================================================================

void Twopdm_driver::save_npdm_text(const int &i, const int &j)
{
  if(!mpigetrank())
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/twopdm.", i, j);
    ofstream ofs(file);
    ofs << twopdm.dim1() << endl;
    for(int k=0;k<twopdm.dim1();++k)
      for(int l=0;l<twopdm.dim2();++l)
        for(int m=0;m<twopdm.dim3();++m)
          for(int n=0;n<twopdm.dim4();++n)
            ofs << boost::format("%d %d %d %d %20.14e\n") % k % l % m % n % twopdm(k,l,m,n);
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_driver::save_spatial_npdm_text(const int &i, const int &j)
{
  //the spatial has a factor of 1/2 in front of it 
  if(!mpigetrank())
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/spatial_twopdm.", i, j);
    ofstream ofs(file);
    ofs << twopdm.dim1()/2 << endl;
    for(int k=0;k<twopdm.dim1()/2;++k)
      for(int l=0;l<twopdm.dim2()/2;++l)
        for(int m=0;m<twopdm.dim3()/2;++m)
          for(int n=0;n<twopdm.dim4()/2;++n) {
	          double pdm = 0.0;
      	    for (int s=0; s<2; s++)
      	      for (int t =0; t<2; t++) {
//                 if ( (k==0) && (l==0) && (m==1) && (n==2) ) 
//                 if ( (k==0) && (l==0) && (m==2) && (n==1) ) 
//                     std::cout << 2*k+s<<","<< 2*l+t<<","<< 2*m+t<<","<< 2*n+s<<"\t\t"<<pdm<<"\t"<<twopdm(2*k+s, 2*l+t, 2*m+t, 2*n+s)*0.5 <<std::endl;
                 pdm += twopdm(2*k+s, 2*l+t, 2*m+t, 2*n+s)*0.5;
               }
		
             ofs << boost::format("%d %d %d %d %20.14e\n") % k % l % m % n % pdm;
	  }
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_driver::save_spatial_npdm_binary(const int &i, const int &j)
{
  //the spatial has a factor of 1/2 in front of it 
  if(!mpigetrank())
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/spatial_binary_twopdm.", i, j);
    FILE* f = fopen(file, "wb");

    int nrows = twopdm.dim1()/2;
    array_4d<double> pdm; pdm.resize(nrows, nrows, nrows, nrows);
    for(int k=0;k<twopdm.dim1()/2;++k)
      for(int l=0;l<twopdm.dim2()/2;++l)
        for(int m=0;m<twopdm.dim3()/2;++m)
          for(int n=0;n<twopdm.dim4()/2;++n) {
	    pdm(k, l, m, n) = 0.0;
	    for (int s=0; s<2; s++)
	      for (int t =0; t<2; t++)
		pdm(k, l, m, n) += twopdm(2*k+s, 2*l+t, 2*m+t, 2*n+s)*0.5;
	  }
    int result = fwrite(&nrows,  1, sizeof(int), f);
    result = fwrite(&pdm(0,0,0,0), pdm.size(), sizeof(double), f);
    fclose(f);
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_driver::save_npdm_binary(const int &i, const int &j)
{
  if(!mpigetrank())
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/twopdm.", i, j);
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << twopdm;
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_driver::load_npdm_binary(const int &i, const int &j)
{
  if(!mpigetrank())
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/twopdm.", i, j);
    std::ifstream ifs(file, std::ios::binary);
    boost::archive::binary_iarchive load(ifs);
    load >> twopdm;
    ifs.close();
  }
#ifndef SERIAL
  mpi::communicator world;
  mpi::broadcast(world,twopdm,0);
  if(mpigetrank())
    twopdm.Clear();
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//
//void Twopdm_driver::save_averaged_twopdm(const int &nroots)
//{
//  if(!mpigetrank())
//  {
//    array_4d<double> twopdm;
//    char file[5000];
//    for(int i=0;i<nroots;++i)
//    {
//      sprintf (file, "%s%s%d.%d", dmrginp.save_prefix().c_str(),"/twopdm.", i, i);
//      ifstream ifs(file);
//      int size = 0;
//      ifs >> size;
//      if(i==0)
//        twopdm.resize(size,size,size,size);
//      int k,l,m,n;
//      double val;
//      while(ifs >> k)
//      {
//        ifs >> l >> m >> n >> val;
//        twopdm(k,l,m,n) += dmrginp.weights()[i]*val;
//      }
//    }
//    sprintf (file, "%s%s", dmrginp.save_prefix().c_str(),"/twopdm");
//    ofstream ofs(file);
//    ofs << twopdm.dim1() << endl;
//    for(int k=0;k<twopdm.dim1();++k)
//      for(int l=0;l<twopdm.dim2();++l)
//        for(int m=0;m<twopdm.dim3();++m)
//          for(int n=0;n<twopdm.dim4();++n)
//            ofs << boost::format("%d %d %d %d %20.14e\n") % k % l % m % n % twopdm(k,l,m,n);
//
//    ofs.close();
//  }
//}
//
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_driver::accumulate_npdm()
{
#ifndef SERIAL
  array_4d<double> tmp_recv;
  mpi::communicator world;
  if (!mpigetrank())
    {
      for(int i=1;i<world.size();++i)
	{
	  world.recv(i, i, tmp_recv);
	  for(int k=0;k<twopdm.dim1();++k)
	    for(int l=0;l<twopdm.dim2();++l)
	      for(int m=0;m<twopdm.dim3();++m)
		for(int n=0;n<twopdm.dim4();++n)
		  if(tmp_recv(k,l,m,n) != 0.)
		    twopdm(k,l,m,n) = tmp_recv(k,l,m,n);
	}
    }
  else
    {
      world.send(0, mpigetrank(), twopdm);
    }
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_driver::assign_npdm_antisymmetric(const int i, const int j, const int k, const int l, const double val)
{

//MAW
//if ( abs(val) > 1e-8 ) pout << "so-twopdm val: i,j,k,l = " << i << "," << j << "," << k << "," << l << "\t\t" << val << endl;
//pout << "so-twopdm val: i,j,k,l = " << i << "," << j << "," << k << "," << l << "\t\t" << val << endl;

  // Test for duplicates
  if ( twopdm(i, j, k, l) != 0.0 && abs(twopdm(i,j,k,l)-val) > 1e-6)
    {
      void *array[10];
      size_t size;
      size = backtrace(array, 10);
      cout << "WARNING: Already calculated "<<i<<" "<<j<<" "<<k<<" "<<l<<endl;
      //backtrace_symbols_fd(array, size, 2);
      cout << "earlier value: "<<twopdm(i,j,k,l)<<endl<< "new value:     "<<val<<endl;
      assert( false );
      return;
    }

  twopdm(i, j, k, l) = val;
  twopdm(i, j, l, k) = -val;
  twopdm(j, i, k, l) = -val;
  twopdm(j, i, l, k) = val;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_driver::assign_npdm_elements( std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements)
{
  for (int i=0; i < new_spin_orbital_elements.size(); ++i) {
    int ix = new_spin_orbital_elements[i].first[0];
    int jx = new_spin_orbital_elements[i].first[1];
    int kx = new_spin_orbital_elements[i].first[2];
    int lx = new_spin_orbital_elements[i].first[3];
    double x = new_spin_orbital_elements[i].second;
    assign_npdm_antisymmetric(ix, jx, kx, lx, x);
  }

//FIXME is the transpose always needed?
  for (int i=0; i < new_spin_orbital_elements.size(); ++i) {
    int ix = new_spin_orbital_elements[i].first[3];
    int jx = new_spin_orbital_elements[i].first[2];
    int kx = new_spin_orbital_elements[i].first[1];
    int lx = new_spin_orbital_elements[i].first[0];
    double x = new_spin_orbital_elements[i].second;
    assign_npdm_antisymmetric(ix, jx, kx, lx, x);
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Twopdm_driver::calcenergy(int state)
{
  using namespace SpinAdapted;
  Matrix onepdm(2*dmrginp.last_site(), 2*dmrginp.last_site()); onepdm = 0.0;
  int nelec = dmrginp.real_particle_number();

  for (int i=0; i<dmrginp.last_site()*2; i++)
  for (int j=0; j<dmrginp.last_site()*2; j++)
  for (int k=0; k<dmrginp.last_site()*2; k++)
    onepdm(i+1, j+1) += twopdm(i, k, k, j);

  onepdm /= (nelec-1);
  double nel = 0.0, sz=0.0;
  for (int i=0; i<dmrginp.last_site(); i++) {
    nel += onepdm(2*i+1, 2*i+1)+onepdm(2*i+2, 2*i+2);
    sz += onepdm(2*i+1, 2*i+1)-onepdm(2*i+2, 2*i+2);
  }

  double energy = 0.0;
  for (int i=0; i<dmrginp.last_site()*2; i++)
  for (int j=0; j<dmrginp.last_site()*2; j++)
  for (int k=0; k<dmrginp.last_site()*2; k++)
  for (int l=0; l<dmrginp.last_site()*2; l++)
    energy += v_2(i,j,k,l)*twopdm(i,j,l,k);

  energy *= 0.5;

  for (int i=0; i<dmrginp.last_site()*2; i++)
  for (int j=0; j<dmrginp.last_site()*2; j++)
    energy += v_1(i,j) * onepdm(i+1,j+1);

  pout << "energy of state "<< state <<" = "<< energy+dmrginp.get_coreenergy()<<endl;

  ofstream out("onepdm_fromtpdm");
  for (int i=0; i<dmrginp.last_site()*2; i++)
  for (int j=0; j<dmrginp.last_site()*2; j++)
    out<<i<<"  "<<j<<"  "<<onepdm(i+1,j+1)<<endl;
  out.close();
  
}

//===========================================================================================================================================================

}

