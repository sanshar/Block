/*                                                                           
Developed by Roberto Olivares-Amaya and Garnet K.-L. Chan, 2013                      
Copyright (c) 2013, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma, Garnet K.-L. Chan and Roberto Olivares-Amaya
*/

#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <newmat.h>
#include <newmatio.h>
#include <newmatap.h>
#include <sortutils.h>
#include "fiedler.h"
#ifdef UNITTEST
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Fiedler
#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>
#endif

Matrix argpermute(const Matrix& m, const int* indices)
{
  Matrix newm(m.Nrows(),m.Ncols());
  newm=0.;
  for (int i=0;i<m.Nrows();++i)
    for (int j=0;j<m.Ncols();++j)
      newm.element(i,j)=m.element(indices[i],indices[j]);
  return newm;
}

Matrix permute(const Matrix& m, const int* indices)
{
  Matrix newm(m.Nrows(),m.Ncols());
  newm=0.;
  for (int i=0;i<m.Nrows();++i){
    for (int j=0;j<m.Ncols();++j){
      newm.element(indices[i],indices[j])=m.element(i,j);
    }
  }
  return newm;
}

std::vector<int> fiedler_reorder(const SymmetricMatrix& m)
{
  SymmetricMatrix absm=m;
  const int nrows=m.Nrows();
  for (int i=0;i<nrows;++i) {
    for (int j=0;j<=i;++j){
       //absolute value
       absm.element(i,j)=std::fabs(absm.element(i,j));
    }
  }

  //laplacian
  SymmetricMatrix lap(nrows);
  lap=0.;
  for (int i=0;i<nrows;++i)
    lap.element(i,i)=absm.Row(i+1).Sum();
  lap-=absm;

  DiagonalMatrix eigs; 
  Matrix vecs;
  
  EigenValues(lap,eigs,vecs);

  ColumnVector fvec=vecs.Column(2);
  std::vector<double> fvec_stl(nrows);
  //copies over fvec to fvec_stl
  std::copy(&fvec.element(0),&fvec.element(0)+nrows,fvec_stl.begin());
  std::vector<int> findices;
  //sorts the data by eigenvalue in ascending order
  sort_data_to_indices(fvec_stl,findices);
  
  return findices;
  /* BLOCK works with findices*/

}

#ifdef UNITTEST
   // General: Evaluates fiedler_reorder for a Huckel case
   // Subroutines tested: permute, fiedler_reorder
   BOOST_AUTO_TEST_CASE(Fiedler){
      std::vector<int> reorderTest;
      int vSize=8;
      Matrix h(vSize,vSize);
      SymmetricMatrix hsym;
      h=0.;
      for (int i=0;i<vSize;++i)
         for (int j=0;j<vSize;++j)
            if (abs(i-j)==1)
               h.element(i,j)=-1.;

      //Setting up a permutation of indices
      int indices[]={1,0,5,7,4,2,3,6};
      //Permuting
      Matrix hper=permute(h,indices);
      hsym << hper;
      reorderTest = fiedler_reorder(hsym);
      std::vector<int> expected(vSize);
      for (int i=0;i<vSize;i++) expected.at(i)=indices[(vSize-1)-i];
      // Running the test
      BOOST_CHECK_EQUAL_COLLECTIONS(reorderTest.begin(), reorderTest.end(), expected.begin(), expected.end()); }
#endif
