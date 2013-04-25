/*                                                                           
Developed by Roberto Olivares-Amaya and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
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
  cout << m << endl;
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

  DiagonalMatrix eigs; Matrix vecs;
  
  EigenValues(lap,eigs,vecs);

  ColumnVector fvec=vecs.Column(2);
  //cout << "fvec: " << endl; //roa
  //cout << fvec << endl; //roa
  std::vector<double> fvec_stl(nrows);
  //copies over fvec to fvec_stl
  std::copy(&fvec.element(0),&fvec.element(0)+nrows,fvec_stl.begin());
  std::vector<int> findices;
  //sorts the data by eigenvalue in ascending order
  sort_data_to_indices(fvec_stl,findices);
  
  return findices;
  /* BLOCK works with findices

  std::vector<int> permindices(nrows);
  for (int i=0;i<nrows;++i){
     //orders the orbitals: orb 0 is in position 6, orb 1 is in position 7 u.s.w.
    permindices[i]=static_cast<int>(std::find(findices.begin(),findices.end(),i)-findices.begin());
    cout << "permindices findices " << permindices[i]+1 << " " << findices[i]+1 << endl;
  }
  */

  //Use permindices for internal use, findices for printout purposes
  //return permindices;

  //ROA Testing purposes
  //Matrix newmat = permute(m, &permindices[0]);
  //cout << "Permuted matrix" << endl;
  //cout << newmat << endl;

  //ROA print
  /*
  for (int i=0;i<findices.size();++i)
     cout << findices[i]+1 << " " << endl;
  cout << endl;
  */

}
