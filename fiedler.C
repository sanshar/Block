#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <newmat.h>
#include <newmatio.h>
#include <newmatap.h>
#include <sortutils.h>

std::vector<int> fiedler_reorder(const SymmetricMatrix& m)
{
  SymmetricMatrix absm=m;
  const int nrows=m.Nrows();
  for (int i=0;i<nrows;++i)
    for (int j=0;j<=i;++j)
       //absolute value
       absm.element(i,j)=std::fabs(absm.element(i,j));

  //laplacian
  SymmetricMatrix lap(nrows);
  lap=0.;
  for (int i=0;i<nrows;++i)
    lap.element(i,i)=absm.Row(i+1).Sum();
  lap-=absm;

  //roa
  cout << "lap " << endl;
  cout << lap; 
  
  DiagonalMatrix eigs; Matrix vecs;
  
  EigenValues(lap,eigs,vecs);

  ColumnVector fvec=vecs.Column(2);
  cout << "fvec: " << endl;
  cout << fvec << endl;
  std::vector<double> fvec_stl(nrows);
  //copies over fvec to fvec_stl
  std::copy(&fvec.element(0),&fvec.element(0)+nrows,fvec_stl.begin());
  std::vector<int> findices;
  //sorts the data by eigenvalue in ascending order
  sort_data_to_indices(fvec_stl,findices);
  for (int i=0;i<nrows;++i){
     cout << "findices " << findices[i] << " "  << fvec_stl[i] << endl;
  }
  

  std::vector<int> permindices(nrows);
  for (int i=0;i<nrows;++i){
     //orders the orbitals: orb 0 is in position 6, orb 1 is in position 7 u.s.w.
    permindices[i]=static_cast<int>(std::find(findices.begin(),findices.end(),i)-findices.begin());
    //cout << "ROA permindices " << i << " " << permindices[i]  << endl;
  }

  return permindices;
}
