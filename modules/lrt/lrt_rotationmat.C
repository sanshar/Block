/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "rotationmat.h"
#include "pario.h"
#include "MatrixBLAS.h"
#include <include/sortutils.h>
#include <boost/serialization/vector.hpp>
#include "pario.h"

#include "modules/lrt/lrt_rotationmat.h"

using namespace boost;
using namespace std;

// compute selected weights and rejected basis and weights in addition to selectedbasis (rotatematrix)
double SpinAdapted::LRT::assign_matrix_by_dm
(std::vector<DiagonalMatrix>& eigenmatrix,
 std::vector<Matrix>& rotatematrix, std::vector< std::vector<double> >& selectedwts,
 std::vector<Matrix>& rejectedbasis, std::vector< std::vector<double> >& rejectedwts,
 SparseMatrix& transformmatrix, vector<pair<int, int> >& inorderwts, vector<vector<int> >& wtsbyquanta,
 int totalstatesbydm, int totalstatesbyquanta, int left_block_size, int right_block_size)
{
  const int min_states = totalstatesbydm;

  if (dmrginp.outputlevel() > 0)
    pout << " \t\t\t assigning a total of " << min_states << " states using the dm alone " << endl;
  double totalnorm = 0.;
  rotatematrix.resize(eigenmatrix.size());
  selectedwts.resize(eigenmatrix.size());

  for (int i = 0; i < totalstatesbydm; ++i)
    {
      int q = inorderwts[i].first;
      int qs = inorderwts[i].second;

      //if(i < min_states && eigenmatrix[q].element(qs, qs) > 1.e-12)
      if( eigenmatrix[q].element(qs, qs) > dmrginp.density_tol())
      {
        if (rotatematrix[q].Ncols() == 0)
        {
	        rotatematrix[q] = transformmatrix(q, q).Column(qs + 1);
        }
        else
        {
	        rotatematrix[q] |= transformmatrix(q, q).Column(qs + 1);      
        }
        vector<int>::iterator findit = find(wtsbyquanta[q].begin(), wtsbyquanta[q].end(), qs);
        if (findit == wtsbyquanta[q].end()) { pout << " error in assign matrix " << endl; abort(); }
        wtsbyquanta[q].erase(findit);
        totalnorm += eigenmatrix[q].element(qs, qs);
        selectedwts[q].push_back(eigenmatrix[q].element(qs, qs));
      }
    }

  if (dmrginp.outputlevel() > 0)
    pout << " \t\t\t assigning a total of " << totalstatesbyquanta << " states using quanta selection " << " for a norm of " << totalnorm << endl;

  int assignedbyq = 0;

  int nquanta = rotatematrix.size();
 
  int totalstatesleft = 0;
  for (int i = 0; i < nquanta; ++i)
    {
      totalstatesleft += wtsbyquanta[i].size();
    }

  if (dmrginp.outputlevel() > 0)
    pout << " \t\t\t a total of " << totalstatesleft << " to be assigned " << endl;
  
  // now sort quanta in order of importance
  vector<double> totalquantaweights(nquanta);
  for (int q = 0; q < totalquantaweights.size(); ++q)
  {
    for (int qs = 0; qs < wtsbyquanta[q].size(); ++qs)
	    totalquantaweights[q] += eigenmatrix[q].element(qs, qs);
  }
  vector<int> inorderquanta(nquanta);
  sort_data_to_indices(totalquantaweights, inorderquanta);
  reverse(inorderquanta.begin(), inorderquanta.end());  


  // reorder modified wtsbyquanta into a usable form
  vector<pair<int, int> > linearwtsbyquanta;
  int qspointer = 0;
  
  while(totalstatesleft)
  {
    for (int i = 0; i < nquanta; ++i)
	  {
	    int q = inorderquanta[i];
	    if (qspointer < wtsbyquanta[q].size())
	    {
	      linearwtsbyquanta.push_back(make_pair(q, wtsbyquanta[q][qspointer]));
	      --totalstatesleft;
	    }
	  }
    ++qspointer;
  }

  for (int i = 0; i < totalstatesbyquanta; ++i)
  {
    int q = linearwtsbyquanta[i].first;
    int qs = linearwtsbyquanta[i].second;
    if( eigenmatrix[q].element(qs, qs) > dmrginp.density_tol())
    {
      if (rotatematrix[q].Ncols() == 0)
      {
	      rotatematrix[q] = transformmatrix(q, q).Column(qs + 1);
      }
      else
      {
	      rotatematrix[q] |= transformmatrix(q, q).Column(qs + 1);      
      }
      vector<int>::iterator findit = find(wtsbyquanta[q].begin(), wtsbyquanta[q].end(), qs);
      if (findit == wtsbyquanta[q].end()) { pout << " error in assign matrix " << endl; abort(); }
      wtsbyquanta[q].erase(findit);
      totalnorm += eigenmatrix[q].element(qs, qs);
      selectedwts[q].push_back(eigenmatrix[q].element(qs, qs));
    }
  }
  
  rejectedbasis.resize(rotatematrix.size());
  rejectedwts.resize(rotatematrix.size());
  for(int q = 0; q < wtsbyquanta.size(); ++q) {
    for(vector<int>::iterator it = wtsbyquanta[q].begin(); it !=  wtsbyquanta[q].end(); ++it) {
      if(rejectedbasis[q].Ncols() == 0)
        rejectedbasis[q] = transformmatrix(q, q).Column(*it + 1);
      else
        rejectedbasis[q] |= transformmatrix(q, q).Column(*it + 1);
      rejectedwts[q].push_back(eigenmatrix[q].element(*it, *it));
    }
  }

  double norm = 0.;
  for(int i=0;i<eigenmatrix.size();++i)
    for(int j=0;j<eigenmatrix[i].Nrows();++j)
      norm += eigenmatrix[i].element(j, j);
  if (dmrginp.outputlevel() > 0)
    pout << " \t\t\t total norm: " << norm <<"  norm after truncation: "<<totalnorm<< endl;

  return norm-totalnorm;

  //return (1. - totalnorm/norm);
}

//from Jon's code
/****************************************************************************************************
double SpinAdapted::LRT::assign_matrix_by_dm_deriv
(const std::vector<Matrix>& rotatematrix, const std::vector< std::vector<double> >& selectedwts,
 const std::vector<Matrix>& rejectedbasis, const std::vector< std::vector<double> >& rejectedwts,
 const SparseMatrix& density_deriv, std::vector<Matrix>& rotatematrix_deriv)
{
  int nquanta=rotatematrix.size();
  rotatematrix_deriv.resize(nquanta);

  for (int q=0;q<nquanta;++q)
  {
    Matrix matrix_elements;
    if(rotatematrix[q].Ncols() > 0 && rejectedbasis[q].Ncols() > 0) {
//    matrix_elements = (rotatematrix[q].t()*perturbeddm(q,q))*rejectedbasis[q];
      Matrix tmp_elements;
      MatrixMultiply(rotatematrix[q], 't', perturbeddm(q, q), 'n', tmp_elements, 1.0);
      MatrixMultiply(tmp_elements, 'n', rejectedbasis[q], 'n', matrix_elements, 1.0);
    }
    for (int i=0;i<rotatematrix[q].Ncols();++i)
    {
      ColumnVector perturbedvec(rotatematrix[q].Nrows());
      perturbedvec = 0.;
      for (int a=0;a<rejectedbasis[q].Ncols();++a)
      {
        ColumnVector rejectedvec=rejectedbasis[q].Column(a+1);
        double matrix_element_ia = matrix_elements.element(i,a);
        if (abs(rejectedwts[q][a]-selectedwts[q][i]) > 1.e-12)
          perturbedvec += rejectedvec*(matrix_element_ia/(selectedwts[q][i]-rejectedwts[q][a]));
      }
      if(rotatematrix_deriv[q].Ncols() == 0)
        rotatematrix_deriv[q] = perturbedvec;
      else
        rotatematrix_deriv[q] |= perturbedvec;
    }
  }

  return 0.0;
}
*****************************************************************************************************/

// alternative code not using rejectedbasis
double SpinAdapted::LRT::assign_matrix_by_dm_deriv
(const std::vector<Matrix>& rotatematrix,
 const std::vector< std::vector<double> >& selectedwts, const std::vector< std::vector<double> >& rejectedwts,
 const SparseMatrix& density_deriv, std::vector<Matrix>& rotatematrix_deriv, bool projection)
{
  int nquanta = rotatematrix.size();
  rotatematrix_deriv.resize(nquanta);
  for(int q = 0; q < nquanta; ++q) {
    if(rotatematrix[q].Ncols() > 0) {
      Matrix matrix_elements(rotatematrix[q]); matrix_elements = 0.0;
      MatrixMultiply(density_deriv(q, q), 'n', rotatematrix[q], 'n', matrix_elements, 1.0);
      if(rejectedwts[q].size() > 0) {
        for(int i = 0; i < rotatematrix[q].Ncols(); ++i) {
          ColumnVector derived_basis(rotatematrix[q].Nrows());
          derived_basis = 0.0;
          // FIXME: in case M is large enough, small selected weights introduce numerical instability

          if (abs(selectedwts[q][i]) > dmrginp.density_tol()) {
            derived_basis = matrix_elements.Column(i + 1)/selectedwts[q][i];
          }
//        else if(dmrginp.outputlevel() > 0) {
          else {
            pout << "Warning: ignored small selected weight (" << scientific << selectedwts[q][i] << ") meaning that 1-st order rotation matrix is not exact" << endl;
          }

          if(rotatematrix_deriv[q].Ncols() == 0)
            rotatematrix_deriv[q] = derived_basis;
          else
            rotatematrix_deriv[q] |= derived_basis;
        }
        if(projection) {
          Matrix projected_matrix(rotatematrix[q].Ncols(), rotatematrix[q].Ncols()); projected_matrix = 0.0;
          MatrixMultiply(rotatematrix[q], 't', rotatematrix_deriv[q], 'n', projected_matrix, 1.0);
          MatrixMultiply(rotatematrix[q], 'n', projected_matrix, 'n', rotatematrix_deriv[q],-1.0);
        }
      }
      else {
        rotatematrix_deriv[q] = matrix_elements;
      }
    }
  }
  return 0.0;
}

void SpinAdapted::LRT::project_onto_rejectedspace
(const Wavefunction& c, const std::vector<Matrix>& rotatematrix, const bool& dot_with_sys, Wavefunction& c_projected)
{
  int nquanta = rotatematrix.size();
  c_projected = c;
  // FIXME: should check for dot_with_sys, supposed to be false with LRT
  for(int i = 0; i < c_projected.nrows(); ++i) {
    for(int j = 0; j < c_projected.ncols(); ++j) {
      if(c_projected.allowed(i, j)) {
        const Matrix& lM = rotatematrix[i];
              Matrix& nM = c_projected.operator_element(i, j);        
              Matrix  tM = nM;
              tM.ReSize(lM.Ncols(), nM.Ncols()); tM = 0.0;
        MatrixMultiply(lM, 't', nM, 'n', tM, 1.0); // T = L^(t) * C
        MatrixMultiply(lM, 'n', tM, 'n', nM,-1.0); // C = C - L * T
      }
    }
  }
}

