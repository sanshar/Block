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
#include "cmath"
using namespace boost;
using namespace std;


void SpinAdapted::SaveRotationMatrix (const std::vector<int>& sites, const std::vector<Matrix>& m1, int state)
{
  Timer disktimer;
  int rank = mpigetrank();
  if (rank == 0)
    {

      char file [5000];
      int first = min(sites[0], *sites.rbegin()), last = max(sites[0], *sites.rbegin());
      if (state == -1)
	sprintf (file, "%s%s%d%s%d%s%d%s", dmrginp.save_prefix().c_str(), "/Rotation-", first, "-", last, ".", mpigetrank(),".state_average.tmp");
      else
	sprintf (file, "%s%s%d%s%d%s%d%s%d%s", dmrginp.save_prefix().c_str(), "/Rotation-", first, "-", last, ".", mpigetrank(),".state",state, ".tmp");
      p1out << "\t\t\t Saving Rotation Matrix :: " << file << endl;
      std::ofstream ofs(file, std::ios::binary);
      boost::archive::binary_oarchive save_mat(ofs);
      save_mat << m1;
      ofs.close();
    }
}

void SpinAdapted::LoadRotationMatrix (const std::vector<int>& sites, std::vector<Matrix>& m1, int state)
{
  Timer disktimer;
  int rank = mpigetrank();
  if (rank == 0)
  {
    char file [5000];
    int first = min(sites[0], *sites.rbegin()), last = max(sites[0], *sites.rbegin());
    if(state == -1)
      sprintf (file, "%s%s%d%s%d%s%d%s", dmrginp.save_prefix().c_str(), "/Rotation-", first, "-", last, ".", mpigetrank(),".state_average.tmp");
    else
      sprintf (file, "%s%s%d%s%d%s%d%s%d%s", dmrginp.save_prefix().c_str(), "/Rotation-", first, "-", last, ".", mpigetrank(),".state",state, ".tmp");
    p1out << "\t\t\t Loading Rotation Matrix :: " << file << endl;
    std::ifstream ifs(file, std::ios::binary);
    boost::archive::binary_iarchive load_mat(ifs);
    load_mat >> m1;
    ifs.close();
  }
}

void SpinAdapted::allocate(const StateInfo& row, const StateInfo& col, std::vector<Matrix>& rotations)
{
  rotations.resize(row.quanta.size());
  for (int i=0; i<row.quanta.size(); i++) {
    int nrows = row.quantaStates[i];
    int ncols = 0;
    for (int j=0; j<col.quanta.size(); j++) 
      if (col.quanta[j] == row.quanta[i])
	ncols += col.quantaStates[j];

    rotations[i].ReSize(nrows, ncols);
  }
}


bool SpinAdapted::can_connect(int n, int spin, int right_block_size)
{
  for(int alpha=0;alpha<=min(right_block_size/2,dmrginp.total_particle_number()-n);++alpha)
  {
    int beta = dmrginp.total_particle_number() - n - alpha;
    if(dmrginp.total_spin_number().getirrep() - (alpha - beta) == spin && beta <= right_block_size/2)
      return true;
  }
  return false;
}

int Binom(int n, int k)
{
  vector<int> b(n+1);
  b[0] = 1;
  for(int i=1;i<=n;++i)
  {
    b[i] = 1;
    for(int j=i-1;j>0;--j)
      if(INT_MAX - b[j-1] > 0)
        b[j] += b[j-1];
      else
        b[j] = INT_MAX;
  }
  return b[k];
}

int SpinAdapted::get_total_states(const int &this_size, const int &other_size)
{
  int maxj = this_size/2+1;
  vector<int> statesWithJ(maxj, 0);
  vector<int> temp(maxj, 0);

  statesWithJ[0] = 2; statesWithJ[1] = 1;
  for (int i=1; i<this_size/2; i++)
  {
    temp = statesWithJ;
    for (int j=0; j<maxj; j++)
      statesWithJ[j] =0;
    for (int j=0; j<maxj; j++)
    {
      if (INT_MAX - 2*temp[j] > statesWithJ[j])
	statesWithJ[j] += 2*temp[j];
      else statesWithJ[j] = INT_MAX;

      if (j+1 <maxj)
      {
	if (INT_MAX - temp[j] > statesWithJ[j+1])
	  statesWithJ[j+1] += temp[j];
	else
	  statesWithJ[j+1] = INT_MAX;
      }
      if (j-1 >= 0)
      {
	if (INT_MAX - temp[j] > statesWithJ[j-1])
	  statesWithJ[j-1] += temp[j];
	else
	  statesWithJ[j-1] = INT_MAX;
      }
    }
  }
  double retval = 0;
  for (int i=0; i<maxj; i++)
  {
    if (INT_MAX - retval > statesWithJ[i])
	retval += statesWithJ[i];
    else
      retval = INT_MAX;
  }
  return retval;
}

double SpinAdapted::assign_matrix_by_dm(std::vector<Matrix>& rotatematrix, std::vector<DiagonalMatrix>& eigenmatrix, SparseMatrix& transformmatrix, vector<pair<int, int> >& inorderwts, vector<vector<int> >& wtsbyquanta, int totalstatesbydm, int totalstatesbyquanta, int left_block_size, int right_block_size)
{
  const int min_states = totalstatesbydm;

  

  p2out << " \t\t\t assigning a total of " << min_states << " states using the dm alone " << endl;
  double totalnorm = 0.;
  rotatematrix.resize(eigenmatrix.size());

  for (int i = 0; i < totalstatesbydm; ++i)
    {
      int q = inorderwts[i].first;
      int qs = inorderwts[i].second;

      //if(i < min_states && eigenmatrix[q].element(qs, qs) > 1.e-12)
      if( eigenmatrix[q].element(qs, qs) > 1.e-13 || dmrginp.do_pdm())
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
      }
    }

  p2out << " \t\t\t assigning a total of " << totalstatesbyquanta << " states using quanta selection " << " for a norm of " << totalnorm << endl;

  int assignedbyq = 0;
  
  int nquanta = rotatematrix.size();
 
  int totalstatesleft = 0;
  for (int i = 0; i < nquanta; ++i)
    {
      totalstatesleft += wtsbyquanta[i].size();
    }

  p2out << " \t\t\t a total of " << totalstatesleft << " to be assigned " << endl;
  
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
    if( eigenmatrix[q].element(qs, qs) > 1.e-13 || dmrginp.do_pdm())
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
    }
  }
  
  double norm = 0.;
  for(int i=0;i<eigenmatrix.size();++i)
    for(int j=0;j<eigenmatrix[i].Nrows();++j)
      norm += eigenmatrix[i].element(j, j);
  p2out << " \t\t\t total norm: " << norm <<"  norm after truncation: "<<totalnorm<< endl;

  return norm-totalnorm;

  //return (1. - totalnorm/norm);
}

void SpinAdapted::diagonalise_dm(SparseMatrix& tracedMatrix, SparseMatrix& transformMatrix, std::vector<DiagonalMatrix>& eigenMatrix)
{
  int nquanta = tracedMatrix.nrows();
  eigenMatrix.resize(nquanta);
  vector<double> totalquantaweights(nquanta);
  for (int tQ = 0; tQ < nquanta; ++tQ)
    {
      int nStates = tracedMatrix.operator_element(tQ, tQ).Nrows ();
      DiagonalMatrix weights (nStates);
#ifdef USELAPACK
      diagonalise(tracedMatrix.operator_element(tQ, tQ), weights, transformMatrix.operator_element(tQ, tQ));
#else
      SymmetricMatrix dM (nStates);
      dM << tracedMatrix.operator_element(tQ,tQ);
      EigenValues (dM, weights, transformMatrix.operator_element(tQ,tQ));
#endif
      for(int i=0;i<weights.Nrows();++i)
        if(weights.element(i,i) < 1.e-14)
          weights.element(i,i) = 0.;
      eigenMatrix [tQ] = weights;
    }  
}

void SpinAdapted::svd_densitymat(SparseMatrix& tracedMatrix, SparseMatrix& transformMatrix, std::vector<DiagonalMatrix>& eigenMatrix) {
  // SVD of matrix M=(A,B,C)=USV^T
  // since MM^T=AA^T+BB^T+CC^T=USS^TU^T, we don't have to explicitly construct M
  int nquanta = tracedMatrix.nrows();
  eigenMatrix.resize(nquanta);
  vector<double> totalquantaweights(nquanta);
  for (int tQ = 0; tQ < nquanta; ++tQ) {
    int nStates = tracedMatrix.operator_element(tQ, tQ).Nrows ();
    DiagonalMatrix weights(nStates);
    Matrix M(nStates, nStates);
    M = 0.;
    for (int sQ = 0; sQ < nquanta; ++sQ)
      if (tracedMatrix.allowed(tQ, sQ))
        M += tracedMatrix.operator_element(tQ, tQ) * tracedMatrix.operator_element(tQ, tQ).t();

#ifdef USELAPACK
    diagonalise(M, weights, transformMatrix.operator_element(tQ, tQ));
#else
    SymmetricMatrix dM (nStates);
    dM << M;
    EigenValues(dM, weights, transformMatrix.operator_element(tQ,tQ));
#endif
    for (int i=0; i<weights.Nrows(); ++i) {
      if(weights.element(i,i) < 1.e-28)
        weights.element(i,i) = 0.;
      else
        weights.element(i,i) = sqrt(weights.element(i,i));
    }
    eigenMatrix[tQ] = weights;    
  }
}

void SpinAdapted::sort_weights(std::vector<DiagonalMatrix>& eigenMatrix, vector<pair<int, int> >& inorderwts, vector<vector<int> >& weightsbyquanta)
{
  // first sort weights with a multimap
  multimap<double, pair<int,int> > weightmap;
  int nquanta = eigenMatrix.size();
  vector<double> totalquantaweights(nquanta);

  for (int q = 0; q < nquanta; ++q)
  {
    //if (q == hfQuantaindex) continue; //the hartree fock quanta are always included first during warmup
    for (int qs = 0; qs < eigenMatrix[q].Nrows(); ++qs)
      {
	weightmap.insert (pair <double, pair<int,int> > (eigenMatrix[q].element(qs, qs), pair<int,int> (q, qs)));
	totalquantaweights[q] += eigenMatrix[q].element(qs, qs);

      }
  }

  multimap<double, pair<int,int> >::reverse_iterator w = weightmap.rbegin();


  // now put all the sorted indices in
  while (w != weightmap.rend())
    {
      inorderwts.push_back(make_pair(w->second.first, w->second.second));
      ++w;
    }

  // sort quantas by weight
  weightsbyquanta.resize(nquanta);
  for (int q = 0; q < nquanta; ++q)
    for (int qs = 0; qs < eigenMatrix[q].Nrows(); ++qs)
      weightsbyquanta[q].push_back(eigenMatrix[q].Nrows() - qs - 1);
}

