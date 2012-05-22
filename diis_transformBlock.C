#include "diis.h"
#include "global.h"
#include "solver.h"
#include "initblocks.h"
#include "MatrixBLAS.h"
#include "density.h"
#include "rotationmat.h"
#include "davidson.h"
#include "linear.h"
#include "guess_wavefunction.h"
#include <include/sortutils.h>

#ifndef SERIAL
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif

using namespace boost;
using namespace std;



void SpinAdapted::DIIS::transformBlock(SpinBlock& big, const std::vector<Wavefunction>& wave_solutions, SweepParams &sweepParams) {

  int sweepiter = sweepParams.get_sweep_iter();
  SpinBlock& newSystem =* big.get_leftBlock();
  vector<Matrix> rotateMatrix;
  StateInfo stateInfo = newSystem.get_stateInfo();
  DensityMatrix tracedMatrix;
  tracedMatrix.allocate(stateInfo);
  double twodotnoise = 0.0;
  
  tracedMatrix.makedensitymatrix(wave_solutions, big, dmrginp.weights(sweepiter), 0.0e0, 0.0e0, false);
  std::vector<Matrix> prevRotation;

  LoadRotationMatrix (newSystem.get_sites(), prevRotation, 0); //check if the rotation block of state 0 should be read?

  if (!mpigetrank())
    {
      // find and sort weight info
      DensityMatrix transformmatrix;
      transformmatrix.allocate(stateInfo);
      std::vector<DiagonalMatrix> eigenMatrix;
      diagonalise_dm(tracedMatrix, transformmatrix, eigenMatrix);
      vector<pair<int, int> > inorderwts;
      vector<vector<int> > wtsbyquanta;

      sort_weights(eigenMatrix, inorderwts, wtsbyquanta);


      // make transformation matrix by various algorithms
      int totalstatesbydm = min(static_cast<int>(inorderwts.size()), sweepParams.get_keep_states());
      int totalstatesbyquanta = min(static_cast<int>(inorderwts.size()), 0) - totalstatesbydm;
      if (totalstatesbyquanta < 0) totalstatesbyquanta = 0;


      std::vector<int> statesperquanta(wave_solutions[0].nrows(),0);
      for (int i=0; i<wave_solutions[0].nrows(); i++)
	for (int j=0; j<wave_solutions[0].ncols(); j++)
	  if (wave_solutions[0].allowed(i,j))
	    statesperquanta.at(i) = wave_solutions[0](i,j).Ncols();

      double error = assign_dm(rotateMatrix, eigenMatrix, transformmatrix, inorderwts, wtsbyquanta, totalstatesbydm, 
				  totalstatesbyquanta, statesperquanta);
      pout << "\t\t\t Total discarded weight "<<error<<endl;

      /*
      cout << "weights in density matrix"<<endl;
      for (int i=0; i<eigenMatrix.size(); i++)
	for (int j=0; j<eigenMatrix.at(i).Ncols(); j++)
	  if (eigenMatrix.at(i)(j+1) >= 1.0e-13)
	    cout << i<<"  "<<j<<"  "<<eigenMatrix.at(i)(j+1)<<endl;
      */
    }

#ifndef SERIAL
  mpi::communicator world;
  broadcast(world, rotateMatrix, 0);
#endif
  maxOverlapPrevRotation(rotateMatrix, newSystem.get_sites(), prevRotation);

  for (int i=0; i<dmrginp.nroots(); i++)
    SaveRotationMatrix (newSystem.get_sites(), rotateMatrix, i);
  newSystem.transform_operators(rotateMatrix);
  mcheck("after rotation and transformation of block");
  pout <<newSystem<<endl;


  //mcheck("After renorm transform");
}

void SpinAdapted::DIIS::maxOverlapPrevRotation(std::vector<Matrix>& rotation, const std::vector<int>& sites, std::vector<Matrix>& prevRotation)
{
  if (buildup && currentIndex == 0)
    return;
  
  assert(rotation.size() == prevRotation.size());
  for (int i=0; i<rotation.size(); i++)
  {
    if(rotation.at(i).Nrows() == 0) continue;

    int mrow, mcol;

    Matrix m(rotation.at(i).Nrows(),rotation.at(i).Ncols());  
    copy(rotation.at(i), m);

    if (m.Ncols() != prevRotation.at(i).Ncols() || m.Nrows() != prevRotation.at(i).Nrows()) {
      cout <<m.Ncols() <<"  "<< prevRotation.at(i).Ncols() <<"  "<< m.Nrows() <<"  "<< prevRotation.at(i).Nrows()<<endl;
      cout << m<<endl;
      cout << prevRotation.at(i)<<endl;

      exit(0);
    }

    Matrix C(m.Ncols(), m.Ncols());
    MatrixMultiply(m, 't', prevRotation.at(i), 'n', C, 1.0, 0.0);

    //perform svd on C
    Matrix U, V;
    DiagonalMatrix d;
    svd(C, d, U, V);

    //the orthogonal matrix
    Matrix Q(m.Ncols(), m.Ncols());
    MatrixMultiply(U, 'n', V, 't', Q, 1.0, 0.0);

    //now rotate the matrix m
    MatrixMultiply(m, 'n', Q, 'n', rotation.at(i), 1.0, 0.0);
    
    copy(rotation.at(i), m);


    std::vector<int> reorder(m.Ncols());
    std::vector<int> sign(m.Ncols(), 1);
    std::vector<double> maxo(m.Ncols(),0);
    Matrix overlapMat(m.Ncols(), prevRotation.at(i).Ncols());
    for (int j=0; j<min(m.Ncols(), prevRotation.at(i).Ncols()); j++)
    {
      double maxoverlap = 0.0;
      double overlap = DDOT(m.Nrows(), &(m(1, j+1)), m.Ncols(), &(prevRotation.at(i)(1, j+1)), m.Ncols());
      maxoverlap = fabs(overlap);
      maxo.at(j) = overlap;
      reorder.at(j) = j;
      sign.at(j) = overlap > 0.0 ? 1 : -1;           
      
    }


    for (int j=0; j<m.Ncols(); j++)
      for (int k=0; k<prevRotation.at(i).Ncols(); k++)
	overlapMat(j+1, k+1) = DDOT(m.Nrows(), &(m(1, j+1)), m.Ncols(), &(prevRotation.at(i)(1, k+1)), m.Ncols());
    for (int k=0; k<overlapMat.Ncols(); k++) {
      if (fabs(overlapMat(k+1, k+1)) <0.8) {
	cout << "matrix "<<i<<endl;
	cout << "rows is present for matrix "<<i<<endl;
	cout << overlapMat<<endl;
      }
    }

  }
}



double SpinAdapted::DIIS::assign_dm(std::vector<Matrix>& rotatematrix, std::vector<DiagonalMatrix>& eigenmatrix, SparseMatrix& transformmatrix, vector<pair<int, int> >& inorderwts, vector<vector<int> >& wtsbyquanta, int totalstatesbydm, int totalstatesbyquanta, std::vector<int>& statesperquanta)
{
  /*const int left_states = get_total_states(left_block_size*2, right_block_size*2);
  pout << " \t\t\t possible left states " << left_states << endl;
  const int right_states = get_total_states(right_block_size*2, left_block_size*2);
  pout << " \t\t\t possible right states " << right_states << endl;
  const int min_states = min(totalstatesbydm, min(left_states, right_states));
  */
  const int min_states = totalstatesbydm;

  

  pout << " \t\t\t assigning a total of " << min_states << " states using the dm alone " << endl;
  double totalnorm = 0.;
  rotatematrix.resize(eigenmatrix.size());

  for (int i = 0; i < inorderwts.size(); ++i)
    {
      int q = inorderwts[i].first;
      int qs = inorderwts[i].second;

      if (rotatematrix[q].Ncols() == statesperquanta.at(q))
	continue;
      if( eigenmatrix[q].element(qs, qs) > 1.e-13 || rotatematrix[q].Ncols() < statesperquanta.at(q))
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

  pout << " \t\t\t assigning a total of " << totalstatesbyquanta << " states using quanta selection " << " for a norm of " << totalnorm << endl;

  int assignedbyq = 0;
  
  int nquanta = rotatematrix.size();
 
  int totalstatesleft = 0;
  for (int i = 0; i < nquanta; ++i)
    {
      totalstatesleft += wtsbyquanta[i].size();
    }

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
    if( eigenmatrix[q].element(qs, qs) > 1.e-13)
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
  pout << " \t\t\t total norm: " << norm <<"  norm after truncation: "<<totalnorm<< endl;

  return norm-totalnorm;

  //return (1. - totalnorm/norm);
}
