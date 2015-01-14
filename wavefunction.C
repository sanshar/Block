/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "wavefunction.h"
#include "spinblock.h"
#include "SpinQuantum.h"
#include "MatrixBLAS.h"
#include <boost/serialization/vector.hpp>
#include "pario.h"

void SpinAdapted::Wavefunction::initialise(const SpinQuantum dQ, const SpinBlock* b, const bool &onedot_)
{
  // initialized a ket wavefunction
  initialised = true;
  fermion = false;
  deltaQuantum.assign(1, dQ);
  onedot = onedot_; 
  resize(b->get_leftBlock()->get_stateInfo().quanta.size (), b->get_rightBlock()->get_stateInfo().quanta.size ());

  const SpinBlock* lBlock = b->get_leftBlock();
  const SpinBlock* rBlock = b->get_rightBlock();

  long totalmemory = 0;

  for (int lQ = 0; lQ < lBlock->get_stateInfo().quanta.size (); ++lQ)
    for (int rQ = 0; rQ < rBlock->get_stateInfo().quanta.size (); ++rQ)
    {
      for (int i = 0; i < deltaQuantum.size(); ++i)
        allowedQuantaMatrix(lQ, rQ) = deltaQuantum[i].allow(lBlock->get_stateInfo().quanta [lQ] , rBlock->get_stateInfo().quanta [rQ]);
      if (allowedQuantaMatrix(lQ, rQ))
        totalmemory += lBlock->get_stateInfo().quantaStates [lQ]* rBlock->get_stateInfo().quantaStates [rQ];
    }
  //double* largeArray = new double[totalmemory];
  long usedindex = 0;

  for (int lQ = 0; lQ < lBlock->get_stateInfo().quanta.size (); ++lQ)
    for (int rQ = 0; rQ < rBlock->get_stateInfo().quanta.size (); ++rQ)
    {
      if (allowedQuantaMatrix(lQ, rQ))
      {
	(*this)(lQ, rQ).ReSize (lBlock->get_stateInfo().quantaStates [lQ], rBlock->get_stateInfo().quantaStates [rQ]);//, &largeArray[usedindex]);
	SpinAdapted::Clear ((*this)(lQ, rQ));
	usedindex += lBlock->get_stateInfo().quantaStates [lQ]* rBlock->get_stateInfo().quantaStates [rQ];
      }
    }
}

void SpinAdapted::Wavefunction::initialisebra(const SpinQuantum dQ, const SpinBlock* b, const bool &onedot_)
{
  // initialized a bra wavefunction
  initialised = true;
  fermion = false;
  deltaQuantum.assign(1, dQ);
  onedot = onedot_; 
  resize(b->get_leftBlock()->get_braStateInfo().quanta.size (), b->get_rightBlock()->get_braStateInfo().quanta.size ());

  const SpinBlock* lBlock = b->get_leftBlock();
  const SpinBlock* rBlock = b->get_rightBlock();

  long totalmemory = 0;

  for (int lQ = 0; lQ < lBlock->get_braStateInfo().quanta.size (); ++lQ)
    for (int rQ = 0; rQ < rBlock->get_braStateInfo().quanta.size (); ++rQ)
    {
      for (int i = 0; i < deltaQuantum.size(); ++i)
        allowedQuantaMatrix(lQ, rQ) = deltaQuantum[i].allow(lBlock->get_braStateInfo().quanta [lQ] , rBlock->get_braStateInfo().quanta [rQ]);
      if (allowedQuantaMatrix(lQ, rQ))
        totalmemory += lBlock->get_braStateInfo().quantaStates [lQ]* rBlock->get_braStateInfo().quantaStates [rQ];
    }
  //double* largeArray = new double[totalmemory];
  long usedindex = 0;

  for (int lQ = 0; lQ < lBlock->get_braStateInfo().quanta.size (); ++lQ)
    for (int rQ = 0; rQ < rBlock->get_braStateInfo().quanta.size (); ++rQ)
    {
      if (allowedQuantaMatrix(lQ, rQ))
      {
	(*this)(lQ, rQ).ReSize (lBlock->get_braStateInfo().quantaStates [lQ], rBlock->get_braStateInfo().quantaStates [rQ]);//, &largeArray[usedindex]);
	SpinAdapted::Clear ((*this)(lQ, rQ));
	usedindex += lBlock->get_braStateInfo().quantaStates [lQ]* rBlock->get_braStateInfo().quantaStates [rQ];
      }
    }
}

void SpinAdapted::Wavefunction::initialise(const vector<SpinQuantum>& dQ, const SpinBlock* b, const bool &onedot_)
{
  initialised = true;
  fermion = false;
  deltaQuantum = dQ;
  onedot = onedot_; 
  resize(b->get_leftBlock()->get_stateInfo().quanta.size (), b->get_rightBlock()->get_stateInfo().quanta.size ());

  const SpinBlock* lBlock = b->get_leftBlock();
  const SpinBlock* rBlock = b->get_rightBlock();

  long totalmemory = 0;
  for (int lQ = 0; lQ < lBlock->get_stateInfo().quanta.size (); ++lQ)
    for (int rQ = 0; rQ < rBlock->get_stateInfo().quanta.size (); ++rQ) {
      allowedQuantaMatrix(lQ, rQ) = false;
      for (int i = 0; i < deltaQuantum.size(); ++i)
        if (deltaQuantum[i].allow(lBlock->get_stateInfo().quanta [lQ] , rBlock->get_stateInfo().quanta [rQ])) {
          allowedQuantaMatrix(lQ, rQ) = true;
          break;
        }
      if (allowedQuantaMatrix(lQ, rQ))
	    totalmemory += lBlock->get_stateInfo().quantaStates [lQ]* rBlock->get_stateInfo().quantaStates [rQ];
    }
  //double* largeArray = new double[totalmemory];
  long usedindex = 0;

  for (int lQ = 0; lQ < lBlock->get_stateInfo().quanta.size (); ++lQ)
    for (int rQ = 0; rQ < rBlock->get_stateInfo().quanta.size (); ++rQ)
    {
      if (allowedQuantaMatrix(lQ, rQ))
      {
	(*this)(lQ, rQ).ReSize (lBlock->get_stateInfo().quantaStates [lQ], rBlock->get_stateInfo().quantaStates [rQ]);//, &largeArray[usedindex]);
	SpinAdapted::Clear ((*this)(lQ, rQ));
	usedindex += lBlock->get_stateInfo().quantaStates [lQ]* rBlock->get_stateInfo().quantaStates [rQ];
      }
    }
}

void SpinAdapted::Wavefunction::initialisebra(const vector<SpinQuantum>& dQ, const SpinBlock* b, const bool &onedot_)
{
  initialised = true;
  fermion = false;
  deltaQuantum = dQ;
  onedot = onedot_; 
  resize(b->get_leftBlock()->get_braStateInfo().quanta.size (), b->get_rightBlock()->get_braStateInfo().quanta.size ());

  const SpinBlock* lBlock = b->get_leftBlock();
  const SpinBlock* rBlock = b->get_rightBlock();

  long totalmemory = 0;
  for (int lQ = 0; lQ < lBlock->get_braStateInfo().quanta.size (); ++lQ)
    for (int rQ = 0; rQ < rBlock->get_braStateInfo().quanta.size (); ++rQ) {
      allowedQuantaMatrix(lQ, rQ) = false;
      for (int i = 0; i < deltaQuantum.size(); ++i)
        if (deltaQuantum[i].allow(lBlock->get_braStateInfo().quanta [lQ] , rBlock->get_braStateInfo().quanta [rQ])) {
          allowedQuantaMatrix(lQ, rQ) = true;
          break;
        }
      if (allowedQuantaMatrix(lQ, rQ))
	    totalmemory += lBlock->get_braStateInfo().quantaStates [lQ]* rBlock->get_braStateInfo().quantaStates [rQ];
    }
  //double* largeArray = new double[totalmemory];
  long usedindex = 0;

  for (int lQ = 0; lQ < lBlock->get_braStateInfo().quanta.size (); ++lQ)
    for (int rQ = 0; rQ < rBlock->get_braStateInfo().quanta.size (); ++rQ)
    {
      if (allowedQuantaMatrix(lQ, rQ))
      {
	(*this)(lQ, rQ).ReSize (lBlock->get_braStateInfo().quantaStates [lQ], rBlock->get_braStateInfo().quantaStates [rQ]);//, &largeArray[usedindex]);
	SpinAdapted::Clear ((*this)(lQ, rQ));
	usedindex += lBlock->get_braStateInfo().quantaStates [lQ]* rBlock->get_braStateInfo().quantaStates [rQ];
      }
    }
}

void SpinAdapted::Wavefunction::FlattenInto (Matrix& C)
{
  int flatIndex = 0;
  for (int lQ = 0; lQ < nrows (); ++lQ)
    for (int rQ = 0; rQ < ncols (); ++rQ)
      if (allowed(lQ, rQ))
	flatIndex += operator_element(lQ,rQ).Nrows()*operator_element(lQ,rQ).Ncols();
  
  C.ReSize(flatIndex,1);

  flatIndex = 0;
  for (int lQ = 0; lQ < nrows (); ++lQ)
    for (int rQ = 0; rQ < ncols (); ++rQ)
      if (allowed(lQ, rQ))
        for (int lQState = 0; lQState < operator_element(lQ, rQ).Nrows (); ++lQState)
          for (int rQState = 0; rQState < operator_element(lQ, rQ).Ncols (); ++rQState)
	    {
	      C.element (flatIndex,0) = operator_element(lQ, rQ).element (lQState, rQState);
	      ++flatIndex;
	    }
}

void SpinAdapted::Wavefunction::CollectFrom (const RowVector& C)
{
  int flatIndex = 0;
  for (int lQ = 0; lQ < nrows (); ++lQ)
    for (int rQ = 0; rQ < ncols (); ++rQ)
      if (allowedQuantaMatrix (lQ, rQ))
        for (int lQState = 0; lQState < operator_element(lQ, rQ).Nrows (); ++lQState)
          for (int rQState = 0; rQState < operator_element(lQ, rQ).Ncols (); ++rQState)
	    {
	      operator_element(lQ, rQ).element (lQState, rQState) = C.element (flatIndex);
	      ++flatIndex;
	    }
}


  
void SpinAdapted::Wavefunction::SaveWavefunctionInfo (const StateInfo &waveInfo, const std::vector<int>& sites, const int wave_num)
{
  char file [5000];
  int first = min(sites[0], *sites.rbegin()), last = max(sites[0], *sites.rbegin());
  sprintf (file, "%s%s%d%s%d%s%d%s%d%s", dmrginp.save_prefix().c_str(), "/wave-", first, "-", last, ".", mpigetrank(), ".", wave_num, ".tmp");
  p1out << "\t\t\t Saving Wavefunction " << file << endl;
  if (mpigetrank() == 0)
    {
      std::ofstream ofs(file, std::ios::binary);
      boost::archive::binary_oarchive save_wave(ofs);
      save_wave << onedot << waveInfo << *waveInfo.leftStateInfo << *(waveInfo.leftStateInfo->leftStateInfo);
      save_wave << *(waveInfo.leftStateInfo->rightStateInfo) << *waveInfo.rightStateInfo;
      if(!onedot)
	save_wave << *(waveInfo.rightStateInfo->leftStateInfo) << *(waveInfo.rightStateInfo->rightStateInfo);
      this->Save (ofs);
      ofs.close();
    }

}

void SpinAdapted::Wavefunction::LoadWavefunctionInfo (StateInfo &waveInfo, const std::vector<int>& sites, const int wave_num)
{
  char file [5000];
  int first = min(sites[0], *sites.rbegin()), last = max(sites[0], *sites.rbegin());
  sprintf (file, "%s%s%d%s%d%s%d%s%d%s", dmrginp.load_prefix().c_str(), "/wave-", first, "-", last, ".", mpigetrank(), ".", wave_num, ".tmp");
  p1out << "\t\t\t Loading Wavefunction " << file << endl;
  waveInfo.Allocate ();
  if (mpigetrank() == 0)
    {
      std::ifstream ifs(file, std::ios::binary);
      boost::archive::binary_iarchive load_wave(ifs);
      load_wave >> onedot >> waveInfo >> *waveInfo.leftStateInfo >> *(waveInfo.leftStateInfo->leftStateInfo)
		>> *(waveInfo.leftStateInfo->rightStateInfo) >> *waveInfo.rightStateInfo;
      if(!onedot)
	load_wave >> *(waveInfo.rightStateInfo->leftStateInfo) >> *(waveInfo.rightStateInfo->rightStateInfo);
      this->Load (ifs);
      ifs.close();
    }
}

void SpinAdapted::Wavefunction::CollectQuantaAlongRows (const StateInfo& sRow, const StateInfo& sCol)
{
  //mdebugcheck("before collectquantaalongrows");
  try
    {
      Timer ctimer;
      StateInfo tmpState = sRow;
      tmpState.CollectQuanta ();
      Wavefunction tmpOper;

      tmpOper.AllowQuantaFor(tmpState, sCol, deltaQuantum);
      ObjectMatrix<Matrix*> matRef;
      for (int i = 0; i < tmpState.quanta.size (); ++i)
	    for (int j = 0; j < sCol.quanta.size (); ++j) {
	      std::vector<int> dum (1); dum [0] = j;
	      if (tmpOper.allowed(i, j)) {
		    OperatorMatrixReference (matRef, tmpState.oldToNewState [i], dum);
		    CatenateProduct (matRef, tmpOper.operator_element(i,j), true);
	      }
	    }
      allowedQuantaMatrix = tmpOper.allowedQuantaMatrix;
      copy(tmpOper.operatorMatrix, operatorMatrix);
    }
  catch (Exception)
    {
      Exception::what ();
      abort ();
    }
  //mdebugcheck("after collectquantaalongrows");
}



void SpinAdapted::Wavefunction::UnCollectQuantaAlongRows (const StateInfo& sRow, const StateInfo& sCol)
{
  try
    {
      Wavefunction tmpOper;
      tmpOper.AllowQuantaFor (*sRow.unCollectedStateInfo, sCol, deltaQuantum);
      for (int i = 0; i < sRow.quanta.size (); ++i)
	{
	  const std::vector<int>& oldToNewStateI = sRow.oldToNewState [i];
	  int firstRow = 0;
	  for (int iSub = 0; iSub < oldToNewStateI.size (); ++iSub)
	    {
	      int unCollectedI = oldToNewStateI [iSub];
	      int lastRowSize = sRow.unCollectedStateInfo->quantaStates [unCollectedI];
	      for (int j = 0; j < sCol.quanta.size (); ++j)
		if (tmpOper.allowedQuantaMatrix (unCollectedI, j))
		  tmpOper.operator_element(unCollectedI, j) = operator_element(i, j).Rows (firstRow + 1, firstRow + lastRowSize);
	      firstRow += lastRowSize;
	    }
	}
      allowedQuantaMatrix = tmpOper.allowedQuantaMatrix;
      copy(tmpOper.operatorMatrix,operatorMatrix);

      //tmpOper.Release ();                                                                                                                 
      //*this = tmpOper;                                                                                                                    
    }
  catch (Exception)
    {
      Exception::what ();
      abort ();
    }
}



void SpinAdapted::Wavefunction::CollectQuantaAlongColumns (const StateInfo& sRow, const StateInfo& sCol)
{
  //mdebugcheck("before collectquantaalongcolumns");
  try
    {
      StateInfo tmpState = sCol;
      tmpState.CollectQuanta ();
      Wavefunction tmpOper;
      tmpOper.AllowQuantaFor (sRow, tmpState, deltaQuantum);
      ObjectMatrix<Matrix*> matRef;
      for (int i = 0; i < sRow.quanta.size (); ++i)
	for (int j = 0; j < tmpState.quanta.size (); ++j)
	  {
	    std::vector<int> dum (1); dum [0] = i;
	    if (tmpOper.allowed(i, j))
	      {
		OperatorMatrixReference (matRef, dum, tmpState.oldToNewState [j]);
		CatenateProduct (matRef, tmpOper.operator_element(i,j));
	      }
	  }
      allowedQuantaMatrix = tmpOper.allowedQuantaMatrix;
      copy(tmpOper.operatorMatrix,operatorMatrix);
    }
  catch (Exception)
    {
      Exception::what ();
      abort ();
    }
  //mdebugcheck("after collectquantaalongcolumns");
}


void  SpinAdapted::Wavefunction::UnCollectQuantaAlongColumns (const StateInfo& sRow, const StateInfo& sCol)
{
  //try                                                                                                                                   
  // {                                                                                                                                    
  Wavefunction tmpOper;
  tmpOper.AllowQuantaFor (sRow, *sCol.unCollectedStateInfo, deltaQuantum);
  for (int i = 0; i < sCol.quanta.size (); ++i)
    {
      const std::vector<int>& oldToNewStateI = sCol.oldToNewState [i];
      int firstCol = 0;
      for (int iSub = 0; iSub < oldToNewStateI.size (); ++iSub)
	{
	  int unCollectedI = oldToNewStateI [iSub];
	  int lastColSize = sCol.unCollectedStateInfo->quantaStates [unCollectedI];
	  for (int j = 0; j < sRow.quanta.size (); ++j)
	    if (tmpOper.allowed(j, unCollectedI)){
	      tmpOper.operatorMatrix (j, unCollectedI) = operatorMatrix (j, i).Columns (firstCol + 1, firstCol + lastColSize);
	    }
	  firstCol += lastColSize;
	}
    }
  allowedQuantaMatrix = tmpOper.allowedQuantaMatrix;
  copy(tmpOper.operatorMatrix,operatorMatrix);
  /*}                                                                                                                                   
  catch (Exception)                                                                                                                       
  {                                                                                                                                       
    Exception::what ();                                                                                                                   
    abort ();                                                                                                                             
    }*/
}

/*
void SpinAdapted::Wavefunction::AllowQuantaFor (const StateInfo& sRow, const StateInfo& sCol, const SpinQuantum q)
{
  deltaQuantum.assign(1, q);
  resize(sRow.quanta.size (), sCol.quanta.size ());
  initialised = true;
  for (int i = 0; i < sRow.quanta.size (); ++i)
    for (int j = 0; j < sCol.quanta.size (); ++j) {
	  allowedQuantaMatrix(i, j) = q.allow(sRow.quanta[i], sCol.quanta[j]);
	  operator_element(i, j).CleanUp ();
	  if (allowedQuantaMatrix(i,j))
	  {
	    operator_element(i, j).ReSize (sRow.quantaStates [i], sCol.quantaStates [j]);
	    SpinAdapted::Clear (operator_element(i, j));
	  }
    }
}
*/

void SpinAdapted::Wavefunction::AllowQuantaFor (const StateInfo& sRow, const StateInfo& sCol, const vector<SpinQuantum>& q)
{
  deltaQuantum = q;
  resize(sRow.quanta.size (), sCol.quanta.size ());
  initialised = true;
  for (int i = 0; i < sRow.quanta.size (); ++i)
    for (int j = 0; j < sCol.quanta.size (); ++j) {
      allowedQuantaMatrix(i, j) = false;
      for (int k = 0; k < q.size(); ++k) {
        if (q[k].allow(sRow.quanta[i], sCol.quanta[j])) {
	      allowedQuantaMatrix(i, j) = true;
          break;
        }
      }
	  operator_element(i, j).CleanUp ();
	  if (allowedQuantaMatrix(i,j))
	  {
	    operator_element(i, j).ReSize (sRow.quantaStates [i], sCol.quantaStates [j]);
	    SpinAdapted::Clear (operator_element(i, j));
	  }
    }
}

SpinAdapted::Wavefunction& SpinAdapted::Wavefunction::operator+=(const Wavefunction& other)
{
  for (int i = 0; i < nrows(); ++i)
    for (int j = 0; j < ncols(); ++j)
      if (allowed(i, j))
        {
          assert(other.allowed(i, j));
          MatrixScaleAdd(1., other.operator_element(i, j), operator_element(i, j));
        }
  return *this;
}
