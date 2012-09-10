/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is free software: you can redistribute it and/or modify         
it under the terms of the GNU General Public License as published by         
the Free Software Foundation, either version 3 of the License, or            
(at your option) any later version.                                          
                                                                             
This program is distributed in the hope that it will be useful,              
but WITHOUT ANY WARRANTY; without even the implied warranty of               
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
GNU General Public License for more details.                                 
                                                                             
You should have received a copy of the GNU General Public License            
along with this program.  If not, see <http://www.gnu.org/licenses/>.        

This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "BaseOperator.h"
#include "MatrixBLAS.h"
#include "spinblock.h"
#include "distribute.h"
#include "tensor_operator.h"
#include "blas_calls.h"

namespace SpinAdapted{

double getCommuteParity(SpinQuantum a, SpinQuantum b, SpinQuantum c)
{
  int aspin = a.get_s(), airrep = a.get_symm().getirrep();
  int bspin = b.get_s(), birrep = b.get_symm().getirrep();
  int cspin = c.get_s(), cirrep = c.get_symm().getirrep();

  int an = a.get_n(), bn = b.get_n();
  int parity = IsFermion(a) && IsFermion(b) ? -1 : 1;
  for (int asz = -aspin; asz<aspin+1; asz+=2)
  for (int bsz = -bspin; bsz<bspin+1; bsz+=2)
  for (int al = 0; al<Symmetry::sizeofIrrep(airrep); al++)
  for (int bl = 0; bl<Symmetry::sizeofIrrep(birrep); bl++)
  {
    double cleb = cleb_(aspin, asz, bspin, bsz, cspin, cspin);
    double clebdinfh = Symmetry::spatial_cg(airrep, birrep, cirrep, al, bl, 0);
    if (fabs(cleb) <= 1.0e-14 || fabs(clebdinfh) <= 1.0e-14)
      continue;
    else
      return parity*cleb*clebdinfh/cleb_(bspin, bsz, aspin, asz, cspin, cspin)/Symmetry::spatial_cg(birrep, airrep, cirrep, bl, al, 0);
  }
  cout << "Major trouble, getCommuteParity asked for three inappropriate operators"<<endl;
  cout << a<<"  "<<b<<"  "<<c<<endl;
  return 1.0;
}



void SparseMatrix::allocate(const SpinBlock& b)
{
  allocate(b.get_stateInfo());
}

void SparseMatrix::allocate(const StateInfo& s)
{
  resize(s.quanta.size(), s.quanta.size());
  long totalmemory = 0;
  for (int i = 0; i < allowedQuantaMatrix.Nrows (); ++i)
    for (int j = 0; j < allowedQuantaMatrix.Ncols (); ++j)
    {
      allowedQuantaMatrix (i,j) = s.quanta[i].allow(deltaQuantum, s.quanta[j]);
      if (allowedQuantaMatrix (i,j)) 
	{
	  operatorMatrix (i,j).ReSize (s.quantaStates.at(i), s.quantaStates.at(j));//, largeArray+usedindex);
	  SpinAdapted::Clear (operatorMatrix (i,j));
	}

    }     
}

void SparseMatrix::allocate(const StateInfo& sr, const StateInfo& sc)
{
  resize(sr.quanta.size(), sc.quanta.size());
  long totalmemory = 0;
  for (int i = 0; i < allowedQuantaMatrix.Nrows (); ++i)
    for (int j = 0; j < allowedQuantaMatrix.Ncols (); ++j)
    {
      allowedQuantaMatrix (i,j) = sr.quanta[i].allow(deltaQuantum, sc.quanta[j]);
      if (allowedQuantaMatrix (i,j)) 
	{
	  operatorMatrix (i,j).ReSize (sr.quantaStates.at(i), sc.quantaStates.at(j));//, largeArray+usedindex);
	  SpinAdapted::Clear (operatorMatrix (i,j));
	}

    }     
}

void SparseMatrix::CleanUp ()
{
  built = false;
  initialised = false;
  fermion = 0;
  deltaQuantum = SpinQuantum (0, 0, IrrepSpace(0));
  orbs.resize(0);
  allowedQuantaMatrix.ReSize (0,0);
  operatorMatrix.ReSize (0,0);
}

const Transposeview Transpose(SparseMatrix& op) { return Transposeview(op); };

ostream& operator<< (ostream& os, const SparseMatrix& a)
{
  assert (a.initialised);
	for (int i = 0; i < a.allowedQuantaMatrix.Nrows (); ++i)
	for (int j = 0; j < a.allowedQuantaMatrix.Ncols (); ++j)
	{
	  if (a.allowed(i, j)) 
	    os << i << " " << j << endl << a.operator_element(i, j) << endl;

	}
      return os;
}

double SparseMatrix::memoryUsed(const SpinBlock& b)
{
  StateInfo stateinfo = b.get_stateInfo();
  double memory = 0.0;

  for (int i=0; i < stateinfo.quanta.size(); i++)
    for (int j=0; j<stateinfo.quanta.size(); j++)
      if (allowedQuantaMatrix(i,j) ) {
	memory += 8.0*operatorMatrix(i,j).Storage();
      }
  return memory;
}
  
void SparseMatrix::buildUsingCsf(const SpinBlock& b, vector< vector<Csf> >& ladders, std::vector< Csf >& s) 
{
  StateInfo stateinfo = b.get_stateInfo();
  built = true;
  allocate(stateinfo);

  for (int i=0; i < stateinfo.quanta.size(); i++)
    for (int j=0; j<stateinfo.quanta.size(); j++)
      if (allowedQuantaMatrix(i,j) )
	for (int jq =stateinfo.unBlockedIndex[j]; jq < stateinfo.unBlockedIndex[j]+stateinfo.quantaStates[j]; jq++) 
	{
	  for (int iq =stateinfo.unBlockedIndex[i]; iq < stateinfo.unBlockedIndex[i]+stateinfo.quantaStates[i]; iq++) 
	    operatorMatrix(i,j)(iq-stateinfo.unBlockedIndex[i]+1, jq-stateinfo.unBlockedIndex[j]+1) = redMatrixElement(s[iq], ladders[jq], &b);
	}

}

void SparseMatrix::Randomise ()
{
  for (int lQ = 0; lQ < nrows(); ++lQ)
    for (int rQ = 0; rQ < ncols(); ++rQ)
      if (allowed(lQ, rQ))
	SpinAdapted::Randomise (operator_element(lQ, rQ));
}

double DotProduct(const SparseMatrix& lhs, const SparseMatrix& rhs)
{
  double result = 0.;
  for (int lQ = 0; lQ < lhs.nrows(); ++lQ)
    for (int rQ = 0; rQ < lhs.ncols (); ++rQ)
      if (lhs.allowed(lQ, rQ) && rhs.allowed(lQ, rQ))
	result += MatrixDotProduct(lhs.operator_element(lQ, rQ), rhs.operator_element(lQ, rQ));

  return result;
}

void Scale(double d, SparseMatrix& a)
{
  for (int lQ = 0; lQ < a.nrows(); ++lQ)
    for (int rQ = 0; rQ < a.ncols(); ++rQ)
      if (a.allowed(lQ, rQ))
        MatrixScale(d, a.operator_element(lQ, rQ));
}

void ScaleAdd(double d, const SparseMatrix& a, SparseMatrix& b)
{
  for (int lQ = 0; lQ < a.nrows(); ++lQ)
    for (int rQ = 0; rQ < a.ncols(); ++rQ)
      if (a.allowed(lQ, rQ))
        {
	  if (!b.allowed(lQ, rQ))
	    cout <<"Not a valid addition"<<endl;
          assert(b.allowed(lQ, rQ));
          MatrixScaleAdd(d, a.operator_element(lQ, rQ), b.operator_element(lQ, rQ));
        }
}

void Normalise(SparseMatrix& a, int* success)
{
  a.Normalise(success);
}

void SparseMatrix::Normalise (int* success)
{
  double normalisation = DotProduct(*this, *this);
  if(normalisation > 1.e-12)
    Scale(1./sqrt(normalisation), *this);
  else {
    pout << "\t\t\t Warning :: Didn't Normalise, because norm is too small" << endl;
    *success = 1; //not successful in normlaising                                                                                  
  }
}

void SparseMatrix::Clear ()
{
  built = false;
  for (int i = 0; i < allowedQuantaMatrix.Nrows (); ++i)
    for (int j = 0; j < allowedQuantaMatrix.Ncols (); ++j)
      if (allowedQuantaMatrix (i,j)) SpinAdapted::Clear (operatorMatrix (i,j));
}

void assignloopblock(SpinBlock*& loopblock, SpinBlock*& otherblock, SpinBlock* leftBlock,
			    SpinBlock* rightBlock)
{
  loopblock = leftBlock;
  otherblock = rightBlock;
  if (!leftBlock->is_loopblock()) {loopblock = rightBlock; otherblock = leftBlock;}
}

void copy(const ObjectMatrix<Matrix>& a, ObjectMatrix<Matrix>& b)
{
  b.resize(a.Nrows(), a.Ncols());
  for (int i = 0; i < a.Nrows(); ++i)
    for (int j = 0; j < a.Ncols(); ++j)
      copy(a(i, j), b(i, j));
}

void copy(const Matrix& a, Matrix& b)
{
  if ((b.Nrows() != a.Nrows()) || (b.Ncols() != a.Ncols()))
    b.ReSize(a.Nrows(), a.Ncols());

#ifdef BLAS
  DCOPY((FORTINT) a.Storage(), a.Store(), (FORTINT) 1, b.Store(), (FORTINT) 1);
#else
  b = a;
#endif
}

void SparseMatrix::OperatorMatrixReference (ObjectMatrix<Matrix*>& m, const std::vector<int>& oldToNewStateI, 
					const std::vector<int>& oldToNewStateJ)
{
  int rows = oldToNewStateI.size ();
  int cols = oldToNewStateJ.size ();
  m.ReSize (rows, cols);
  for (int i = 0; i < rows; ++i)
    for (int j = 0; j < cols; ++j)
      {
	assert (allowedQuantaMatrix (oldToNewStateI [i], oldToNewStateJ [j]));
	m (i,j) = &operatorMatrix (oldToNewStateI [i], oldToNewStateJ [j]);
      }
}

//Renormalization functions for core and virtual operators                                                                                
void SparseMatrix::renormalise_transform(const std::vector<Matrix>& rotate_matrix, const StateInfo *stateinfo)
{
  ObjectMatrix<Matrix> tmp = operatorMatrix; //cannot instantiate a SparseMatrix and so instantiating a Cre

  this->allocate(*stateinfo); // new allocations                                                                                           

  int newQ = 0;
  for (int Q = 0; Q < rotate_matrix.size (); ++Q)
    if (rotate_matrix[Q].Ncols () != 0)
      {
	int newQPrime = 0;
	for (int QPrime = 0; QPrime < rotate_matrix.size (); ++QPrime)
	  if (rotate_matrix[QPrime].Ncols () != 0)
	    {
	      if (this->allowedQuantaMatrix (newQ, newQPrime))
		MatrixRotate (rotate_matrix[Q], tmp(Q, QPrime), rotate_matrix[QPrime],
			      this->operatorMatrix (newQ, newQPrime) );
	      ++newQPrime;
	    }
	++newQ;
      }

}

void SparseMatrix::build_and_renormalise_transform(SpinBlock *big, const opTypes &ot, const std::vector<Matrix>& rotate_matrix, 
					       const StateInfo *newStateInfo)
{
  
  boost::shared_ptr<SparseMatrix> tmp;
  if (orbs.size() == 0)
    tmp =   big->get_op_rep(ot, deltaQuantum);
  if (orbs.size() == 1)
    tmp =   big->get_op_rep(ot, deltaQuantum, orbs[0]);
  if (orbs.size() == 2)
    tmp =   big->get_op_rep(ot, deltaQuantum, orbs[0], orbs[1]);

  tmp->built = true;

  this->allocate(*newStateInfo);
  this->built = true;

  int newQ = 0;
  for (int Q = 0; Q < rotate_matrix.size (); ++Q)
    if (rotate_matrix[Q].Ncols () != 0)
      {
	int newQPrime = 0;
	for (int QPrime = 0; QPrime < rotate_matrix.size (); ++QPrime)
	  if (rotate_matrix[QPrime].Ncols () != 0)
	    {
	      if (this->allowedQuantaMatrix (newQ, newQPrime)) {
		MatrixRotate (rotate_matrix[Q], tmp->operatorMatrix(Q, QPrime), rotate_matrix[QPrime],
			      this->operatorMatrix (newQ, newQPrime) );
	      }
	      ++newQPrime;
	    }
	++newQ;
      }

}

SparseMatrix& SparseMatrix::operator+=(const SparseMatrix& other)
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
}
