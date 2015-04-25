/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "BaseOperator.h"
#include "MatrixBLAS.h"
#include "spinblock.h"
#include "distribute.h"
#include "tensor_operator.h"
#include "blas_calls.h"
#include "pario.h"

namespace SpinAdapted{


double getCommuteParity(SpinQuantum a, SpinQuantum b, SpinQuantum c)
{
  int aspin = a.get_s().getirrep(), airrep = a.get_symm().getirrep();
  int bspin = b.get_s().getirrep(), birrep = b.get_symm().getirrep();
  int cspin = c.get_s().getirrep(), cirrep = c.get_symm().getirrep();

  //int an = a.get_n(), bn = b.get_n();
  int parity = IsFermion(a) && IsFermion(b) ? -1 : 1;

  if (!dmrginp.spinAdapted()) return parity;

  for (int asz = -aspin; asz<aspin+1; asz+=2)
  for (int bsz = -bspin; bsz<bspin+1; bsz+=2)
  for (int al = 0; al<Symmetry::sizeofIrrep(airrep); al++)
  for (int bl = 0; bl<Symmetry::sizeofIrrep(birrep); bl++)
  {
    //double cleb = cleb_(aspin, asz, bspin, bsz, cspin, cspin);
    double cleb = clebsch(aspin, asz, bspin, bsz, cspin, cspin);
    double clebspatial = Symmetry::spatial_cg(airrep, birrep, cirrep, al, bl, 0);
    if (fabs(cleb) <= NUMERICAL_ZERO || fabs(clebspatial) <= NUMERICAL_ZERO)
      continue;
    else {
      //return parity*cleb*clebdinfh/cleb_(bspin, bsz, aspin, asz, cspin, cspin)/Symmetry::spatial_cg(birrep, airrep, cirrep, bl, al, 0);
      double spinscale = cleb/clebsch(bspin, bsz, aspin, asz, cspin, cspin);
      double spatscale = clebspatial/Symmetry::spatial_cg(birrep, airrep, cirrep, bl, al, 0);

      return parity*spinscale*spatscale;
    }
  }
  pout << "Major trouble, getCommuteParity asked for three inappropriate operators"<<endl;
  pout << a<<"  "<<b<<"  "<<c<<endl;
  abort();
  return 1.0;
}


double Transposeview::get_scaling(SpinQuantum leftq, SpinQuantum rightq) const 
{
  if(!dmrginp.spinAdapted()) return 1.0;
  if (conjugacy() == 'n') {return 1.0;}

  int lspin = leftq.get_s().getirrep(), lirrep = leftq.get_symm().getirrep();
  int rspin = rightq.get_s().getirrep(), rirrep = rightq.get_symm().getirrep();
  int cspin = opdata->get_deltaQuantum(0).get_s().getirrep(), cirrep = opdata->get_deltaQuantum(0).get_symm().getirrep();

  int cirrepTranspose = (-opdata->get_deltaQuantum(0)).get_symm().getirrep();
  for (int lsz = -lspin; lsz<lspin+1; lsz+=2)
  for (int rsz = -rspin; rsz<rspin+1; rsz+=2)
  for (int ll = 0; ll<Symmetry::sizeofIrrep(lirrep); ll++)
  for (int rl = 0; rl<Symmetry::sizeofIrrep(rirrep); rl++)
  {
    double cleb = clebsch(lspin, lsz, cspin, -cspin, rspin, rsz);
    double clebspatial = Symmetry::spatial_cg(lirrep, cirrep, rirrep, ll, 0, rl);
    if (fabs(cleb) <= NUMERICAL_ZERO || fabs(clebspatial) <= NUMERICAL_ZERO)
      continue;
    else {
      ///CHANGE THE SPATIAL_CG cirrep,1 to cirrep,0 depending on how the transpose works out!!!
      double spinscale = pow(-1.0,cspin) * cleb/clebsch(rspin, rsz, cspin, cspin, lspin, lsz);
      double spatscale =  clebspatial/Symmetry::spatial_cg(rirrep, cirrepTranspose, lirrep, rl, Symmetry::sizeofIrrep(cirrep)-1, ll);  

      return spinscale*spatscale;
    }
  }
  pout << "Major trouble, inappropriate arguments to get_scaling!!!"<<endl;
  pout << leftq<<"  ";
  for (int i = 0; i < get_deltaQuantum_size(); ++i) {
    pout <<get_deltaQuantum(i)<<"  ";
  }
  pout <<rightq<<endl;
  exit(0);
  return 1.0;
}


void SparseMatrix::makeIdentity(const StateInfo& s) 
{
  set_deltaQuantum(1, SpinQuantum (0, SpinSpace(0), IrrepSpace(0)));
  allocate(s);
  built = true;
  initialised = true;
  for (int i = 0; i < allowedQuantaMatrix.Nrows (); ++i)
    for (int j = 0; j < allowedQuantaMatrix.Ncols (); ++j)
      {
	if (allowedQuantaMatrix (i,j)) {
	  operatorMatrix(i,j) = 0.0;
	  for (int k=0; k<s.quantaStates.at(i); k++)
	    operatorMatrix(i,j)(k+1, k+1) = 1.0;
	}
      }
}


void SparseMatrix::allocate(const SpinBlock& b)
{
  allocate(b.get_stateInfo());
}


void SparseMatrix::allocate(const StateInfo& s)
{
  allocate(s, s);
}

void SparseMatrix::allocate(const StateInfo& sr, const StateInfo& sc)
{
  resize(sr.quanta.size(), sc.quanta.size());
  long totalmemory = 0;
  for (int i = 0; i < allowedQuantaMatrix.Nrows (); ++i)
    for (int j = 0; j < allowedQuantaMatrix.Ncols (); ++j)
    {
      bool allowedcoupling = false;
      for (int k = 0; k < deltaQuantum.size(); ++k) {
        if (sr.quanta[i].allow(deltaQuantum[k], sc.quanta[j])) {
          allowedcoupling = true;
          break;
        }
      }
      allowedQuantaMatrix (i,j) = allowedcoupling;
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
  deltaQuantum.clear();
  orbs.resize(0);
  allowedQuantaMatrix.ReSize (0,0);
  operatorMatrix.ReSize (0,0);
}

const Transposeview Transpose(SparseMatrix& op) { return Transposeview(op); };

ostream& operator<< (ostream& os, const SparseMatrix& a)
{
  if (!a.initialised){
    os <<" not initialised"<<endl;
    //return os;
  };
  os<<"indices : ";
  for(int i=0; i<a.orbs.size(); i++)
    os<<a.orbs[i]<<"  ";
  os <<endl;
  for (int i = 0; i < a.get_deltaQuantum_size(); ++i) {
    os<<a.get_deltaQuantum(i)<<endl;
  }
  for (int i = 0; i < a.nrows (); ++i)
	for (int j = 0; j < a.ncols (); ++j)
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
	    for (int jq =stateinfo.unBlockedIndex[j]; jq < stateinfo.unBlockedIndex[j]+stateinfo.quantaStates[j]; jq++) {
	      for (int iq =stateinfo.unBlockedIndex[i]; iq < stateinfo.unBlockedIndex[i]+stateinfo.quantaStates[i]; iq++) 
	        operatorMatrix(i,j)(iq-stateinfo.unBlockedIndex[i]+1, jq-stateinfo.unBlockedIndex[j]+1) = redMatrixElement(s[iq], ladders[jq], &b);
	    }

}

void SparseMatrix::Randomise ()
{
  for (int lQ = 0; lQ < nrows(); ++lQ)
    for (int rQ = 0; rQ < ncols(); ++rQ)
      if (allowed(lQ, rQ))
	    SpinAdapted::Randomise(operator_element(lQ, rQ));
}

void SparseMatrix::SymmetricRandomise ()
{
  for (int lQ = 0; lQ < nrows(); ++lQ)
    for (int rQ = 0; rQ < ncols(); ++rQ)
      if (allowed(lQ, rQ)) 
	SpinAdapted::SymmetricRandomise(operator_element(lQ, rQ));
}

double trace(const SparseMatrix& lhs)
{
  assert(lhs.nrows() == lhs.ncols());
  double trace = 0.0;

  for(int lQ=0;lQ<lhs.nrows();++lQ)
    if(lhs.allowed(lQ,lQ))
      for(int i=0;i<(lhs)(lQ,lQ).Nrows();++i)
	trace += (lhs)(lQ,lQ)(i+1,i+1);

  return trace;
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
	      pout <<"Not a valid addition"<<endl;
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

  //if the norm is really small then dont normalize??
  if(normalisation > NUMERICAL_ZERO)
    Scale(1./sqrt(normalisation), *this);
  else {
    pout << "\t\t\t Warning :: Didn't Normalise, because norm is " << normalisation<<endl;
    *success = 1; //not successful in normlaising                                                                                  
  }
}

void SparseMatrix::Clear ()
{
  built = false;
  for (int i = 0; i < allowedQuantaMatrix.Nrows (); ++i)
    for (int j = 0; j < allowedQuantaMatrix.Ncols (); ++j)
      if (allowedQuantaMatrix (i,j)) SpinAdapted::Clear(operatorMatrix (i,j));
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

//Renormalization functions for core and virtual operators                                                                                
void SparseMatrix::renormalise_transform(const std::vector<Matrix>& leftrotate_matrix, const StateInfo *leftstateinfo, const std::vector<Matrix>& rightrotate_matrix, const StateInfo *rightstateinfo)
{
  ObjectMatrix<Matrix> tmp = operatorMatrix; //cannot instantiate a SparseMatrix and so instantiating a Cre

  this->allocate(*leftstateinfo, *rightstateinfo); // new allocations                                                                                           

  int lQPrime = 0;
  for (int lQ = 0; lQ < leftrotate_matrix.size (); ++lQ)
    if (leftrotate_matrix[lQ].Ncols () != 0)
    {
      int rQPrime = 0;
      for (int rQ = 0; rQ < rightrotate_matrix.size (); ++rQ)
	if (rightrotate_matrix[rQ].Ncols () != 0)
	{
	  if (this->allowedQuantaMatrix (lQPrime, rQPrime))
	    MatrixRotate (leftrotate_matrix[lQ], tmp(lQ, rQ), rightrotate_matrix[rQ],
			  this->operatorMatrix (lQPrime, rQPrime) );
	  ++rQPrime;
	}
      ++lQPrime;
    }
  
}

void SparseMatrix::build_and_renormalise_transform(SpinBlock *big, const opTypes &ot, const std::vector<Matrix>& leftrotate_matrix, const StateInfo *newleftStateInfo, 
						   const std::vector<Matrix>& rightrotate_matrix, const StateInfo *newrightStateInfo)
{
  
  boost::shared_ptr<SparseMatrix> tmp;
  if (orbs.size() == 0)
    tmp =   big->get_op_rep(ot, deltaQuantum);
  if (orbs.size() == 1)
    tmp =   big->get_op_rep(ot, deltaQuantum, orbs[0]);
  if (orbs.size() == 2)
    tmp =   big->get_op_rep(ot, deltaQuantum, orbs[0], orbs[1]);

  tmp->built = true;

  this->allocate(*newleftStateInfo, *newrightStateInfo);
  this->built = true;

  int lQPrime = 0;
  for (int lQ = 0; lQ < leftrotate_matrix.size (); ++lQ)
    if (leftrotate_matrix[lQ].Ncols () != 0)
    {
      int rQPrime = 0;
      for (int rQ = 0; rQ < rightrotate_matrix.size (); ++rQ)
	if (rightrotate_matrix[rQ].Ncols () != 0)
	{
	  if (this->allowedQuantaMatrix (lQPrime, rQPrime))
	    MatrixRotate (leftrotate_matrix[lQ], tmp->operator_element(lQ, rQ), rightrotate_matrix[rQ],
			  this->operatorMatrix (lQPrime, rQPrime) );
	  ++rQPrime;
	}
      ++lQPrime;
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

SubSparseMatrix::SubSparseMatrix(SparseMatrix& op, int sec, const StateInfo& s): section(sec) {
  opdata = boost::shared_ptr<SparseMatrix>(&op, boostutils::null_deleter());
  SpinQuantum q = opdata->get_deltaQuantum(section);
  SuballowedQuantaMatrix.ReSize(s.quanta.size(), s.quanta.size());
  for (int i = 0; i < SuballowedQuantaMatrix.Nrows (); ++i)
    for (int j = 0; j < SuballowedQuantaMatrix.Ncols (); ++j) {
      if (s.quanta[i].allow(q, s.quanta[j])) {
        SuballowedQuantaMatrix(i,j) = true;
      }
    }
}

SubSparseMatrix::SubSparseMatrix(const boost::shared_ptr<SparseMatrix>& opptr, int sec, const StateInfo& s): opdata(opptr), section(sec) {
  SpinQuantum q = opdata->get_deltaQuantum(section);
  SuballowedQuantaMatrix.ReSize(s.quanta.size(), s.quanta.size());
  for (int i = 0; i < SuballowedQuantaMatrix.Nrows (); ++i)
    for (int j = 0; j < SuballowedQuantaMatrix.Ncols (); ++j) {
      if (s.quanta[i].allow(q, s.quanta[j])) {
        SuballowedQuantaMatrix(i,j) = true;
      }
    }
}

SubSparseMatrix::SubSparseMatrix(const boost::shared_ptr<SparseMatrix>& opptr, int sec, const StateInfo& sr, const StateInfo& sc): opdata(opptr), section(sec) {
  SpinQuantum q = opdata->get_deltaQuantum(section);
  SuballowedQuantaMatrix.ReSize(sr.quanta.size(), sc.quanta.size());
  for (int i = 0; i < SuballowedQuantaMatrix.Nrows (); ++i)
    for (int j = 0; j < SuballowedQuantaMatrix.Ncols (); ++j) {
      if (sr.quanta[i].allow(q, sc.quanta[j])) {
        SuballowedQuantaMatrix(i,j) = true;
      }
    }
}

}
