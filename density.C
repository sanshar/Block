/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "density.h"
#include "wavefunction.h"
#include "operatorloops.h"
#include "operatorfunctions.h"
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0 
#endif
#include "guess_wavefunction.h"
#include "distribute.h"
#include <boost/format.hpp>
#include "pario.h"


namespace SpinAdapted{
using namespace operatorfunctions;

void DensityMatrix::makedensitymatrix(const std::vector<Wavefunction>& wave_solutions, SpinBlock &big, 
				      const std::vector<double> &wave_weights, const double noise, const double additional_noise, bool warmup)
{

  for(int i=0;i<wave_weights.size()&& mpigetrank() == 0;++i) {
    makedensitymatrix(wave_solutions[i], big, wave_weights[i]);
  }

#ifndef SERIAL
  boost::mpi::communicator world;
  boost::mpi::broadcast(world, *this, 0);
#endif

  if(noise > NUMERICAL_ZERO) {
    
    /* check normalisation */
    double norm = 0.0;
    for(int lQ=0;lQ<nrows();++lQ)
      if(allowed(lQ,lQ))
	for(int i=0;i<(*this)(lQ,lQ).Nrows();++i)
	  norm += (*this)(lQ,lQ)(i+1,i+1);
    p2out << "\t\t\t norm before modification " << norm << endl;
    
    
    int nroots = wave_solutions.size();
    
#ifndef SERIAL
    boost::mpi::communicator world;
    boost::mpi::broadcast(world, nroots, 0);
#endif
    
    for (int i=0; i<nroots; i++) {
      this->add_onedot_noise(wave_solutions[i], big, (1.0*noise)/nroots);
    }

    norm = 0.0;
    for(int lQ=0;lQ<nrows();++lQ)
      if(this->allowed(lQ,lQ))
	for(int i=0;i<(*this)(lQ,lQ).Nrows();++i)
	  norm += (*this)(lQ,lQ)(i+1,i+1);
    p2out << "\t\t\t norm after modification " << norm << endl;
    
  }

  if(additional_noise > NUMERICAL_ZERO) {
    int nroots = wave_solutions.size();
    
#ifndef SERIAL
    boost::mpi::communicator world;
    boost::mpi::broadcast(world, nroots, 0);
#endif
    for (int i=0; i<nroots; i++) {
      this->add_twodot_noise(big, (1.0*noise)/nroots);
    }
  }
}
  
void DensityMatrix::makedensitymatrix(const Wavefunction& wave_solution, SpinBlock &big, 
				      const double &wave_weight)
{
  
  MultiplyProduct (wave_solution, Transpose(const_cast<Wavefunction&> (wave_solution)), *this, wave_weight);
  
}

void DensityMatrix::add_twodot_noise(const SpinBlock &big, const double noise)
{
  p1out << "\t\t\t adding noise " << noise << endl;
  double norm = 0.0;
  for(int lQ=0;lQ<this->nrows();++lQ)
    for(int rQ=0;rQ<this->ncols();++rQ)
      if(this->allowed(lQ,rQ))
        for(int i=0;i<(*this)(lQ,rQ).Nrows();++i)
          norm += (*this)(lQ,rQ)(i+1,i+1);
  p2out << "\t\t\t norm before modification " << norm << endl;

  Wavefunction noiseMatrix;
  double reweight = 0.;

  DensityMatrix noisedm = *this;
  noisedm.Clear();

  vector<SpinQuantum> toadd;

  {
    int particlenumber = dmrginp.total_particle_number();
    if (dmrginp.hamiltonian() == BCS) {
      particlenumber /= 2;
    } // FIXME do I need to use all particle numbers when doing BCS?
    const int spinnumber = dmrginp.total_spin_number().getirrep();
    const IrrepSpace& symmetrynumber = dmrginp.total_symmetry_number();
    toadd.push_back(SpinQuantum(particlenumber + 1, SpinSpace(spinnumber + 1), symmetrynumber));
    toadd.push_back(SpinQuantum(particlenumber - 1, SpinSpace(spinnumber + 1), symmetrynumber));
    if (spinnumber >= 1) {
      toadd.push_back(SpinQuantum(particlenumber + 1, SpinSpace(spinnumber - 1), symmetrynumber));
      toadd.push_back(SpinQuantum(particlenumber - 1, SpinSpace(spinnumber - 1), symmetrynumber));
    }
    toadd.push_back(SpinQuantum(particlenumber + 2, SpinSpace(spinnumber + 2), symmetrynumber));
    toadd.push_back(SpinQuantum(particlenumber, SpinSpace(spinnumber + 2), symmetrynumber));
    toadd.push_back(SpinQuantum(particlenumber - 2, SpinSpace(spinnumber + 2), symmetrynumber));

    toadd.push_back(SpinQuantum(particlenumber + 2, SpinSpace(spinnumber), symmetrynumber));
    toadd.push_back(SpinQuantum(particlenumber - 2, SpinSpace(spinnumber), symmetrynumber));

    if (spinnumber >= 2) {
      toadd.push_back(SpinQuantum(particlenumber + 2, SpinSpace(spinnumber - 2), symmetrynumber));
      toadd.push_back(SpinQuantum(particlenumber, SpinSpace(spinnumber - 2), symmetrynumber));
      toadd.push_back(SpinQuantum(particlenumber - 2, SpinSpace(spinnumber - 2), symmetrynumber));
    }
  }

  for (int q = 0; q < toadd.size(); ++q)
  {
    noiseMatrix.initialise(toadd[q], &big, false);
    noiseMatrix.Randomise();
    double norm = DotProduct(noiseMatrix, noiseMatrix);
    if (abs(norm) > NUMERICAL_ZERO)
    {
      Scale(1./sqrt(norm), noiseMatrix);
      MultiplyProduct(noiseMatrix, Transpose(noiseMatrix), noisedm, noise/toadd.size());
      noiseMatrix.CleanUp();
    }
    else
    {
      noiseMatrix.CleanUp();
      //pout << "\t\t\t no noise for quantum " << toadd[q] << endl;
    }
  }
  //Scale(1. - reweight, *this);
  *this += noisedm;
  norm = 0.0;
  for(int lQ=0;lQ<this->nrows();++lQ)
    for(int rQ=0;rQ<this->ncols();++rQ)
      if(this->allowed(lQ,rQ))
        for(int i=0;i<(*this)(lQ,rQ).Nrows();++i)
          norm += (*this)(lQ,rQ)(i+1,i+1);
  p2out << "\t\t\t norm after modification " << norm << endl;

}

DensityMatrix& DensityMatrix::operator+=(const DensityMatrix& other)
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



class onedot_noise_f 
{
private:
  const Wavefunction& wavefunction;
  vector<DensityMatrix>& dm;
  const SpinBlock& big; 
  const double scale;
  const int num_threads;
  opTypes optype, optype2;
  bool distributed;
  bool synced;
public:
  onedot_noise_f(vector<DensityMatrix>& dm_, const Wavefunction& wavefunction_, const SpinBlock& big_, const double scale_, const int num_threads_)
    : distributed(false), synced(true), wavefunction(wavefunction_), dm(dm_), big(big_), scale(scale_), num_threads(num_threads_) { }
  
  void set_opType(const opTypes &optype_)
  {
    optype = optype_;
    distributed = !big.get_leftBlock()->get_op_array(optype).is_local();
    if(distributed) synced = false;
  }
  void operator()(const std::vector<boost::shared_ptr<SparseMatrix> >& opvec) const
  {
    if ((mpigetrank() == 0 || distributed))// && op.get_deltaQuantum().particleNumber != 0)
    {
      for (int opind=0; opind<opvec.size(); opind++) {
	SparseMatrix& op = *opvec[opind];
#ifndef SERIAL
	boost::mpi::communicator world;
	if (op.get_orbs().size() == 1 && op.get_orbs()[0]%world.size() != mpigetrank())
	  continue;
#endif	  
	vector<SpinQuantum> wQ = wavefunction.get_deltaQuantum();
	vector<SpinQuantum> oQ = op.get_deltaQuantum();
	vector<IrrepSpace> vec = wQ[0].get_symm() + oQ[0].get_symm();
	vector<SpinSpace> spinvec = wQ[0].get_s()+oQ[0].get_s();
    for (int k=0; k<wQ.size(); ++k)
    for (int l=0; l<oQ.size(); ++l)
	for (int j=0; j<vec.size(); j++)
	for (int i=0; i<spinvec.size(); i++)
	{
	  SpinQuantum q = SpinQuantum(wQ[k].get_n()+oQ[l].get_n(), spinvec[i], vec[j]);
	  const boost::shared_ptr<SparseMatrix> fullop = op.getworkingrepresentation(big.get_leftBlock());      
	  if (dmrginp.hamiltonian() != BCS || q.get_n() <= dmrginp.effective_molecule_quantum().get_n()) {
	    Wavefunction opxwave = Wavefunction(q, &big, wavefunction.get_onedot());
	    opxwave.Clear();
	    TensorMultiply(big.get_leftBlock(), *fullop, &big, const_cast<Wavefunction&> (wavefunction), opxwave, dmrginp.molecule_quantum(), 1.0);
	    double norm = DotProduct(opxwave, opxwave);
	    if (abs(norm) > NUMERICAL_ZERO) {
	      Scale(1./sqrt(norm), opxwave);
	      MultiplyProduct(opxwave, Transpose(opxwave), dm[omp_get_thread_num()], scale);
	    }
	  }
	  q = SpinQuantum(wQ[k].get_n()-oQ[l].get_n(), spinvec[i], vec[j]);
	  if (dmrginp.hamiltonian() != BCS || q.get_n() >= 0) {
	    Wavefunction opxwave2 = Wavefunction(q, &big, wavefunction.get_onedot());
	    opxwave2.Clear();
	    TensorMultiply(big.get_leftBlock(),Transpose(*fullop),&big, const_cast<Wavefunction&> (wavefunction), opxwave2, dmrginp.molecule_quantum(), 1.0);
	    double norm = DotProduct(opxwave2, opxwave2);
	    if (abs(norm) >NUMERICAL_ZERO) {
	      Scale(1./sqrt(norm), opxwave2);
	      MultiplyProduct(opxwave2, Transpose(opxwave2), dm[omp_get_thread_num()], scale);
	    } 
	  }   
	 }
      } 
    }
  }

  void syncaccumulate(int toproc = 0)
  {
    for(int i=1;i<num_threads;++i)
      dm[0] += dm[i];

    distributedaccumulate(dm[0]);
    synced = true;
  }
};


class onedot_noise_f_compression 
{
private:
  const Wavefunction& wavefunction;
  vector<DensityMatrix>& dm;
  const SpinBlock& big; 
  const double scale;
  const int num_threads;
  opTypes optype, optype2;
  bool distributed;
  bool synced;
public:
  onedot_noise_f_compression(vector<DensityMatrix>& dm_, const Wavefunction& wavefunction_, const SpinBlock& big_, const double scale_, const int num_threads_)
    : distributed(false), synced(true), wavefunction(wavefunction_), dm(dm_), big(big_), scale(scale_), num_threads(num_threads_) { }
  
  void set_opType(const opTypes &optype_)
  {
    optype = optype_;
    distributed = !big.get_leftBlock()->get_op_array(optype).is_local();
    if(distributed) synced = false;
  }
  void operator()(const std::vector<boost::shared_ptr<SparseMatrix> >& opvec) const
  {
    if ((mpigetrank() == 0 || distributed))// && op.get_deltaQuantum().particleNumber != 0)
    {
      for (int opind=0; opind<opvec.size(); opind++) {
	SparseMatrix& op = *opvec[opind];
#ifndef SERIAL
	boost::mpi::communicator world;
	if (op.get_orbs().size() == 1 && op.get_orbs()[0]%world.size() != mpigetrank())
	  continue;
#endif	  
	vector<SpinQuantum> wQ = wavefunction.get_deltaQuantum();
	vector<SpinQuantum> oQ = op.get_deltaQuantum();
	vector<IrrepSpace> vec = wQ[0].get_symm() + oQ[0].get_symm();
	vector<SpinSpace> spinvec = wQ[0].get_s()+oQ[0].get_s();
	for (int k=0; k<wQ.size(); ++k)
	  for (int l=0; l<oQ.size(); ++l)
	    for (int j=0; j<vec.size(); j++)
	      for (int i=0; i<spinvec.size(); i++)
		{
		  SpinQuantum q = SpinQuantum(wQ[k].get_n()+oQ[l].get_n(), spinvec[i], vec[j]);
		  const boost::shared_ptr<SparseMatrix> fullop = op.getworkingrepresentation(big.get_leftBlock());      
		  if (dmrginp.hamiltonian() != BCS || q.get_n() <= dmrginp.effective_molecule_quantum().get_n()) {
		    Wavefunction opxwave;
		    opxwave.AllowQuantaFor(*big.get_braStateInfo().leftStateInfo, *big.get_ketStateInfo().rightStateInfo, std::vector<SpinQuantum>(1,q));opxwave.set_onedot(wavefunction.get_onedot());
		    
		    opxwave.Clear();
		    TensorMultiply(big.get_leftBlock(), *fullop, &big, const_cast<Wavefunction&> (wavefunction), opxwave, oQ[l], 1.0);
		    double norm = DotProduct(opxwave, opxwave);
		    if (abs(norm) > NUMERICAL_ZERO) {
		      Scale(1./sqrt(norm), opxwave);
		      MultiplyProduct(opxwave, Transpose(opxwave), dm[omp_get_thread_num()], scale);
		    }
		  }
		}
      }
    }
  }

  void syncaccumulate(int toproc = 0)
  {
    for(int i=1;i<num_threads;++i)
      dm[0] += dm[i];

    distributedaccumulate(dm[0]);
    synced = true;
  }
};

// accumulates into dm
void DensityMatrix::add_onedot_noise(const Wavefunction& wave_solution, SpinBlock& big, const double noise, bool act2siteops)
{

  SpinBlock* leftBlock = big.get_leftBlock();
  p1out << "\t\t\t Modifying density matrix " << endl;
  //int maxt = 1;
  vector<DensityMatrix> dmnoise(MAX_THRD, DensityMatrix(big.get_leftBlock()->get_stateInfo()));
  for(int j=0;j<MAX_THRD;++j)
    dmnoise[j].allocate(big.get_leftBlock()->get_stateInfo());

  for(int j=0;j<MAX_THRD;++j)
    dmnoise[j].Clear();

  Wavefunction wave;
  Wavefunction *wvptr = &wave;
  if(mpigetrank() == 0)
    wvptr = const_cast<Wavefunction*> (&wave_solution);
  else
    wvptr = &wave;

#ifndef SERIAL
  boost::mpi::communicator world;
  boost::mpi::broadcast(world, *wvptr, 0);
#endif


    onedot_noise_f onedot_noise(dmnoise, *wvptr, big, 1., MAX_THRD);

    if (leftBlock->has(CRE)) {
      onedot_noise.set_opType(CRE);
      for_all_multithread(leftBlock->get_op_array(CRE), onedot_noise);
    }
    
    if (dmrginp.hamiltonian() != HUBBARD) {

      if (leftBlock->has(CRE_CRE)) {
        onedot_noise.set_opType(CRE_CRE);
        for_all_multithread(leftBlock->get_op_array(CRE_CRE), onedot_noise);

        onedot_noise.set_opType(CRE_DES);
        for_all_multithread(leftBlock->get_op_array(CRE_DES), onedot_noise);
      } else if (leftBlock->has(DES_DESCOMP)) {
        onedot_noise.set_opType(DES_DESCOMP);
        for_all_multithread(leftBlock->get_op_array(DES_DESCOMP), onedot_noise);

        onedot_noise.set_opType(CRE_DESCOMP);
        for_all_multithread(leftBlock->get_op_array(CRE_DESCOMP), onedot_noise);

      }
      
      
      onedot_noise.syncaccumulate();
      double norm = 0.0;
      for(int lQ=0;lQ<dmnoise[0].nrows();++lQ)
	if(this->allowed(lQ,lQ))
	  for(int i=0;i<(dmnoise[0])(lQ,lQ).Nrows();++i)
	    norm += (dmnoise[0])(lQ,lQ)(i+1,i+1);
      if (norm > 1.0)
	ScaleAdd(noise/norm, dmnoise[0], *this);
    }
    
    double norm = 0.0;
    for(int lQ=0;lQ<nrows();++lQ)
      if(this->allowed(lQ,lQ))
        for(int i=0;i<(*this)(lQ,lQ).Nrows();++i)
          norm += (*this)(lQ,lQ)(i+1,i+1);
    p2out << "\t\t\t norm after modification " << norm << endl;
}

// accumulates into dm
void DensityMatrix::add_onedot_noise_forCompression(const Wavefunction& wave_solution, SpinBlock& big, const double noise)
{

  SpinBlock* leftBlock = big.get_leftBlock();
  p1out << "\t\t\t Modifying density matrix " << endl;
  //int maxt = 1;
  vector<DensityMatrix> dmnoise(MAX_THRD, DensityMatrix(big.get_leftBlock()->get_braStateInfo()));
  for(int j=0;j<MAX_THRD;++j)
    dmnoise[j].allocate(big.get_leftBlock()->get_braStateInfo());
  
  for(int j=0;j<MAX_THRD;++j)
    dmnoise[j].Clear();
  
  Wavefunction wave;
  Wavefunction *wvptr = &wave;
  if(mpigetrank() == 0)
    wvptr = const_cast<Wavefunction*> (&wave_solution);
  else
    wvptr = &wave;

#ifndef SERIAL
  boost::mpi::communicator world;
  boost::mpi::broadcast(world, *wvptr, 0);
#endif
  
  onedot_noise_f_compression onedot_noise(dmnoise, *wvptr, big, 1., MAX_THRD);
  
  if (leftBlock->has(CRE)) {
    onedot_noise.set_opType(CRE);
    for_all_multithread(leftBlock->get_op_array(CRE), onedot_noise);
  }
  if (leftBlock->has(DES)) {
    onedot_noise.set_opType(DES);
    for_all_multithread(leftBlock->get_op_array(DES), onedot_noise);
  }
  if (leftBlock->has(OVERLAP)) {
    onedot_noise.set_opType(OVERLAP);
    for_all_multithread(leftBlock->get_op_array(OVERLAP), onedot_noise);
  }
  
  if (dmrginp.hamiltonian() != HUBBARD) {
    
    if (leftBlock->has(CRE_CRE)) {
      onedot_noise.set_opType(CRE_CRE);
      for_all_multithread(leftBlock->get_op_array(CRE_CRE), onedot_noise);
      
      onedot_noise.set_opType(CRE_DES);
      for_all_multithread(leftBlock->get_op_array(CRE_DES), onedot_noise);
    }
    else if (leftBlock->has(DES_DESCOMP)) {
      onedot_noise.set_opType(DES_DESCOMP);
      for_all_multithread(leftBlock->get_op_array(DES_DESCOMP), onedot_noise);
      
      onedot_noise.set_opType(CRE_DESCOMP);
      for_all_multithread(leftBlock->get_op_array(CRE_DESCOMP), onedot_noise);
      
    }

    if (leftBlock->has(DES_DES)) {
      onedot_noise.set_opType(DES_DES);
      for_all_multithread(leftBlock->get_op_array(DES_DES), onedot_noise);
      
      onedot_noise.set_opType(DES_CRE);
      for_all_multithread(leftBlock->get_op_array(DES_CRE), onedot_noise);
    }
    else if (leftBlock->has(CRE_CRECOMP)) {
      onedot_noise.set_opType(CRE_CRECOMP);
      for_all_multithread(leftBlock->get_op_array(CRE_CRECOMP), onedot_noise);
      
      onedot_noise.set_opType(DES_CRECOMP);
      for_all_multithread(leftBlock->get_op_array(DES_CRECOMP), onedot_noise);
      
    }

  }
  
  onedot_noise.syncaccumulate();
  double norm = 0.0;
  for(int lQ=0;lQ<dmnoise[0].nrows();++lQ)
    if(this->allowed(lQ,lQ))
      for(int i=0;i<(dmnoise[0])(lQ,lQ).Nrows();++i)
	norm += (dmnoise[0])(lQ,lQ)(i+1,i+1);
  if (norm > 1.0)
    ScaleAdd(noise/norm, dmnoise[0], *this);


}

}
  
