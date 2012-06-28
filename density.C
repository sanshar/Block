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
*/

#include "density.h"
#include "wavefunction.h"
#include "operatorloops.h"
#include "operatorfunctions.h"
#include <omp.h>
#include "guess_wavefunction.h"
#include "distribute.h"
#include <boost/format.hpp>

namespace SpinAdapted{
using namespace operatorfunctions;

void DensityMatrix::makedensitymatrix(const std::vector<Wavefunction>& wave_solutions, SpinBlock &big, 
				      const std::vector<double> &wave_weights, const double noise, const double additional_noise, bool warmup)
{
  for(int i=0;i<wave_weights.size();++i) {
    MultiplyProduct (wave_solutions[i], Transpose(const_cast<Wavefunction&> (wave_solutions[i])), *this, wave_weights[i]);
  }

  
  if(noise > 1.0e-14)
    this->add_onedot_noise(wave_solutions, big, noise);

}

void DensityMatrix::add_twodot_noise(const SpinBlock &big, const double noise)
{
  if (dmrginp.outputlevel() != 0) 
    pout << "\t\t\t adding noise " << noise << endl;
  double norm = 0.0;
  for(int lQ=0;lQ<this->nrows();++lQ)
    for(int rQ=0;rQ<this->ncols();++rQ)
      if(this->allowed(lQ,rQ))
        for(int i=0;i<(*this)(lQ,rQ).Nrows();++i)
          norm += (*this)(lQ,rQ)(i+1,i+1);
  if (dmrginp.outputlevel() != 0) 
    pout << "\t\t\t norm before modification " << norm << endl;

  Wavefunction noiseMatrix;
  double reweight = 0.;

  DensityMatrix noisedm = *this;
  noisedm.Clear();

  vector<SpinQuantum> toadd;

  {
    const int particlenumber = dmrginp.total_particle_number();
    const int spinnumber = dmrginp.total_spin_number();
    const IrrepSpace& symmetrynumber = dmrginp.total_symmetry_number();
    toadd.push_back(SpinQuantum(particlenumber + 1, spinnumber + 1, symmetrynumber));
    toadd.push_back(SpinQuantum(particlenumber - 1, spinnumber + 1, symmetrynumber));
    if (spinnumber >= 1) {
      toadd.push_back(SpinQuantum(particlenumber + 1, spinnumber - 1, symmetrynumber));
      toadd.push_back(SpinQuantum(particlenumber - 1, spinnumber - 1, symmetrynumber));
    }
    toadd.push_back(SpinQuantum(particlenumber + 2, spinnumber + 2, symmetrynumber));
    toadd.push_back(SpinQuantum(particlenumber, spinnumber + 2, symmetrynumber));
    toadd.push_back(SpinQuantum(particlenumber - 2, spinnumber + 2, symmetrynumber));

    toadd.push_back(SpinQuantum(particlenumber + 2, spinnumber, symmetrynumber));
    toadd.push_back(SpinQuantum(particlenumber - 2, spinnumber, symmetrynumber));

    if (spinnumber >= 2) {
      toadd.push_back(SpinQuantum(particlenumber + 2, spinnumber - 2, symmetrynumber));
      toadd.push_back(SpinQuantum(particlenumber, spinnumber - 2, symmetrynumber));
      toadd.push_back(SpinQuantum(particlenumber - 2, spinnumber - 2, symmetrynumber));
    }
  }

  for (int q = 0; q < toadd.size(); ++q)
  {
    noiseMatrix.initialise(toadd[q], &big, false);
    noiseMatrix.Randomise();
    double norm = DotProduct(noiseMatrix, noiseMatrix);
    if (abs(norm) > 1.e-20)
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
  if (dmrginp.outputlevel() != 0) 
    pout << "\t\t\t norm after modification " << norm << endl;

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
  DensityMatrix* dm;
  const SpinBlock& big; 
  const double scale;
  const int num_threads;
  opTypes optype, optype2;
  bool distributed;
  bool synced;
public:
  onedot_noise_f(DensityMatrix* dm_, const Wavefunction& wavefunction_, const SpinBlock& big_, const double scale_, const int num_threads_)
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
	SpinQuantum wQ = wavefunction.get_deltaQuantum();
	SpinQuantum oQ = op.get_deltaQuantum();
	vector<IrrepSpace> vec = wQ.get_symm() + oQ.get_symm();
	for (int j=0; j<vec.size(); j++)
	for (int i=abs(wQ.get_s()-oQ.get_s()); i<= wQ.get_s()+oQ.get_s(); i+=2)
	{
	  SpinQuantum q = SpinQuantum(wQ.get_n()+oQ.get_n(), i, vec[j]);

	  Wavefunction opxwave = Wavefunction(q, &big, wavefunction.get_onedot());
	  opxwave.Clear();
	  const boost::shared_ptr<SparseMatrix> fullop = op.getworkingrepresentation(big.get_leftBlock());
	  TensorMultiply(big.get_leftBlock(), *fullop, &big, const_cast<Wavefunction&> (wavefunction), opxwave, dmrginp.molecule_quantum(), 1.0);
	  double norm = DotProduct(opxwave, opxwave);
	  if (abs(norm) > 1e-14) {
	    Scale(1./sqrt(norm), opxwave);
	    MultiplyProduct(opxwave, Transpose(opxwave), dm[omp_get_thread_num()], scale);
	  }
	  q = SpinQuantum(wQ.get_n()-oQ.get_n(), i, vec[j]);

	  Wavefunction opxwave2 = Wavefunction(q, &big, wavefunction.get_onedot());
	  opxwave2.Clear();
	  TensorMultiply(big.get_leftBlock(),Transpose(*fullop),&big, const_cast<Wavefunction&> (wavefunction), opxwave2, dmrginp.molecule_quantum(), 1.0);
	  norm = DotProduct(opxwave2, opxwave2);
	  if (abs(norm) >1e-14) {
	    Scale(1./sqrt(norm), opxwave2);
	    MultiplyProduct(opxwave2, Transpose(opxwave2), dm[omp_get_thread_num()], scale);
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
// This is the noise that you set out in the schedule
void DensityMatrix::add_onedot_noise(const std::vector<Wavefunction>& wave_solutions, SpinBlock& big, const double noise, bool act2siteops)
{
/* check normalisation */
  double norm = 0.0;
  for(int lQ=0;lQ<this->nrows();++lQ)
    for(int rQ=0;rQ<this->ncols();++rQ)
      if(this->allowed(lQ,rQ))
        for(int i=0;i<(*this)(lQ,rQ).Nrows();++i)
          norm += (*this)(lQ,rQ)(i+1,i+1);
  if (dmrginp.outputlevel() != 0) 
    pout << "\t\t\t norm before modification " << norm << endl;

  SpinBlock* leftBlock = big.get_leftBlock();
  if (dmrginp.outputlevel() != 0) 
    pout << "\t\t\t Modifying density matrix " << endl;
    pout << "\t\t\t noise before modification " << noise << endl;
  //int maxt = 1;
  DensityMatrix* dmnoise = new DensityMatrix[MAX_THRD];
  for(int j=0;j<MAX_THRD;++j)
    dmnoise[j].allocate(big.get_leftBlock()->get_stateInfo());
  for(int i=0;i<wave_solutions.size();++i)
  {
    for(int j=0;j<MAX_THRD;++j)
      dmnoise[j].Clear();
    //Is it here where the noise gets added?
    onedot_noise_f onedot_noise(dmnoise, wave_solutions[i], big, 1., MAX_THRD);

    if (leftBlock->has(CRE))
    {
      onedot_noise.set_opType(CRE);
      for_all_multithread(leftBlock->get_op_array(CRE), onedot_noise);
    }
    
    if (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY) {
    if (leftBlock->has(CRE_CRE))
    {
      onedot_noise.set_opType(CRE_CRE);
      for_all_multithread(leftBlock->get_op_array(CRE_CRE), onedot_noise);

      onedot_noise.set_opType(CRE_DES);
      for_all_multithread(leftBlock->get_op_array(CRE_DES), onedot_noise);

    }

    else if (leftBlock->has(DES_DESCOMP))
    {
      onedot_noise.set_opType(DES_DESCOMP);
      for_all_multithread(leftBlock->get_op_array(DES_DESCOMP), onedot_noise);

      onedot_noise.set_opType(CRE_DESCOMP);
      for_all_multithread(leftBlock->get_op_array(CRE_DESCOMP), onedot_noise);

    }
    }

    onedot_noise.syncaccumulate();
    norm = 0.0;
    for(int lQ=0;lQ<dmnoise[0].nrows();++lQ)
      for(int rQ=0;rQ<dmnoise[0].ncols();++rQ)
	if(this->allowed(lQ,rQ))
	  for(int i=0;i<(dmnoise[0])(lQ,rQ).Nrows();++i)
	    norm += (dmnoise[0])(lQ,rQ)(i+1,i+1);
    if (norm > 1.0)
      ScaleAdd(noise/norm, dmnoise[0], *this);
  }
  delete[] dmnoise;  

  norm = 0.0;
  for(int lQ=0;lQ<this->nrows();++lQ)
    for(int rQ=0;rQ<this->ncols();++rQ)
      if(this->allowed(lQ,rQ))
        for(int i=0;i<(*this)(lQ,rQ).Nrows();++i)
          norm += (*this)(lQ,rQ)(i+1,i+1);
  if (dmrginp.outputlevel() != 0) 
    pout << "\t\t\t norm after modification " << norm << endl;

}

}
  
