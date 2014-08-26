/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef SPIN_DENSITY_HEADER
#define SPIN_DENSITY_HEADER
#include "Operators.h"
#include "spinblock.h"
#include "wavefunction.h"

namespace SpinAdapted{
class DensityMatrix : public SpinAdapted::SparseMatrix
{
private:
  void add_onedot_noise(const Wavefunction& wave_solutions, SpinBlock& big, const double noise, bool act2siteops = true);
public:
  void add_twodot_noise(const SpinBlock &big, const double noise);
  void add_onedot_noise_forCompression(const Wavefunction& wave_solutions, SpinBlock& big, const double noise);
  DensityMatrix()
  {
    orbs = std::vector<int>();
    initialised = true;
    fermion = false;
    deltaQuantum.assign(1, SpinQuantum(0, SpinSpace(0), IrrepSpace(0)));
  }
  DensityMatrix(const StateInfo& s) {
    orbs = std::vector<int>();
    initialised = true;
    fermion = false;
    deltaQuantum.assign(1, SpinQuantum(0, SpinSpace(0), IrrepSpace(0)));   
    if (dmrginp.hamiltonian() == BCS) {
      int n_min = 100000, n_max = 0;
      for (int i = 0; i < s.quanta.size(); ++i) {
        if (s.quanta[i].get_n() > n_max)
          n_max = s.quanta[i].get_n();
        if (s.quanta[i].get_n() < n_min)
          n_min = s.quanta[i].get_n();
      }
      for (int dn=2; dn<=n_max-n_min; dn+=2) {
        deltaQuantum.push_back(SpinQuantum(dn, SpinSpace(0), IrrepSpace(0)));
        deltaQuantum.push_back(SpinQuantum(-dn, SpinSpace(0), IrrepSpace(0)));
      }
    }
  }
  void makedensitymatrix(const std::vector<Wavefunction>& wave_solutions, SpinBlock &big, const std::vector<double> &wave_weights,
			 const double noise, const double additional_noise, bool warmup);
  void makedensitymatrix(const Wavefunction& wave_solutions, SpinBlock &big, const double &wave_weight);
  DensityMatrix& operator+=(const DensityMatrix& other);

  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block) {return boost::shared_ptr<SparseMatrix> (this);}
  void build(const SpinBlock& b){};
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b){return 0.0;}
};
}

#endif
