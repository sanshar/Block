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
  void add_twodot_noise(const SpinBlock &big, const double noise);
  void add_onedot_noise(const std::vector<Wavefunction>& wave_solutions, SpinBlock& big, const double noise, bool act2siteops = true);
  void add_simulatedtwodot_noise(const std::vector<Wavefunction>& wave_solutions, SpinBlock& big, const double noise);
public:
  DensityMatrix()
  {
    orbs = std::vector<int>();
    initialised = true;
    fermion = false;
    deltaQuantum = SpinQuantum(0, SpinSpace(0), IrrepSpace(0));
  }
  void makedensitymatrix(const std::vector<Wavefunction>& wave_solutions, SpinBlock &big, const std::vector<double> &wave_weights,
			 const double noise, const double additional_noise, bool warmup);
  DensityMatrix& operator+=(const DensityMatrix& other);

  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block) {return boost::shared_ptr<SparseMatrix> (this);}
  void build(const SpinBlock& b){};
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b){return 0.0;}
};
}

#endif
