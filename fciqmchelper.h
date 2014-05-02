#ifndef FCIQMC_HELPER_HEADER_H
#define FCIQMC_HELPER_HEADER_H

#include "global.h"
#include <vector>
#include "boost/shared_ptr.hpp"
#include "StateInfo.h"
#include "Operators.h"
#include "wavefunction.h"

namespace SpinAdapted{

//the MPS is stored in the left canonical form
//LLLLL..LC 
class MPS{

 private:
  friend class boost::serialization::access;
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & SiteTensors \
	& w ;
    }

  std::vector< std::vector<Matrix> > SiteTensors; //these are the L matrices
  Wavefunction w; //the last wavefunction

  void Init(std::vector<bool>& occ);
 public:
  static int sweepIters;
  static bool spinAdapted;
  static std::vector<SpinBlock> siteBlocks;

  MPS() {};
  MPS(int stateindex); 
  MPS(std::vector<bool>& occ);
  MPS(long* occnum, int length);
  std::vector<Matrix>& getSiteTensors(int i) {return SiteTensors[i];}
  const std::vector<Matrix>& getSiteTensors(int i) const {return SiteTensors[i];}
  const Wavefunction& getw() const {return w;}
  void normalize() {int success; w.Normalise(&success);}
};


 //statea is multiplied with Operator O|Mpsa> and then we compress it to get stateb
 //void compressOperatorTimesMPS(const MPS& statea, MPS& stateb);

 //calculate overlap between a and b <Mpsa|Mpsb>
 double calculateOverlap (const MPS& a, const MPS& b);

 //calculate hamiltonian matrix between a and b <Mpsa|H|Mpsb>
 void calcHamiltonianAndOverlap(const MPS& statea, const MPS& stateb, double& h, double& o) ;

}

#endif
