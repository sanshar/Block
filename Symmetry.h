#ifndef SPIN_SYMMETRY_HEADER
#define SPIN_SYMMETRY_HEADER
#include <fstream>
#include <cstdio>
#include <multiarray.h>
#include <boost/serialization/serialization.hpp>
#include <vector>
#include "IrrepVector.h"

using namespace boost;
/** Non-Abelian and Abelian Symmetry class */

namespace SpinAdapted {
class IrrepVector;

class Symmetry
{
 public:
  static void InitialiseTable(string sym);

  static std::vector<int> add(int irrepl, int irrepr);
  static int sizeofIrrep(int irrep);

  static double spatial_cg(int a, int b, int c, int la, int lb, int lc) ;
  static double spatial_sixj(int a, int b, int c, int d, int e, int f) ;
  static double spatial_ninej(int j1, int j2, int j12, int j3, int j4, int j34, int j13, int j24, int j);

};
}
#endif
