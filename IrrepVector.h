#ifndef IRREP_VECTOR_HEADER
#define IRREP_VECTOR_HEADER
#include <fstream>
#include <cstdio>
#include <multiarray.h>
#include <boost/serialization/serialization.hpp>
#include <vector>
#include "Symmetry.h"
#include <string>

using namespace boost;
/** Non-Abelian and Abelian Symmetry class */

namespace SpinAdapted {
  extern string sym;
class IrrepVector
{
  int irrep; //this is the irrep of the vector space this vector belongs to
  int row; //this is the row of the vector space that this vector transforms as
  
 public:
  IrrepVector() : irrep(0), row(0) {}
  IrrepVector(int pirrep, int prow) : irrep(pirrep), row(prow) {}
  int getirrep() const {return irrep;}
  int getrow() const {return row;}

  friend ostream& operator<< (ostream& os, const IrrepVector& iv) 
  {
    if (sym == "dinfh") {
      char goru = iv.irrep%2 == 0 ? 'g' : 'u';
      os<< max(0, (iv.irrep-2)/2)<<goru;
      if (iv.irrep <2) os<< '+';
      else if (iv.irrep >=2 && iv.irrep <4 ) os<< '-';
      else if (iv.row == 0) os<< '-' ;
      else os<<'+';
    }
    else {
      os<< iv.irrep+1;
    }
    return os;
  }
      

};
}

#endif
