/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef COUPLING_COEFFS_H
#define COUPLING_COEFFS_H
#include "math.h"
//#include "global.h"
#include <vector>
#include "include/newmatutils.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <multiarray.h>
//#include "anglib.h"
#include "new_anglib.h"
#include "Symmetry.h"

using namespace std;

namespace SpinAdapted{
typedef array_3d<double> array3d;

class ninejCoeffs{
  //initialise all the coeffs
  /*
    middle row represents the spin number of the operators based on these
    numbers different coefficient arrays are definted;
    1. 0 0 0
    2. 1 1 0
    3. 2 2 0
    4. 1 1 2
    5. 2 0 2
    6. 0 2 2
    7. 1 2 1
    8. 2 1 1
    9. 0 1 1
   10. 1 0 1
  */

 private:
  private:
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
      ar & maxj;
      ar & coeffs;
    }
  vector< array3d > coeffs;
  int maxj;
 
 public:
  ninejCoeffs():maxj(0){}
  ninejCoeffs(int maxj_);
  static ninejCoeffs& getinstance();
  void init(int maxj_);
  void buildArray();
  double operator()(int ja, int jb, int jc, int jd, int je, int jf, int jg, int jh, int ji) const;
  void initarray(int arrayindex, int a, int b, int c);
  int arrayindex (int a, int b, int c) const;
};


inline double cg(int two_ja, int two_jb, int two_jc, int two_ma, int two_mb, int two_mc)
{
  //double rval = cleb_(two_ja, two_ma, two_jb, two_mb, two_jc, two_mc);
  double rval = clebsch(two_ja, two_ma, two_jb, two_mb, two_jc, two_mc);
  return rval;
}


inline double sixj(int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf)
{
  //return sixj_(two_ja, two_jb, two_jc, two_jd, two_je, two_jf);
  return six_j(two_ja, two_jb, two_jc, two_jd, two_je, two_jf);
}


inline double racah(int a, int b, int c, int d, int e, int f)
{
  double rval = pow(-1.0, static_cast<int>((a+b+c+d)/2))*sixj(a, b, e, d, c, f);
  return rval;
}

inline double Ninej(int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf, int two_jg, int two_jh, int two_ji)
{

  double rval= pow((two_jg+1)*(two_jh+1)*(two_jc+1)*(two_jf+1),0.5)
    *nine_j(two_ja, two_jb, two_jc, two_jd, two_je, two_jf, two_jg, two_jh, two_ji);
    //*ninej_ (two_ja, two_jb, two_jc, 
	  //   two_jd, two_je, two_jf, 
	   //  two_jg, two_jh, two_ji);
  return rval;
}
}
#endif
