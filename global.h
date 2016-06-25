/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_GLOBAL_HEADER
#define SPIN_GLOBAL_HEADER
 
#define WANT_MATH
#define WANT_STREAM
// Matrix library header files


#include "timer.h"
//#include <malloc.h>
#include <newmat.h>
#include <newmatap.h>
#include <newmatio.h>
#include <IntegralMatrix.h>
// STL headers
#include <iostream>
using namespace std;
#include <map>
#include <algorithm>
#include <vector>
#include <numeric>
// C headers
//#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
//#include "dyall_integral.h"
#include "Symmetry.h"
#include "input.h"
#include "pario.h"
#include <ctime>
#include <iostream>
// Global variables

namespace SpinAdapted{
class OneElectronArray;
class TwoElectronArray;
class PairArray;
class CCCCArray;
class CCCDArray;
//Perturb subspace in nevpt2.

extern int MAX_THRD;

extern Timer globaltimer;

extern std::vector<OneElectronArray> v_1;
extern std::vector<TwoElectronArray> v_2;
extern std::map<TwoPerturbType,PerturbTwoElectronArray> vpt_2;
extern OneElectronArray vpt_1;
//extern OneElectronArray vpt_1;
//TODO
//Now, just use 4d vector to store perturb integral. Change it inot a better one latter. 
//extern std::map<OnePerturbType,PerturbOneElectronArray> vpt_1;
//extern std::map<TwoPerturbType,PerturbTwoElectronArray> vpt_2;
//extern std::vector<OnePerturbArray> vpt_1;
//extern std::vector<TwoPerturbArray> vpt_2;
extern OneElectronArray fock;
extern std::vector<double> coreEnergy;
extern double BWPTenergy;
extern PairArray v_cc;
extern CCCCArray v_cccc;
extern CCCDArray v_cccd;

extern Input dmrginp;

extern bool SHOW_MORE;
extern bool DEBUG_MEMORY;
extern bool RESTART;
extern bool FULLRESTART;
extern bool BACKWARD;
extern bool reset_iter;
extern bool restartwarm;
extern string sym;
extern bool NonabelianSym;
extern std::vector<int> NPROP;
extern int PROPBITLEN;
extern double NUMERICAL_ZERO;
extern double tcpu,twall,ecpu,ewall;
//extern ifstream* coutbuf;
}
#endif
