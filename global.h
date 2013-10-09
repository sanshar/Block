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
#ifndef KCC
using namespace std;
#endif
#include "timer.h"
#include "pario.h"
//#include <malloc.h>
#include <newmat.h>
#include <newmatap.h>
#include <newmatio.h>
#include <IntegralMatrix.h>
// STL headers
#include <iostream>
#include <map>
#include <algorithm>
#include <vector>
#include <numeric>
// C headers
//#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "Symmetry.h"
#include "input.h"
#include <ctime>
#include <iostream>
// Global variables

namespace SpinAdapted{
class OneElectronArray;
class TwoElectronArray;

extern int MAX_THRD;

extern Timer globaltimer;

extern OneElectronArray v_1;
extern TwoElectronArray v_2;

extern Input dmrginp;

extern bool SHOW_MORE;
extern bool DEBUG_MEMORY;
extern bool RESTART;
extern bool FULLRESTART;
extern bool reset_iter;
extern bool restartwarm;
extern string sym;
extern bool NonabelianSym;
extern std::vector<int> NPROP;
extern int PROPBITLEN;
//extern ifstream* coutbuf;
}
#endif
