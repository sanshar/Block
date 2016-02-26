#ifndef PARIO_HEADER_H
#define PARIO_HEADER_H
#include <communicate.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "global.h"

using namespace std;
using namespace SpinAdapted;


#ifdef MOLPRO
#define pm1out if (mpigetrank() == 0 && dmrginp.outputlevel() >= -1) bout
#define pm1err if (mpigetrank() == 0 && dmrginp.outputlevel() >= -1) berr
#endif

#define pout if (mpigetrank() == 0 && dmrginp.outputlevel() >= 0) bout
#define perr if (mpigetrank() == 0 && dmrginp.outputlevel() >= 0) berr

#define p1out if (mpigetrank() == 0 && dmrginp.outputlevel() >= 1) bout
#define p1err if (mpigetrank() == 0 && dmrginp.outputlevel() >= 1) berr

#define p2out if (mpigetrank() == 0 && dmrginp.outputlevel() >= 2) bout
#define p2err if (mpigetrank() == 0 && dmrginp.outputlevel() >= 2) berr

#define p3out if (mpigetrank() == 0 && dmrginp.outputlevel() >= 3) bout
#define p3err if (mpigetrank() == 0 && dmrginp.outputlevel() >= 3) berr

extern ostream &bout, &berr;

class blockout {
   public:
      ostream *outstream;
      char* output;
      blockout(ostream *outstream_ = &cout, char* output_=0): outstream(outstream_),output(output_)
      {
       if(output!=0) {
        ofstream file(output);
        outstream->rdbuf(file.rdbuf());
       }
      }
};

class blockerr {
   public:
      ostream *errstream;
      char* output;
      blockerr(ostream *errstream_ = &cerr, char* output_=0):errstream(errstream_), output(output_)
      {
       if(output!=0) {
        ofstream file(output);
        errstream->rdbuf(file.rdbuf());
       }
      }
};

#endif

