#ifndef SPIN_PARALLEL_IO_HEADER_H
#define SPIN_PARALLEL_IO_HEADER_H
#include <communicate.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#ifdef MOLPRO
#include "global/CxOutputStream.h"
#endif

class blockout {
   public:
      std::ostream *outstream;
      char* output;
      blockout(char* output_=0, std::ostream *outstream_ = &std::cout): output(output_),outstream(outstream_)
      {
       if(output!=0) {
        std::ofstream file(output);
        outstream->rdbuf(file.rdbuf());
       }
      }
};

class blockerr {
   public:
      std::ostream *errstream;
      char* output;
      blockerr(char* output_=0, std::ostream *errstream_ = &std::cerr): output(output_),errstream(errstream_)
      {
       if(output!=0) {
        std::ofstream file(output);
        errstream->rdbuf(file.rdbuf());
       }
      }
};


extern std::ostream &bout, &berr;

#define pout if (mpigetrank() == 0 && dmrginp.outputlevel() >= 0) bout
#define perr if (mpigetrank() == 0 && dmrginp.outputlevel() >= 0) berr

#endif

