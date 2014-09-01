#ifndef SPIN_PARALLEL_IO_HEADER_H
#define SPIN_PARALLEL_IO_HEADER_H
#include <communicate.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#ifdef MOLPRO
#include "global/CxOutputStream.h"
#endif

#define pout if (mpigetrank() == 0) bout
#define perr if (mpigetrank() == 0) berr

extern std::ostream &bout, &berr;

class blockout {
   public:
      std::ostream *outstream;
      char* output;
      blockout(std::ostream *outstream_ = &std::cout, char* output_=0): output(output_),outstream(outstream_)
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
      blockerr(std::ostream *errstream_ = &std::cerr, char* output_=0): output(output_),errstream(errstream_)
      {
       if(output!=0) {
        std::ofstream file(output);
        errstream->rdbuf(file.rdbuf());
       }
      }
};

#endif

