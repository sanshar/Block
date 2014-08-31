#ifndef SPIN_PARALLEL_IO_HEADER_H
#define SPIN_PARALLEL_IO_HEADER_H
#include <communicate.h>
#include <stdio.h>
#include <iostream>

#ifndef MOLPRO
#define pout if (mpigetrank() == 0 && dmrginp.outputlevel() >= 0) cout
#define dout if (mpigetrank() == 0 && dmrginp.outputlevel() >= 0) cout
#else
#include "global/CxOutputStream.h"
//#define pout if (mpigetrank() == 0) xout
#define pout if ((mpigetrank() == 0) && (dmrginp.outputlevel() != 0)) xout
#define dout if ((mpigetrank() == 0) && (dmrginp.outputlevel() != 0)) xout
extern std::ostream &xout, &xerr;
#endif

#endif
