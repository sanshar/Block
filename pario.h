#ifndef SPIN_PARALLEL_IO_HEADER_H
#define SPIN_PARALLEL_IO_HEADER_H
#include <communicate.h>
#include <stdio.h>
#include <iostream>

#define pout if (mpigetrank() == 0) cout
#define dout if (mpigetrank() == 0) cout

#endif
