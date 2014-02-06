/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef NPDM_HEADER
#define NPDM_HEADER

#ifndef SERIAL
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif

#include "sweepgenblock.h"
#include "sweep_params.h"
#include "sweeponepdm.h"
#include "npdm_driver.h"
#include "twopdm_driver.h"
#include "threepdm_driver.h"
#include "fourpdm_driver.h"
#include "nevpt2_pdm_driver.h"

//FIXME
#include <boost/mpi.hpp>
#include <vector>
#include <multiarray.h>
#include "spinblock.h"
#include "wavefunction.h"
#include "BaseOperator.h"

#include "global.h"
#include "solver.h"
#include "initblocks.h"
#include "rotationmat.h"
#include "davidson.h"
#include "linear.h"
#include "guess_wavefunction.h"
#include "density.h"
#include "davidson.h"
#include "pario.h"

using namespace boost;
using namespace std;

namespace SpinAdapted{
//FIXME Npdm namespace
namespace Npdm{

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  void npdm(int npdm_order);
  void npdm_restart(int npdm_order);

  void npdm_block_and_decimate( Npdm_driver& npdm_driver, SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, 
                                const bool &useSlater, const bool& dot_with_sys, const int state);

  double npdm_do_one_sweep(Npdm_driver& npdm_driver, SweepParams &sweepParams, const bool &warmUp, const bool &forward, 
                           const bool &restart, const int &restartSize, const int state);

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

}
}

#endif
