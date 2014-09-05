/* 
 * File:   sweep_ripdm.h
 * Author: roemelt
 *
 * Created on April 5, 2013, 5:45 PM
 */

#ifndef SWEEP_NEVPT2_H
#define	SWEEP_NEVPT2_H

#include "global.h"
#include "solver.h"
#include "initblocks.h"
#include "rotationmat.h"
#include "davidson.h"
#include "linear.h"
#include "guess_wavefunction.h"
#include "nevpt2_operators.h"
#include "nevpt2_info.h"

#ifndef SERIAL
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif

using namespace boost;
using namespace std;

namespace SpinAdapted{
  namespace nevpt2{
    double do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize);
    //void BlockAndDecimate (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys);
    void BlockAndDecimate (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem,
                           const bool &useSlater, const bool& dot_with_sys,
                           ThreeIndOpArray &Ti, ThreeIndOpArray &Ta, IntegralContainer &IKJL,
                           IntegralContainer &IKJA, NEVPT2Info &Info);
    void nevpt2();
    void nevpt2_restart();
  };
}




#endif	/* SWEEP_NEVPT2_H */

