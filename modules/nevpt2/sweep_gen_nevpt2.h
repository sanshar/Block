/* 
 * File:   sweep_gen_nevpt2.h
 * Author: roemelt
 *
 * Created on February 13, 2014, 3:14 PM
 */

#ifndef SWEEP_GEN_NEVPT2_H
#define	SWEEP_GEN_NEVPT2_H

#include "initblocks.h"
#include "nevpt2_info.h"


namespace SpinAdapted{
  namespace nevpt2{
    void BlockAndDecimate_(SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys,NEVPT2Info &Info);
    double do_one_(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize);
  }
}

#endif	/* SWEEP_GEN_NEVPT2_H */

