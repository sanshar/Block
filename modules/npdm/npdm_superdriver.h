///*                                                                           
//Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
//Copyright (c) 2012, Garnet K.-L. Chan                                        
//                                                                             
//This program is integrated in Molpro with the permission of 
//Sandeep Sharma and Garnet K.-L. Chan
//*/
//
//#ifndef NPDM_SUPERDRIVER_HEADER_H
//#define NPDM_SUPERDRIVER_HEADER_H
//
//#include <vector>
//#include <multiarray.h>
//#include "spinblock.h"
//#include "wavefunction.h"
//#include "BaseOperator.h"
//
//namespace SpinAdapted{
//
////===========================================================================================================================================================
//// We wrote this base class because it only has two methods of its own, which are very similar to many in the legacy Block code
//
//class Npdm_superdriver {
//
//  public:
//    Npdm_superdriver() : use_sparse_npdm_(true) {}
//    double do_one_sweep(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize, int const state);
//
//  protected:
//    bool use_sparse_npdm_;
//
//  private:
//    void npdm_block_and_decimate(SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, 
//                                 const bool &useSlater, const bool& dot_with_sys, const int state);
//
//    virtual void compute_npdm_elements(std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos) = 0;
//    virtual void save_npdm_text(const int &i, const int &j) = 0;
//    virtual void save_npdm_binary(const int &i, const int &j) = 0;
//    virtual void save_spatial_npdm_text(const int &i, const int &j) = 0;
//    virtual void save_spatial_npdm_binary(const int &i, const int &j) = 0;
//    virtual void load_npdm_binary(const int &i, const int &j) = 0;
//    virtual void npdm_resize_array(int dim) = 0;
//    virtual void npdm_clear_array() = 0;
//    virtual void accumulate_npdm() = 0;
//
//};
//  
////===========================================================================================================================================================
//
//}
//
//#endif
//
