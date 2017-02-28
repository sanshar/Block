/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NPDM_HEADER
#define NPDM_HEADER

namespace SpinAdapted{
enum NpdmOrder{NPDM_NEVPT2, NPDM_ONEPDM, NPDM_TWOPDM, NPDM_THREEPDM, NPDM_FOURPDM, NPDM_PAIRMATRIX, NPDM_OVERLAP, NPDM_EMPTY, NPDM_DS0, NPDM_DS1};
namespace Npdm{

  void npdm(NpdmOrder npdm_order, bool restartpdm=false, bool transitionpdm=false, bool dS = false);

}
}

#endif
