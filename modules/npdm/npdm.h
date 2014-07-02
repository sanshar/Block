/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NPDM_HEADER
#define NPDM_HEADER

namespace SpinAdapted{
namespace Npdm{

  void npdm(int npdm_order);
  void npdm_restart(int npdm_order);
  void transition_pdm( int npdm_order );
  void transition_pdm_restart( int npdm_order );

}
}

#endif
