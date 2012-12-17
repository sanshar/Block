/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef BOOSTUTILS_HEADER_H
#define BOOSTUTILS_HEADER_H
namespace boostutils
{
  struct null_deleter
  {
    void operator()(void const*) {}
  };
}
#endif
