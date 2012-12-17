/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/



#ifndef PRINT_UTILS_HEADER_H
#define PRINT_UTILS_HEADER_H
#include <iostream>
#include <vector>
using namespace std;
template<class T> void print(ostream& os, const vector<T>& v)
{
  for (int i = 0; i < v.size(); ++i)
    os << v[i] << " ";
  os << endl;
}
#endif
