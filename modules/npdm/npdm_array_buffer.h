/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NPDM_SPARSE_ARRAY_HEADER_H
#define NPDM_SPARSE_ARRAY_HEADER_H

#include <map>
#include <vector>

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

class Npdm_array_buffer {

  public:
    Npdm_array_buffer( const int num_indices, const int dim ) : bufferID_(0), dim_(dim) {};
    ~Npdm_array_buffer() {};

    double operator()(int i, int j, int k, int l, int m, int n) const;
    double& operator()(int i, int j, int k, int l, int m, int n);
    void close_buffer();
    int dim1() const { return dim_; }

  private:
    int bufferID_;
    int dim_;
    std::map< std::vector<int>, double > data_;

};
  
//===========================================================================================================================================================

}
}

#endif

