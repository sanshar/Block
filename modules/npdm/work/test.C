#include <iostream>
#include <iomanip>

#include "npdm_symmetric_array.hpp"
#include "npdm_symmetric_spatial_array.hpp"

int main ()
{
  using namespace SpinAdapted::Npdm;

  const int& N = 4;

  symmetric_array<int, 3> thrpdm(N*2); thrpdm.fill(1);

  symmetric_spatial_array<int, 3> spatial_thrpdm(N); spatial_thrpdm.fill(0);

//for(int i = 0; i < N; ++i)
//  for(int j = 0; j < N; ++j)
//    for(int k = 0; k < N; ++k)
//      for(int l = 0; l < N; ++l)
//        for(int m = 0; m < N; ++m)
//          for(int n = 0; n < N; ++n) {
//            int i2 = i*2;
//            int j2 = j*2;
//            int k2 = k*2;
//            int l2 = l*2;
//            int m2 = m*2;
//            int n2 = n*2;
//            int value = thrpdm(i2  ,j2  ,k2  ,l2  ,m2  ,n2  )
//                      + thrpdm(i2+1,j2  ,k2  ,l2  ,m2  ,n2+1)
//                      + thrpdm(i2  ,j2+1,k2  ,l2  ,m2+1,n2  )
//                      + thrpdm(i2  ,j2  ,k2+1,l2+1,m2  ,n2  )
//                      + thrpdm(i2  ,j2+1,k2+1,l2+1,m2+1,n2  )
//                      + thrpdm(i2+1,j2  ,k2+1,l2+1,m2  ,n2+1)
//                      + thrpdm(i2+1,j2+1,k2  ,l2  ,m2+1,n2+1)
//                      + thrpdm(i2+1,j2+1,k2+1,l2+1,m2+1,n2+1);
//            if(value != 0)
//              std::cout << i << "," << j << "," << k << "," << l << "," << m << "," << n << " :: " << value << std::endl;

//            spatial_thrpdm(i,j,k,l,m,n) = value;
//          }

//for(int i = 0; i < N; ++i)
//  for(int j = 0; j < N; ++j)
//    for(int k = 0; k < N; ++k)
//      for(int l = 0; l < N; ++l)
//        for(int m = 0; m < N; ++m)
//          for(int n = 0; n < N; ++n) {
//            int value = spatial_thrpdm(i,j,k,l,m,n);
//            if(value != 0)
//              std::cout << i << "," << j << "," << k << "," << l << "," << m << "," << n << " :: " << value << std::endl;
//          }

  for(int i = 0; i < N; ++i)
    for(int j = 0; j < N; ++j)
      for(int k = 0; k < N; ++k)
        for(int l = 0; l < N; ++l)
          for(int m = 0; m < N; ++m)
            for(int n = 0; n < N; ++n) {
              int& value = spatial_thrpdm(i,j,k,l,m,n);
              if(value == 0) value = n+10*m+100*l+1000*k+10000*j+100000*i;
              std::cout << i << "," << j << "," << k << "," << l << "," << m << "," << n << " :: "
                        << std::setfill('0') << std::setw(6) << value << std::endl;
            }

  return 0;
}
