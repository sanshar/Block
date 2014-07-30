#include <iostream>

#include "npdm_symmetric_index.hpp"

int main ()
{
  using namespace SpinAdapted::Npdm;

  size_t i,j,k;

  while(1) {
    std::cin >> i >> j >> k;
    std::cout << __SF_get_unique_type(i,j,k) << std::endl;
  }

  return 0;
}
