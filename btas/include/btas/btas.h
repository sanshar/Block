#ifndef _BTAS_CXX11_H
#define _BTAS_CXX11_H 1

#include <vector>
#include <cassert>

// TODO:
// C++11 involves these classes, however, they're not supported by boost serialization.
// They should be replaced by STL library in the future.
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>

#include <iostream>
#define BTAS_DEBUG(msg)\
{ std::cout << "BTAS_DEBUG: " << msg << std::endl; }

#include <stdexcept>
#define BTAS_THROW(truth, msg)\
{ if (!(truth)) { throw std::runtime_error(msg); } }

namespace btas {

typedef unsigned int  uint;
typedef unsigned long wint;

//! Alias to dense shape
typedef std::vector<int> Dshapes;

// Enables scope of boost function and smart pointer
using boost::shared_ptr;
using boost::weak_ptr;
using boost::function;
using boost::bind;

// For full C++11 support
// using std::shared_ptr;
// using std::weak_ptr;
// using std::function;
// using std::bind;
// using std::placeholders;

//! Null-deleter
struct null_deleter {
  void operator() (void const *) const { }
};

}; // namespace btas

#endif // _BTAS_CXX11_H
