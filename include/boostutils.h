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
