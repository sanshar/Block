#ifndef _BTAS_DREINDEX_H
#define _BTAS_DREINDEX_H 1

#include <btas/btas.h>
#include <btas/TVector.h>

namespace btas {

template<size_t N>
void Dreindex(const double* x, double* y, const IVector<N>& xstr, const IVector<N>& yshape) {
  BTAS_THROW(false, "btas::Dreindex does not support rank > 12");
}

template<> void Dreindex< 1>(const double* x, double* y, const IVector< 1>& xstr, const IVector< 1>& yshape);
template<> void Dreindex< 2>(const double* x, double* y, const IVector< 2>& xstr, const IVector< 2>& yshape);
template<> void Dreindex< 3>(const double* x, double* y, const IVector< 3>& xstr, const IVector< 3>& yshape);
template<> void Dreindex< 4>(const double* x, double* y, const IVector< 4>& xstr, const IVector< 4>& yshape);
template<> void Dreindex< 5>(const double* x, double* y, const IVector< 5>& xstr, const IVector< 5>& yshape);
template<> void Dreindex< 6>(const double* x, double* y, const IVector< 6>& xstr, const IVector< 6>& yshape);
template<> void Dreindex< 7>(const double* x, double* y, const IVector< 7>& xstr, const IVector< 7>& yshape);
template<> void Dreindex< 8>(const double* x, double* y, const IVector< 8>& xstr, const IVector< 8>& yshape);
template<> void Dreindex< 9>(const double* x, double* y, const IVector< 9>& xstr, const IVector< 9>& yshape);
template<> void Dreindex<10>(const double* x, double* y, const IVector<10>& xstr, const IVector<10>& yshape);
template<> void Dreindex<11>(const double* x, double* y, const IVector<11>& xstr, const IVector<11>& yshape);
template<> void Dreindex<12>(const double* x, double* y, const IVector<12>& xstr, const IVector<12>& yshape);

}; // namespace btas

#endif // _BTAS_DREINDEX_H
