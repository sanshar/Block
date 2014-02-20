//
/*! \file TSubArray.h
 *  \brief Dense sub-array class and its copy semantics
 */

#ifndef _BTAS_CXX11_TARRAY_H
#include <btas/DENSE/TArray.h>
#endif

#ifndef _BTAS_CXX11_TSUBARRAY_H
#define _BTAS_CXX11_TSUBARRAY_H 1

#include <functional>
#include <numeric>

namespace btas {

//####################################################################################################
// Sub-array class for dense array
//####################################################################################################

//! Sub-array of dense array TArray
/*! Expression template providing copy semantics between TArray and TSubArray
 *  \param T value type
 *  \param N array rank */
template<typename T, size_t N>
class TSubArray : protected TArray<T, N> {
private:
  //! Any TArray classes being friend of TSubArray<T, N>
  template<typename U, size_t M>
  friend class TArray;

  //! Not default-constructible
  TSubArray();

public:
  //! Constructor, from dense array object and lower/upper boundary indices
  TSubArray
  (const TArray<T, N>& a, const IVector<N>& lbound, const IVector<N>& ubound) {
    this->reference(a);
    m_lower_bound = lbound;
    m_upper_bound = ubound;
  }

  //! Copy from array to sub-array
  template<size_t M>
  void copy(const TArray<T, M>& a) {
    // Calc. sub-array shape
    IVector<N> t_shape;
    for(int i = 0; i < N; ++i) t_shape[i] = m_upper_bound[i]-m_lower_bound[i]+1;
    int t_size = std::accumulate(t_shape.begin(), t_shape.end(), 1, std::multiplies<int>());
    assert(a.size() == t_size);
    // If 0-dim. array
    if(t_size == 0) return;
    // Striding
    const IVector<N>& t_stride = this->m_stride;
    int lda = t_shape[N-1];
    // Get bare pointers
    const T* a_ptr = a.data();
          T* t_ptr = this->data();
    // Copying elements
    IVector<N> index(m_lower_bound);
    int nrows = t_size / lda;
    for(int j = 0; j < nrows; ++j, a_ptr += lda) {
      int offset = dot(t_stride, index);
      _fast_copy(lda, a_ptr, t_ptr+offset);
      for(int i = static_cast<int>(N)-2; i >= 0; --i) {
        if(++index[i] <= m_upper_bound[i]) break;
        index[i] = m_lower_bound[i];
      }
    }
  }

  //! Copy assignment operator from TArray
  template<size_t M>
  void operator= (const TArray<T, M>& other) { copy(other); }

private:
  //! lower boundary indices
  IVector<N>
    m_lower_bound;
  //! upper boundary indices
  IVector<N>
    m_upper_bound;
};

}; // namespace btas

#endif // _BTAS_CXX11_TSUBARRAY_H
