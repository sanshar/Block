#ifndef _BTAS_CXX11_INDEXED_LOOP_H
#ifndef _BTAS_CXX11_INDEXED_LOOP_H 1

#include <btas/TVector.h>

namespace btas {

//! Indexed loop for array object.
/*! \par Purpose:
 *  Let multiple for loop, e.g.
 *
 *  	for(int i = 0; i < Ni; ++i) {
 *  	  for(int j = 0; j < Nj; ++j) {
 *  	    for(int k = 0; k < Nk; ++k) {
 *  	      ... do something ...
 *  	    }
 *  	  }
 *  	}
 *
 *  be a single loop, i.e.
 *
 *  	for(indexed_loop<3> loop(make_array(0,0,0), make_array(Ni,Nj,Nk)); !loop.end(); ++loop) {
 *  	  IVector<3> index = loop.index();
 *  	  int i = index[0];
 *  	  int j = index[1];
 *  	  int k = index[2];
 *  	  ... do something ...
 *  	}
 */
template<size_t N>
class indexed_loop {
private:
  friend class indexed_loop<N+1>;

  template<size_t M>
  void make_index(IVector<M>& index) const {
    m_oloop.make_index(index);
    index[N-1] = m_index;
  }

public:
  template<size_t M>
  indexed_loop(const IVector<M>& first, const IVector<M>& last)
  : m_oloop(first, last), m_index(first[N-1]), m_first(first[N-1]), m_last(last[N-1]) { }

  template<size_t M>
  void reset(const IVector<M>& first, const IVector<M>& last) {
    m_oloop.reset(first, last);
    m_index = first[N-1];
    m_first = first[N-1];
    m_last  = last [N-1];
  }

  IVector<N> index() const {
    IVector<N> _index;
    make_index(_index);
    return std::move(_index);
  }

  indexed_loop& operator++ () {
    if(++m_index == m_last) {
      m_index = m_first;
      ++m_oloop;
    }
    return *this;
  }

  bool end() const { return m_oloop.end(); }

private:
  indexed_loop<N-1>
    m_oloop; // outer loop control
  int
    m_index;
  int
    m_first;
  int
    m_last;
};

template<>
class indexed_loop<1> {
private:
  friend class indexed_loop<2>;

  template<size_t M>
  void make_index(IVector<M>& index) const { index[0] = m_index; }

public:
  template<size_t M>
  indexed_loop(const IVector<M>& first, const IVector<M>& last)
  : m_index(first[0]), m_first(first[0]), m_last(last[0]), m_end(false) { }

  template<size_t M>
  void reset(const IVector<M>& first, const IVector<M>& last) {
    m_index = first[0];
    m_first = first[0];
    m_last  = last [0];
    m_end   = false;
  }

  IVector<1> index() const {
    IVector<1> _index = { m_index };
    return std::move(_index);
  }

  indexed_loop& operator++ () {
    if(++m_index == m_last) {
      m_index = m_first;
      m_end   = true;
    }
    return *this;
  }

  bool end() const { return m_end; }

private:
  int
    m_index;
  int
    m_first;
  int
    m_last;
  bool
    m_end;
};

}; // namespace btas

#endif // _BTAS_CXX11_INDEXED_LOOP_H
