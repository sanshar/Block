#ifndef __NPDM_SYMMETRIC_ITERATOR_HPP
#define __NPDM_SYMMETRIC_ITERATOR_HPP

namespace SpinAdapted {
namespace Npdm {

template<typename Iter, size_t N> class symmetric_iterator;

template<typename Iter>
class symmetric_iterator<Iter, 4> : public std::forward_iterator_tag
{

public:

  typedef std::array<int, 2*N> index_type;

  const index_type& index () const { return index_; }


  // ++it
  symmetric_iterator& operator++ ()
  {
    if(r_index_ < l_index_)
    {
      ++l_index_;
      size_t i = N-1; for(; i > 0; --i)
      {
        if(++index_[N+i] < index_[N+i-1]) break;
        index_[N+i] = N-i-1;
      }
      if(i == 0) ++index_[N];
    }
    else
    {
      l_index_ = 0; ++r_index_;
      size_t i = N-1; for(; i > 0; --i)
      {
        if(++index_[i] < index_[i-1]) break;
        index_[i] = N-i-1;
      }
      if(i == 0) ++index_[0];
    }
  }

private:

  Iter curr_;

  size_t extent_;

  size_t l_index_;

  size_t r_index_;

  index_type index_;

};

} // namespace Npdm
} // namespace SpinAdapted

#endif // __NPDM_SYMMETRIC_ITERATOR_HPP
