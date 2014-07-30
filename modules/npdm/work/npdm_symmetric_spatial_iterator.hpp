#ifndef __NPDM_SYMMETRIC_SPATIAL_ITERATOR_HPP
#define __NPDM_SYMMETRIC_SPATIAL_ITERATOR_HPP

namespace SpinAdapted {
namespace Npdm {

template<typename Iter, size_t N> class symmetric_spatial_iterator;

template<typename Iter>
class symmetric_spatial_iterator<Iter, 3>
{

public:

  typedef std::array<int, 6> index_type;

  symmetric_spatial_iterator (const size_t& n, Iter start)
  : extent_(n), It_current_(start)
  {
    // symmetric 2d extent for convenience
    size_t extent_2m_ = extent_*(extent_-1)/2;
    // symmetric 3d extent for convenience
    size_t extent_3m_ = extent_*(extent_-1)*(extent_-2)/6;

    It_offset_ijk_lmn_ = st;
    // T(ijk,lmn), T(ijk,lnm), T(ijk,mln), T(ijk,mnl), T(ijk,nlm), T(ijk,nml) : i > j > k, l > m > n, ijk > lmn
    It_offset_ijk_ijk_ = It_offset_ijk_lmn_ + 6*extent_3m_*(extent_3m_-1)/2;
    // T(ijk,ijk), T(ijk,ikj), T(ijk,jik), T(ijk,jki), T(ijk,kji) : i > j > k
    It_offset_ijk_lmm_ = It_offset_ijk_ijk_ + 5*extent_3m_;
    // T(ijk,lmm), T(ijk,mlm), T(ijk,mml) : i > j > k, l > m
    It_offset_ijk_llm_ = It_offset_ijk_lmm_ + 3*extent_3m_*extent_2m_;
    // T(ijk,llm), T(ijk,lml), T(ijk,mll) : i > j > k, l > m
    It_offset_ijk_lll_ = It_offset_ijk_llm_ + 3*extent_3m_*extent_2m_;
    // T(ijk,lll) : i > j > k
    It_offset_ijj_lmm_ = It_offset_ijk_lll_ + extent_3m_*extent_;
    // T(ijj,lmm), T(ijj,mml) : i > j, l > m, ij >= lm
    It_offset_ijj_llm_ = It_offset_ijj_lmm_ + extent_2m_*(extent_2m_+1);
    // T(ijj,llm), T(ijj,mll) : i > j, l > m, ij >  lm
    It_offset_ijj_iij_ = It_offset_ijj_llm_ + extent_2m_*(extent_2m_-1);
    // T(ijj,iij), T(ijj,jii) : i > j
    It_offset_ijj_lll_ = It_offset_ijj_iij_ + 2*extent_2m_;
    // T(ijj,lll) : i > j
    It_offset_iij_lmm_ = It_offset_ijj_lll_ + extent_2m_*extent_;
    // T(iij,lmm), T(iij,mml) : i > j, l > m, ij >  lm
    It_offset_iij_llm_ = It_offset_iij_lmm_ + extent_2m_*(extent_2m_-1);
    // T(iij,llm), T(iij,lmm) : i > j, l > m, ij >= lm
    It_offset_iij_lll_ = It_offset_iij_llm_ + extent_2m_*(extent_2m_+1);
    // T(iij,lll) : i > j
    It_offset_iii_lll_ = It_offset_iij_lll_ + extent_2m_*extent_;
    // T(iii,lll) : i >= l
    It_offset_end_     = It_offset_iii_lll_ + extent_*(extent_+1)/2;

    incrementum_.reset(new SMC_incrementum_ijk_lmn());
  }

  const index_type& index () const { return index_; }

  // ++it
  symmetric_iterator& operator++ () {
    if(incrementum_->increment()) incrementum_ = incrementum_->next();
  }

private:

  struct SMC_incrementum_base
  {
    SMC_incrementum_base (size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) : ib_{i,j,k,l,m,n} { }

    virtual ~SMC_incrementum_base () { }

    virtual void increment () = 0;

    virtual boost::shared_ptr<SMC_incrementum_base> next () = 0;

    std::array<size_t, 6> ib_;
  };

  struct SMC_incrementum_ijk_lmn : public SMC_incrementum_base
  {
    SMC_incrementum_ijk_lmn () : offset_(0), SMC_incrementum_base(3,1,0,2,1,0) { }

   ~SMC_incrementum_ijk_lmn () { }

    void increment () {
      if(offset_ < 5) {
        ++offset_;
      }
      else {
        size_t i = 5;
        for(; i > 0; --i) {
          if(++ib_[i] < ib_[i-1]) break;
          ib_[i] = 0;
        }
        if(i == 0) {
          ++ib_[0];
        }
        offset_ = 0;
      }

      if(ib_[0] == extent_) return true;

      switch(offset_) {
        case 0:
          index_ = { ib_[0], ib_[1], ib_[2], ib_[3], ib_[4], ib_[5] }; break;
        case 1:
          index_ = { ib_[0], ib_[1], ib_[2], ib_[3], ib_[5], ib_[4] }; break;
        case 2:
          index_ = { ib_[0], ib_[1], ib_[2], ib_[4], ib_[3], ib_[5] }; break;
        case 3:
          index_ = { ib_[0], ib_[1], ib_[2], ib_[4], ib_[5], ib_[3] }; break;
        case 4:
          index_ = { ib_[0], ib_[1], ib_[2], ib_[5], ib_[3], ib_[4] }; break;
        case 5:
          index_ = { ib_[0], ib_[1], ib_[2], ib_[5], ib_[4], ib_[3] }; break;
        default:
          abort();
      }

      return false;
    }

    boost::shared_ptr<SMC_incrementum_base> next () {
      return boost::shared_ptr<SMC_incrementum_base>(new SMC_incrementum_ijk_ijk());
    }

    // member

    size_t offset_;
  };

  struct SMC_incrementum_ijk_ijk : public SMC_incrementum_base
  {
    SMC_incrementum_ijk_ijk () : offset_(0), SMC_incrementum_base() { }
   ~SMC_incrementum_ijk_ijk () { }

    void increment () {
      if(offset_ < 4) {
        ++offset_;
      }
      else {
        size_t i = 2;
        for(; i > 0; --i) {
          if(++ib_[i] < ib_[i-1]) break;
          ib_[i] = 0;
        }
        if(i == 0) {
          ++ib_[0];
        }
        offset_ = 0;
      }

      if(ib_[0] == extent_) return true;

      switch(offset_) {
        case 0:
          index_ = { ib_[0], ib_[1], ib_[2], ib_[0], ib_[1], ib_[2] }; break;
        case 1:
          index_ = { ib_[0], ib_[1], ib_[2], ib_[0], ib_[2], ib_[1] }; break;
        case 2:
          index_ = { ib_[0], ib_[1], ib_[2], ib_[1], ib_[0], ib_[2] }; break;
        case 3:
          index_ = { ib_[0], ib_[1], ib_[2], ib_[1], ib_[2], ib_[0] }; break;
        case 4:
          index_ = { ib_[0], ib_[1], ib_[2], ib_[2], ib_[1], ib_[0] }; break;
        default:
          abort();
      }

      return false;
    }

    boost::shared_ptr<SMC_incrementum_base> next () {
      return boost::shared_ptr<SMC_incrementum_base>(new SMC_incrementum_ijk_lmm());
    }

    // member

    size_t offset_;
  };

  boost::shared_ptr<SMC_incrementum_base> incrementum_;


/*    V_offset_ijk_lmn_ = 0ul;
      // T(ijk,lmn), T(ijk,lnm), T(ijk,mln), T(ijk,mnl), T(ijk,nlm), T(ijk,nml) : i > j > k, l > m > n, ijk > lmn
      V_offset_ijk_ijk_ = V_offset_ijk_lmn_ + 6*extent_3m_*(extent_3m_-1)/2;
      // T(ijk,ijk), T(ijk,ikj), T(ijk,jik), T(ijk,jki), T(ijk,kji) : i > j > k
      V_offset_ijk_lmm_ = V_offset_ijk_ijk_ + 5*extent_3m_;
      // T(ijk,lmm), T(ijk,mlm), T(ijk,mml) : i > j > k, l > m
      V_offset_ijk_llm_ = V_offset_ijk_lmm_ + 3*extent_3m_*extent_2m_;
      // T(ijk,llm), T(ijk,lml), T(ijk,mll) : i > j > k, l > m
      V_offset_ijk_lll_ = V_offset_ijk_llm_ + 3*extent_3m_*extent_2m_;
      // T(ijk,lll) : i > j > k
      V_offset_ijj_lmm_ = V_offset_ijk_lll_ + extent_3m_*extent_;
      // T(ijj,lmm), T(ijj,mml) : i > j, l > m, ij >= lm
      V_offset_ijj_llm_ = V_offset_ijj_lmm_ + extent_2m_*(extent_2m_+1);
      // T(ijj,llm), T(ijj,mll) : i > j, l > m, ij >  lm
      V_offset_ijj_iij_ = V_offset_ijj_llm_ + extent_2m_*(extent_2m_-1);
      // T(ijj,iij), T(ijj,jii) : i > j
      V_offset_ijj_lll_ = V_offset_ijj_iij_ + 2*extent_2m_;
      // T(ijj,lll) : i > j
      V_offset_iij_lmm_ = V_offset_ijj_lll_ + extent_2m_*extent_;
      // T(iij,lmm), T(iij,mml) : i > j, l > m, ij >  lm
      V_offset_iij_llm_ = V_offset_iij_lmm_ + extent_2m_*(extent_2m_-1);
      // T(iij,llm), T(iij,lmm) : i > j, l > m, ij >= lm
      V_offset_iij_lll_ = V_offset_iij_llm_ + extent_2m_*(extent_2m_+1);
      // T(iij,lll) : i > j
      V_offset_iii_lll_ = V_offset_iij_lll_ + extent_2m_*extent_;
      // T(iii,lll) : i >= l
      return              V_offset_iii_lll_ + extent_*(extent_+1)/2; */

  Iter It_offset_ijk_lmn_;
  Iter It_offset_ijk_ijk_;

  Iter It_offset_ijk_lmm_;
  Iter It_offset_ijk_llm_;
  Iter It_offset_ijk_lll_;

  Iter It_offset_ijj_lmm_;
  Iter It_offset_ijj_llm_;
  Iter It_offset_ijj_iij_;
  Iter It_offset_ijj_lll_;

  Iter It_offset_iij_lmm_;
  Iter It_offset_iij_llm_;
  Iter It_offset_iij_lll_;

  Iter It_offset_iii_lll_;

  Iter It_offset_end_;

  Iter It_current_;

  size_t extent_;

  index_type index_;

};

} // namespace Npdm
} // namespace SpinAdapted

#endif // __NPDM_SYMMETRIC_SPATIAL_ITERATOR_HPP
