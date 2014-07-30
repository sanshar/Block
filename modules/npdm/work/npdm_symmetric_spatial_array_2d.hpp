#ifndef __NPDM_SYMMETRIC_SPATIAL_ARRAY_2D_HPP
#define __NPDM_SYMMETRIC_SPATIAL_ARRAY_2D_HPP

#include <vector>
#include <algorithm>

#include <npdm_symmetric_index.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>

namespace SpinAdapted {
namespace Npdm {

/// supportive class to get actual index
struct __SC_index_info_2pdm
{
   unsigned int type_;

   unsigned int offset_;

   size_t ordinal_;

   __SC_index_info_2pdm (int i, int j)
   {
      offset_ = 0u;
      switch(__SF_get_unique_type(i,j)) // gives unique number for each case
      {
         case C2_IJ:
            offset_ = 0u; type_ = T2_IJ; ordinal_ = i*(i-1)/2+j; break;
         case C2_JI:
            offset_ = 1u; type_ = T2_IJ; ordinal_ = j*(j-1)/2+i; break;
         case C2_II:
            offset_ = 0u; type_ = T2_II; ordinal_ = i; break;

         default:
            abort();
      }
   }

   bool operator<  (const __SC_index_info_2pdm& x) const
   {
      switch (type_)
      {
         case T2_IJ:
            if(x.type_ == T2_IJ)
               return ordinal_ < x.ordinal_;
            else
               return false;
         case T2_II:
            if(x.type_ == T2_II)
               return ordinal_ < x.ordinal_;
            else
               return true;
         default:
            abort();
      }
   }

   void swap (__SC_index_info_2pdm& x)
   {
      std::swap(type_, x.type_);
      std::swap(offset_, x.offset_);
      std::swap(ordinal_, x.ordinal_);
   }
}; // struct __SC_index_info_2pdm

// ==================================================================================================== 
// ==================================================================================================== 
// ==================================================================================================== 

/// rank-6 array for spatial 3RDM with full-permutation symmetry
/// def. < C_i C_j C_k D_l D_m D_n > : C_p & D_q are resp. creation & destruction operators
template<typename T>
class symmetric_spatial_array<T, 2>
{

public:

   explicit
   symmetric_spatial_array (size_t n = 0)
   { this->resize(n); }

   symmetric_spatial_array (const symmetric_spatial_array& x)
   :  extent_        (x.extent_),
      extent_2m_     (x.extent_2m_),
      V_offset_ij_kl_(x.V_offset_ij_kl_),
      V_offset_ij_kk_(x.V_offset_ij_kk_),
      V_offset_ii_kk_(x.V_offset_ii_kk_),
      store_         (x.store_)
   { }

   symmetric_spatial_array& operator= (const symmetric_spatial_array& x)
   {
      extent_         = x.extent_;
      extent_2m_      = x.extent_2m_;
      V_offset_ij_kl_ = x.V_offset_ij_kl_;
      V_offset_ij_kk_ = x.V_offset_ij_kk_;
      V_offset_ii_kk_ = x.V_offset_ii_kk_;
      store_          = x.store_;
      return *this;
   }

   void resize (size_t n)
   {
      extent_    = n;
      extent_2m_ = extent_*(extent_-1)/2;
      store_.resize(__SF_set_storage_size());
      std::fill(store_.begin(), store_.end(), static_cast<T>(0));
   }

   T& operator() (int i, int j, int k, int l)
   { return store_[__SF_get_ordinal(i,j,k,l)]; }

   const T& operator() (int i, int j, int k, int l) const
   { return store_[__SF_get_ordinal(i,j,k,l)]; }

   void fill (const T& value) { std::fill(store_.begin(), store_.end(), value); }

   void clear () { extent_ = 0; store_.swap(std::vector<T>()); }

   const size_t& extent () const
   { return extent_; }

private:

   size_t extent_; ///< extent of each dimension, # of MOs

   size_t extent_2m_; ///< N*(N-1)/2 where N = extent_

   size_t V_offset_ij_kl_;
   size_t V_offset_ij_kk_;
   size_t V_offset_ii_kk_;

   std::vector<T> store_;

   // BOOST::SERIALIZATION

   friend class boost::serialization::access;

   template<class Archive>
   void serialize (Archive & ar, const unsigned int version)
   {
      ar & extent_ & extent_2m_ & V_offset_ij_kl_ & V_offset_ij_kk_ & V_offset_ii_kk_;
   }

   // SUPPORTIVE CLASS & FUNCTION

   // return total # of unique elements to be stored
   size_t __SF_set_storage_size ()
   {
      V_offset_ij_kl_ = 0ul;
      // T(ij,kl), T(ij,lk) : i > j, k > l, ij >= kl
      V_offset_ij_kk_ = V_offset_ij_kl_ + extent_2m_*(extent_2m_+1);
      // T(ij,kk) : i > j
      V_offset_ii_kk_ = V_offset_ij_kk_ + extent_2m_*extent_;
      // T(ii,kk) : i > k
      return            V_offset_ii_kk_ + extent_*(extent_+1)/2;
   }

   // algorithm to determine non-zero element?
   // given i,j,k,l,m,n find < C_i C_j C_k D_l D_m D_n >
   size_t __SF_get_ordinal (int i, int j, int k, int l) const
   {
      // check {i,j,k}
      __SC_index_info_2pdm ij(i,j);

      // check {l,m,n}
      __SC_index_info_2pdm kl(l,k);

      if(ij < kl) std::swap(ij,kl);

      unsigned int ijOrd = ij.ordinal_;
      unsigned int ijOff = ij.offset_;

      unsigned int klOrd = kl.ordinal_;
      unsigned int klOff = kl.offset_;

      size_t iaddr = 0ul;
      switch ((ij.type_ << 1) | kl.type_)
      {
         case T2_IJKL:
            iaddr = V_offset_ij_kl_+2*(ijOrd*(ijOrd+1)/2+klOrd)+(ijOff^klOff);
            break;
         case T2_IJKK:
            iaddr = V_offset_ij_kk_+ijOrd*extent_+klOrd;
            break;
         case T2_IIKK:
            iaddr = V_offset_ii_kk_+ijOrd*(ijOrd+1)/2+klOrd;
      }

      return iaddr;
   }

}; // class symmetric_spatial_array<T, 3>

} // namespace Npdm
} // namespace SpinAdapted

#endif // __NPDM_SYMMETRIC_SPATIAL_ARRAY_2D_HPP
