#ifndef __NPDM_SYMMETRIC_SPATIAL_ARRAY_3D_H
#define __NPDM_SYMMETRIC_SPATIAL_ARRAY_3D_H

#include <vector>
#include <algorithm>

#include <npdm_symmetric_index.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>

namespace SpinAdapted {
namespace Npdm {

/// supportive class to get actual index
struct __SC_index_info_3pdm
{
   unsigned int type_;

   unsigned int offset_;

   size_t ordinal_;

   __SC_index_info_3pdm (int i, int j, int k)
   {
      offset_ = 0u;
      switch(__SF_get_unique_type(i,j,k)) // gives unique number for each case
      {
         case C3_IJK:
            offset_ = 0u; type_ = T3_IJK; ordinal_ = i*(i-1)*(i-2)/6+j*(j-1)/2+k; break;
         case C3_IKJ:
            offset_ = 1u; type_ = T3_IJK; ordinal_ = i*(i-1)*(i-2)/6+k*(k-1)/2+j; break;
         case C3_JIK:
            offset_ = 2u; type_ = T3_IJK; ordinal_ = j*(j-1)*(j-2)/6+i*(i-1)/2+k; break;
         case C3_JKI:
            offset_ = 3u; type_ = T3_IJK; ordinal_ = k*(k-1)*(k-2)/6+i*(i-1)/2+j; break;
         case C3_KIJ:
            offset_ = 4u; type_ = T3_IJK; ordinal_ = j*(j-1)*(j-2)/6+k*(k-1)/2+i; break;
         case C3_KJI:
            offset_ = 5u; type_ = T3_IJK; ordinal_ = k*(k-1)*(k-2)/6+j*(j-1)/2+i; break;

         case C3_IJJ:
            offset_ = 0u; type_ = T3_IJJ; ordinal_ = i*(i-1)/2+j; break;
         case C3_JIJ:
            offset_ = 1u; type_ = T3_IJJ; ordinal_ = j*(j-1)/2+k; break;
         case C3_JJI:
            offset_ = 2u; type_ = T3_IJJ; ordinal_ = k*(k-1)/2+i; break;

         case C3_IIJ:
            offset_ = 0u; type_ = T3_IIJ; ordinal_ = i*(i-1)/2+k; break;
         case C3_IJI:
            offset_ = 1u; type_ = T3_IIJ; ordinal_ = k*(k-1)/2+j; break;
         case C3_JII:
            offset_ = 2u; type_ = T3_IIJ; ordinal_ = j*(j-1)/2+i; break;

         case C3_III:
            offset_ = 0u; type_ = T3_III; ordinal_ = i; break;

         default: abort();
      }
   }

   bool operator<  (const __SC_index_info_3pdm& x) const
   {
      switch (type_)
      {
         case T3_IJK:
            if(x.type_ == T3_IJK)
               return ordinal_ < x.ordinal_;
            else
               return false;
         case T3_IJJ:
         case T3_IIJ:
            if(x.type_ == T3_IJJ || x.type_ == T3_IIJ)
               return ordinal_ < x.ordinal_;
            else
               return type_ < x.type_;
         case T3_III:
            if(x.type_ == T3_III)
               return ordinal_ < x.ordinal_;
            else
               return true;
         default:
            abort();
      }
   }

   void swap (__SC_index_info_3pdm& x)
   {
      std::swap(type_, x.type_);
      std::swap(offset_, x.offset_);
      std::swap(ordinal_, x.ordinal_);
   }
}; // struct __SC_index_info_3pdm

// ==================================================================================================== 
// ==================================================================================================== 
// ==================================================================================================== 

template<>
struct symmetric_spatial_array_base<3>
{
   /// sorting {ijk},{lmn} -> ijk,{lmn}
   /// i > j > k, l > m > n, ijk > lmn
   /// get index offset as
   /// offset = C_offset_ijk_lmn[ijk.case()][lmn.case()];
   static const unsigned int
   C_offset_ijk_lmn[6][6];

   /// sorting {ijk},{ijk} -> ijk,{ijk}
   /// i > j > k : a special case ijk = lmn above
   /// removed conjugate pair, i.e. (ijk,jki) = (ijk,kij)^t
   static const unsigned int
   C_offset_ijk_ijk[6][6];

   static const unsigned int
   C_offset_ijk_lmm[6][3];

   static const unsigned int
   C_offset_ijk_llm[6][3];

// static const unsigned int
// C_offset_iii_lmn[6][1];

   static const unsigned int
   C_offset_ijj_lmm[3][3];

   static const unsigned int
   C_offset_ijj_llm[3][3];

   // note that (ijj,iij) = (iij,ijj)
   static const unsigned int
   C_offset_ijj_iij[3][3];

// static const unsigned int
// C_offset_ijj_lll[3][1];

   static const unsigned int
   C_offset_iij_lmm[3][3];

   static const unsigned int
   C_offset_iij_llm[3][3];

// static const unsigned int
// C_offset_iij_lmm[3][1];

// static const unsigned int
// C_offset_iii_lll[1][1];
};

/// rank-6 array for spatial 3RDM with full-permutation symmetry
/// def. < C_i C_j C_k D_l D_m D_n > : C_p & D_q are resp. creation & destruction operators
template<typename T>
class symmetric_spatial_array<T, 3> : public symmetric_spatial_array_base<3>
{

public:

   explicit
   symmetric_spatial_array (size_t n = 0)
   { this->resize(n); }

   symmetric_spatial_array (const symmetric_spatial_array& x)
   :  extent_          (x.extent_),
      extent_2m_       (x.extent_2m_),
      extent_3m_       (x.extent_3m_),
      V_offset_ijk_lmn_(x.V_offset_ijk_lmn_),
      V_offset_ijk_ijk_(x.V_offset_ijk_ijk_),
      V_offset_ijk_lmm_(x.V_offset_ijk_lmm_),
      V_offset_ijk_llm_(x.V_offset_ijk_llm_),
      V_offset_ijk_lll_(x.V_offset_ijk_lll_),
      V_offset_ijj_lmm_(x.V_offset_ijj_lmm_),
      V_offset_ijj_llm_(x.V_offset_ijj_llm_),
      V_offset_ijj_iij_(x.V_offset_ijj_iij_),
      V_offset_ijj_lll_(x.V_offset_ijj_lll_),
      V_offset_iij_lmm_(x.V_offset_iij_lmm_),
      V_offset_iij_llm_(x.V_offset_iij_llm_),
      V_offset_iij_lll_(x.V_offset_iij_lll_),
      V_offset_iii_lll_(x.V_offset_iii_lll_),
      store_           (x.store_)
   { }

   symmetric_spatial_array& operator= (const symmetric_spatial_array& x)
   {
      extent_           = x.extent_;
      extent_2m_        = x.extent_2m_;
      extent_3m_        = x.extent_3m_;
      V_offset_ijk_lmn_ = x.V_offset_ijk_lmn_;
      V_offset_ijk_ijk_ = x.V_offset_ijk_ijk_;
      V_offset_ijk_lmm_ = x.V_offset_ijk_lmm_;
      V_offset_ijk_llm_ = x.V_offset_ijk_llm_;
      V_offset_ijk_lll_ = x.V_offset_ijk_lll_;
      V_offset_ijj_lmm_ = x.V_offset_ijj_lmm_;
      V_offset_ijj_llm_ = x.V_offset_ijj_llm_;
      V_offset_ijj_iij_ = x.V_offset_ijj_iij_;
      V_offset_ijj_lll_ = x.V_offset_ijj_lll_;
      V_offset_iij_lmm_ = x.V_offset_iij_lmm_;
      V_offset_iij_llm_ = x.V_offset_iij_llm_;
      V_offset_iij_lll_ = x.V_offset_iij_lll_;
      V_offset_iii_lll_ = x.V_offset_iii_lll_;
      store_            = x.store_;
      return *this;
   }

   void resize (size_t n)
   {
      extent_    = n;
      extent_2m_ = extent_*(extent_-1)/2;
      extent_3m_ = extent_*(extent_-1)*(extent_-2)/6;
      store_.resize(__SF_set_storage_size());
      std::fill(store_.begin(), store_.end(), static_cast<T>(0));
// DEBUG */ std::cout << "Spatial 3PDM :: # of non-redundant elements = " << store_.size() << std::endl;
   }

   T& operator() (int i, int j, int k, int l, int m, int n)
   { return store_[__SF_get_ordinal(i,j,k,l,m,n)]; }

   const T& operator() (int i, int j, int k, int l, int m, int n) const
   { return store_[__SF_get_ordinal(i,j,k,l,m,n)]; }

   void fill (const T& value) { std::fill(store_.begin(), store_.end(), value); }

   void clear () { extent_ = 0; store_.swap(std::vector<T>()); }

   const size_t& extent () const
   { return extent_; }

private:

   size_t extent_; ///< extent of each dimension, # of MOs

   size_t extent_2m_; ///< N*(N-1)/2 where N = extent_

   size_t extent_3m_; ///< N*(N-1)*(N-2)/6 where N = extent_

   size_t V_offset_ijk_lmn_;
   size_t V_offset_ijk_ijk_;

   size_t V_offset_ijk_lmm_;
   size_t V_offset_ijk_llm_;
   size_t V_offset_ijk_lll_;

   size_t V_offset_ijj_lmm_;
   size_t V_offset_ijj_llm_;
   size_t V_offset_ijj_iij_;
   size_t V_offset_ijj_lll_;

   size_t V_offset_iij_lmm_;
   size_t V_offset_iij_llm_;
   size_t V_offset_iij_lll_;

   size_t V_offset_iii_lll_;

   std::vector<T> store_;

   // BOOST::SERIALIZATION

   friend class boost::serialization::access;

   template<class Archive>
   void serialize (Archive & ar, const unsigned int version)
   {
      ar & extent_ & extent_2m_ & extent_3m_ &
      V_offset_ijk_lmn_ & V_offset_ijk_ijk_ & V_offset_ijk_lmm_ & V_offset_ijk_llm_ &
      V_offset_ijk_lll_ & V_offset_ijj_lmm_ & V_offset_ijj_llm_ & V_offset_ijj_iij_ &
      V_offset_ijj_lll_ & V_offset_iij_lmm_ & V_offset_iij_llm_ & V_offset_iij_lll_ &
      V_offset_iii_lll_ & store_;
   }

   // SUPPORTIVE CLASS & FUNCTION

   // return total # of unique elements to be stored
   size_t __SF_set_storage_size ()
   {
      V_offset_ijk_lmn_ = 0ul;
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
      return              V_offset_iii_lll_ + extent_*(extent_+1)/2;
   }

   // algorithm to determine non-zero element?
   // given i,j,k,l,m,n find < C_i C_j C_k D_l D_m D_n >
   size_t __SF_get_ordinal (int i, int j, int k, int l, int m, int n) const
   {
      // check {i,j,k}
      __SC_index_info_3pdm ijk(i,j,k);

      // check {l,m,n}
      __SC_index_info_3pdm lmn(n,m,l); // (i,n),(j,m),(k,l)

      if(ijk < lmn) std::swap(ijk,lmn);

      unsigned int ijkOrd = ijk.ordinal_;
      unsigned int ijkOff = ijk.offset_;

      unsigned int lmnOrd = lmn.ordinal_;
      unsigned int lmnOff = lmn.offset_;

      size_t iaddr = 0ul;
      switch ((ijk.type_ << 2) | lmn.type_)
      {
         // N1 = extent_
         // N2 = N1*(N1-1)/2
         // N3 = N1*(N1-1)*(N1-2)/6
         case T3_IJKLMN:
            if(ijkOrd == lmnOrd)
            {
               // w/ 5 permutations in T(ijk,ijk)
               // size: 5*N3
               iaddr = V_offset_ijk_ijk_
                     + 5*(ijkOrd)
                     + C_offset_ijk_ijk[ijkOff][lmnOff];
            }
            else
            {
               // w/ 6 permutations in T(ijk,lmn)
               // size: 6*N3*(N3-1)/2
               iaddr = V_offset_ijk_lmn_
                     + 6*(ijkOrd*(ijkOrd-1)/2+lmnOrd)
                     + C_offset_ijk_lmn[ijkOff][lmnOff];
            }
            break;
         case T3_IJKLMM:
            // size: 3*N3*N2
            iaddr = V_offset_ijk_lmm_
                  + 3*(ijkOrd*extent_2m_+lmnOrd)
                  + C_offset_ijk_lmm[ijkOff][lmnOff];
            break;
         case T3_IJKLLM:
            // size: 3*N3*N2
            iaddr = V_offset_ijk_llm_
                  + 3*(ijkOrd*extent_2m_+lmnOrd)
                  + C_offset_ijk_llm[ijkOff][lmnOff];
            break;
         case T3_IJKLLL:
            // size: N3*N1
            iaddr = V_offset_ijk_lll_
                  + ijkOrd*extent_+lmnOrd;
            break;
         case T3_IJJLMM:
            // size: 2*N2*(N2+1)/2
            iaddr = V_offset_ijj_lmm_
                  + 2*(ijkOrd*(ijkOrd+1)/2+lmnOrd)
                  + C_offset_ijj_lmm[ijkOff][lmnOff];
            break;
         case T3_IJJLLM:
            if(ijkOrd == lmnOrd)
            {
               // size: 2*N2
               iaddr = V_offset_ijj_iij_
                     + 2*ijkOrd
                     + C_offset_ijj_iij[ijkOff][lmnOff];
            }
            else
            {
               // size: 2*N2*(N2-1)/2
               iaddr = V_offset_ijj_llm_
                     + 2*(ijkOrd*(ijkOrd-1)/2+lmnOrd)
                     + C_offset_ijj_llm[ijkOff][lmnOff];
            }
            break;
         case T3_IJJLLL:
            // size: N2*N1
            iaddr = V_offset_ijj_lll_
                  + ijkOrd*extent_+lmnOrd;
            break;
         case T3_IIJLMM:
            if(ijkOrd == lmnOrd)
            {
               // size: 2*N2 -- iij_ijj is the same as ijj_iij (only transposed)
               iaddr = V_offset_ijj_iij_
                     + 2*ijkOrd
                     + C_offset_ijj_iij[lmnOff][ijkOff];
            }
            else
            {
               // size: 2*N2*(N2-1)/2
               iaddr = V_offset_iij_lmm_
                     + 2*(ijkOrd*(ijkOrd-1)/2+lmnOrd)
                     + C_offset_iij_lmm[ijkOff][lmnOff];
            }
            break;
         case T3_IIJLLM:
            // size: 2*N2*(N2+1)/2
            iaddr = V_offset_iij_llm_
                  + 2*(ijkOrd*(ijkOrd+1)/2+lmnOrd)
                  + C_offset_iij_llm[ijkOff][lmnOff];
            break;
         case T3_IIJLLL:
            // size: N2*N1
            iaddr = V_offset_iij_lll_
                  + ijkOrd*extent_+lmnOrd;
            break;
         case T3_IIILLL:
            // size: N1*(N1+1)/2
            iaddr = V_offset_iii_lll_
                  + ijkOrd*(ijkOrd+1)/2+lmnOrd;
      }

// DEBUG */ std::cout << "Spatial-array 3d :: " << i << "," << j << "," << k << "," << l << "," << m << "," << n << " :: " << iaddr << "/" << store_.size() << std::endl;
      return iaddr;
   }

}; // class symmetric_spatial_array<T, 3>

} // namespace Npdm
} // namespace SpinAdapted

#endif // __NPDM_SYMMETRIC_SPATIAL_ARRAY_3D_H
