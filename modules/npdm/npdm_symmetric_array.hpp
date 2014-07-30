#ifndef __NPDM_SYMMETRIC_ARRAY_HPP
#define __NPDM_SYMMETRIC_ARRAY_HPP

#include <algorithm>
#include <npdm_symmetric_index.hpp>
#include <npdm_symmetric_element.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>

namespace SpinAdapted {
namespace Npdm {

template<typename T, size_t N> class symmetric_array;

/// Array class for 1PDM
template<typename T>
class symmetric_array<T, 1>
{

   typedef T& element_type;

   typedef const T& const_element_type;

public:

   explicit
   symmetric_array (size_t n = 0) { this->resize(n); }

   element_type operator() (int i, int j)
   {
      return store_[get_ordinal_index(i,j)];
   }

   const_element_type operator() (int i, int j) const
   {
      return store_[get_ordinal_index(i,j)];
   }

   void resize (size_t n)
   {
      extent_ = n;
      store_.clear();
      store_.resize(n*(n+1)/2, static_cast<T>(0));
   }

   void fill (const T& value) { std::fill(store_.begin(), store_.end(), value); }

   void clear () { extent_ = 0; store_.swap(std::vector<T>()); }

   size_t extent () const { return extent_; }

private:

   friend class boost::serialization::access;

   template<class Archive>
   void serialize (Archive & ar, const unsigned int version)
   {
      ar & extent_ & store_;
   }

   inline size_t get_ordinal_index(int i, int j) const
   {
      return (i > j) ? i*(i+1)/2+j : j*(j+1)/2+i;
   }

   /// extent of array
   size_t extent_;

   /// storage for non-redundant elements
   std::vector<T> store_;

};

/// Array class for 2PDM
template<typename T>
class symmetric_array<T, 2>
{

   typedef symmetric_element<T> element_type;

   typedef const_symmetric_element<T> const_element_type;

public:

   explicit
   symmetric_array (size_t n = 0) { this->resize(n); }

   element_type operator() (int i, int j, int k, int l)
   {
      return element_type(store_[get_ordinal_index(i,j,k,l)], get_parity(i,j)^get_parity(k,l));
   }

   const_element_type operator() (int i, int j, int k, int l) const
   {
      return const_element_type(store_[get_ordinal_index(i,j,k,l)], get_parity(i,j)^get_parity(k,l));
   }

   void resize (size_t n)
   {
      extent_ = n;
      // FIXME: some zero elements are involved
      //        for really non-zero elements, use n*(n-1)/2 instead
      //        but, this makes additional if statement at element access above...
      size_t n2 = n*(n+1)/2;
      store_.clear();
      store_.resize(n2*(n2+1)/2, static_cast<T>(0));
   }

   void fill (const T& value) { std::fill(store_.begin(), store_.end(), value); }

   void clear () { extent_ = 0; store_.swap(std::vector<T>()); }

   size_t extent () const { return extent_; }

private:

   friend class boost::serialization::access;

   template<class Archive>
   void serialize (Archive & ar, const unsigned int version)
   {
      ar & extent_ & store_;
   }

   inline size_t get_ordinal_index_half(int i, int j) const
   {
      return (i > j) ? i*(i+1)/2+j : j*(j+1)/2+i;
   }

   inline size_t get_ordinal_index(int i, int j, int k, int l) const
   {
      int ij = get_ordinal_index_half(i,j);
      int kl = get_ordinal_index_half(k,l);

      return (ij > kl) ? ij*(ij+1)/2+kl : kl*(kl+1)/2+ij;
   }

   /// extent of array
   size_t extent_;

   /// storage for non-redundant elements
   std::vector<T> store_;

};

/// Array class for 3PDM
template<typename T>
class symmetric_array<T, 3>
{

   typedef symmetric_element<T> element_type;

   typedef const_symmetric_element<T> const_element_type;

public:

   explicit
   symmetric_array (size_t n = 0) { this->resize(n); }

   element_type operator() (int i, int j, int k, int l, int m, int n)
   {
      return element_type(store_[get_ordinal_index(i,j,k,l,m,n)], get_parity(i,j,k)^get_parity(l,m,n));
   }

   const_element_type operator() (int i, int j, int k, int l, int m, int n) const
   {
      return const_element_type(store_[get_ordinal_index(i,j,k,l,m,n)], get_parity(i,j,k)^get_parity(l,m,n));
   }

   void resize (size_t n)
   {
      extent_ = n;
      // FIXME: some zero elements are involved
      //        for really non-zero elements, use n*(n-1)*(n-2)/6 instead
      //        but, this makes additional if statement at element access above...
      size_t n3 = n*(n+1)*(n+2)/6;
      store_.clear();
      store_.resize(n3*(n3+1)/2, static_cast<T>(0));
   }

   void fill (const T& value) { std::fill(store_.begin(), store_.end(), value); }

   void clear () { extent_ = 0; store_.swap(std::vector<T>()); }

   size_t extent () const { return extent_; }

private:

   friend class boost::serialization::access;

   template<class Archive>
   void serialize (Archive & ar, const unsigned int version)
   {
      ar & extent_ & store_;
   }

   inline size_t get_ordinal_index_half(int i, int j, int k) const
   {
      int idx[3] = { i, j, k }; std::sort(idx, idx+3);

      return idx[2]*(idx[2]+1)*(idx[2]+2)/6+idx[1]*(idx[1]+1)/2+idx[0];
   }

   inline size_t get_ordinal_index(int i, int j, int k, int l, int m, int n) const
   {
      int ijk = get_ordinal_index_half(i,j,k);
      int lmn = get_ordinal_index_half(l,m,n);

      return (ijk > lmn) ? ijk*(ijk+1)/2+lmn : lmn*(lmn+1)/2+ijk;
   }

   /// extent of array
   size_t extent_;

   /// storage for non-redundant elements
   std::vector<T> store_;

};

/// Array class for 4PDM
template<typename T>
class symmetric_array<T, 4>
{

   typedef symmetric_element<T> element_type;

   typedef const_symmetric_element<T> const_element_type;

public:

   explicit
   symmetric_array (size_t n = 0) { this->resize(n); }

   element_type operator() (int i, int j, int k, int l, int m, int n, int p, int q)
   {
      return element_type(store_[get_ordinal_index(i,j,k,l,m,n,p,q)], get_parity(i,j,k,l)^get_parity(m,n,p,q));
   }

   const_element_type operator() (int i, int j, int k, int l, int m, int n, int p, int q) const
   {
      return const_element_type(store_[get_ordinal_index(i,j,k,l,m,n,p,q)], get_parity(i,j,k,l)^get_parity(m,n,p,q));
   }

   void resize (size_t n)
   {
      extent_ = n;
      // FIXME: some zero elements are involved
      //        for really non-zero elements, use n*(n-1)*(n-2)*(n-3)/24 instead
      //        but, this makes additional if statement at element access above...
      size_t n4 = n*(n+1)*(n+2)*(n+3)/24;
      store_.clear();
      store_.resize(n4*(n4+1)/2, static_cast<T>(0));
   }

   void fill (const T& value) { std::fill(store_.begin(), store_.end(), value); }

   void clear () { extent_ = 0; store_.swap(std::vector<T>()); }

   size_t extent () const { return extent_; }

private:

   friend class boost::serialization::access;

   template<class Archive>
   void serialize (Archive & ar, const unsigned int version)
   {
      ar & extent_ & store_;
   }

   inline size_t get_ordinal_index_half(int i, int j, int k, int l) const
   {
      int idx[4] = { i, j, k, l }; std::sort(idx, idx+4);

      return idx[3]*(idx[3]+1)*(idx[3]+2)*(idx[3]+3)/24+idx[2]*(idx[2]+1)*(idx[2]+2)/6+idx[1]*(idx[1]+1)/2+idx[0];
   }

   inline size_t get_ordinal_index(int i, int j, int k, int l, int m, int n, int p, int q) const
   {
      int ijkl = get_ordinal_index_half(i,j,k,l);
      int mnpq = get_ordinal_index_half(m,n,p,q);

      return (ijkl > mnpq) ? ijkl*(ijkl+1)/2+mnpq : mnpq*(mnpq+1)/2+ijkl;
   }

   /// extent of array
   size_t extent_;

   /// storage for non-redundant elements
   std::vector<T> store_;

};

} // namespace Npdm
} // namespace SpinAdapted

#endif // __NPDM_SYMMETRIC_ARRAY_HPP
