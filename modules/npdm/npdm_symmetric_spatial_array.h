#ifndef __NPDM_SYMMETRIC_SPATIAL_ARRAY_H
#define __NPDM_SYMMETRIC_SPATIAL_ARRAY_H

#include <algorithm>
#include <npdm_symmetric_index.hpp>
#include <npdm_symmetric_array.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>

namespace SpinAdapted {
namespace Npdm {

template<size_t N> struct symmetric_spatial_array_base;

template<typename T, size_t N> class symmetric_spatial_array;

/// Array class for 1PDM w/ spatial orbitals :: Xij = Xji
template<typename T>
class symmetric_spatial_array<T, 1>
{

public:

   explicit
   symmetric_spatial_array (size_t nMO = 0) { this->resize(nMO); }

   /// C'tor from spin-orbital 1PDM
   explicit
   symmetric_spatial_array (const symmetric_array<T, 1>& spin_array)
   {
      size_t nSO = spin_array.extent();
      size_t nMO = nSO/2;

      this->resize(nMO);

      for(size_t i = 0; i < nMO; ++i)
         for(size_t j = 0; j <= i; ++j)
         {
            for(size_t s = 0; s < 2; ++s)
            {
               store_[get_ordinal_index(i,j)] += spin_array(i+s,j+s);
            }
         }
   }

   T& operator() (int i, int j)
   {
      return store_[get_ordinal_index(i,j)];
   }

   const T& operator() (int i, int j) const
   {
      return store_[get_ordinal_index(i,j)];
   }

   void resize (size_t nMO)
   {
      extent_ = nMO;
      store_.clear();
      store_.resize(nMO*(nMO+1)/2, static_cast<T>(0));
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

} // namespace Npdm
} // namespace SpinAdapted

// Array class for 3PDM
#include <npdm_symmetric_spatial_array_2d.h>

// Array class for 3PDM
#include <npdm_symmetric_spatial_array_3d.h>

// Array class for 4PDM
//#include <npdm_symmetric_spatial_array_4d.hpp>

#endif // __NPDM_SYMMETRIC_SPATIAL_ARRAY_H
