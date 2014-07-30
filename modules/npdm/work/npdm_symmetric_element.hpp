#ifndef __NPDM_SYMMETRIC_ELEMENT_HPP
#define __NPDM_SYMMETRIC_ELEMENT_HPP

namespace SpinAdapted {
namespace Npdm {

//  Forward decralation
template<typename T> struct const_symmetric_element;

/// Element reference with parity operation
template<typename T>
struct symmetric_element
{
   explicit
   symmetric_element (T& value, bool parity = false) : parity_(parity), ref_(&value) { }

   symmetric_element (const symmetric_element& x) : parity_(x.parity_), ref_(x.ref_) { }

   void operator= (const T& value) { *ref_ = parity_ ? -value : value; }

   symmetric_element& operator= (const symmetric_element& x)
   {
      *ref_ = (parity_ ^ x.parity_) ? -(*x.ref_) : (*x.ref_);
      return *this;
   }

   symmetric_element& operator= (const const_symmetric_element<T>& x)
   {
      *ref_ = (parity_ ^ x.parity_) ? -(*x.ref_) : (*x.ref_);
      return *this;
   }

   operator T () const { return parity_ ? -(*ref_) : (*ref_); }

private:

   friend class const_symmetric_element<T>;

   bool parity_;

   T* ref_;
};

/// Element const reference with parity operation
template<typename T>
struct const_symmetric_element
{
   explicit
   const_symmetric_element (const T& value, bool parity = false) : parity_(parity), ref_(&value) { }

   const_symmetric_element (const const_symmetric_element& x) : parity_(x.parity_), ref_(x.ref_) { }

   const_symmetric_element (const symmetric_element<T>& x) : parity_(x.parity_), ref_(x.ref_) { }

   operator T () const { return parity_ ? -(*ref_) : (*ref_); }

private:

   friend class symmetric_element<T>;

   bool parity_;

   const T* ref_;
};

} // namespace Npdm
} // namespace SpinAdapted

#endif // __NPDM_SYMMETRIC_ELEMENT_HPP
