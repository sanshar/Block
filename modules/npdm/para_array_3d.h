/*
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012
Copyright (c) 2012, Garnet K.-L. Chan

This program is integrated in Molpro with the permission of
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef PARA_ARRAY_3D_H
#define PARA_ARRAY_3D_H

#include <tuple>
#include <vector>
//#include <iostream>
//#include <utility>
//#ifndef SERIAL
//#include <boost/mpi/communicator.hpp>
//#endif
//#include <communicate.h>
//#include <multiarray.h>
//#include <boost/serialization/serialization.hpp>


//-----------------------------------------------------------------------------------------------------------------------------------------------------------

template<class T> class para_array_3d;

// Parallel 3d array class
// -----------------------
// Mimics triang_2d array class

template<class T> class para_array_3d : public para_sparse_vector<T>
{
public:
  para_array_3d() : stored_local_(true), upper_triangular_(false) {}

  bool is_upper() const { return upper_triangular_; }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  /// clears all elements
  void clear()
  {
    global_indices_.clear();
    global_indices_map_.clear();
    local_indices_.clear();
    local_indices_map_.clear();
    global_index_tuple_.clear();
    local_index_tuple_.clear();
    store_.clear();
    length_ = 0;
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  /// set storage flags
  bool& set_local() { return stored_local_; }
  /// is storage distributed?
  bool is_distributed() const { return !stored_local_; }
  /// is storage local?
  bool is_local() const { return stored_local_; }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  /// number of non-null elements in local storage
  int local_nnz() const { return local_indices_.size(); }

  /// number of non-null elements in global storage
  int global_nnz() const { return global_indices_.size(); }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

   inline int tristore(int i)
   {
     return i * (i+1) * (i+2) / 6;
   }
   
   //inline int trimap(int i, int j, int length_, bool ut = false)
   //inline int trimap(int i, int j, int k) const { return ::trimap(i,j,k, length_, upper_triangular_); }

   // Returns 1d index from i,j,k
   // e.g. length_ = 3; i=j=k=0 => 0; i=j=k=2 => 9
   inline int trimap(int i, int j, int k) const
   {
     assert ( i<=j );
     assert ( j<=k );
     int base_k = k*(k+1)*(k+2)/6;
     int base_j = j*(j+1)/2;
     return i + base_k + base_j;
   } 
//     if (i>=j)
//       {
//         //if (ut)
//         //return tristore(i) + j;
//   
//         int halflen = length_/2;
//         //there are three slots
//         //slot1 i >= halflen and j>= halflen
//         //slot2 i<halflen and j < halflen
//         //slot3 i >= halflen and j<halflen
//   
//         //first check if our case is in slot 1
//         if (i>=halflen && j >= halflen)
//           return tristore(length_ - j - 1) + length_ - i - 1;
//         else if (i < halflen && j <halflen)
//           return tristore(length_ - halflen - 1) + length_ - halflen + tristore(i) + j;
//         else {
//           int base= tristore(length_ - halflen - 1) + length_ - halflen + tristore(halflen);
//           return base + (i-halflen)*(halflen) + (j);
//         }
//         //return tristore(length_ - j - 1) + length_ - i - 1;
//       }
//     else
//       assert(false);
//     return 0;
//   }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//FIXME change all .at(ix) to [ix]

  /// query whether elements are non-null
  bool has(int i, int j, int k) const { return has_global_index( trimap(i,j,k) ); }
  bool has(const std::vector<int>& orbs) const { assert(orbs.size() == 3); return has(orbs[0], orbs[1], orbs[2]); }

  bool has_global_index(int i, int j, int k) const { return has_global_index( trimap(i,j,k) ); }
  bool has_local_index(int i, int j, int k) const { return has_local_index( trimap(i,j,k) ); }
  bool has_global_index(int i) const { return (global_indices_map_.at(i) != -1); }
  bool has_local_index(int i) const { return (local_indices_map_.at(i) != -1); }

  const std::vector<int>& get_indices() const { return global_indices_; }
  std::vector<int>& get_local_indices_() { return local_indices_; }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
  /// returns elements at orbs
  T& operator()(const std::vector<int>& orbs)
  {
    int i = orbs[0];
    int j = orbs[1];
    int k = orbs[2];
    return (*this)(i,j,k);
  }

  const T& operator()(const std::vector<int>& orbs) const
  {
    int i = orbs[0];
    int j = orbs[1];
    int k = orbs[2];
    return (*this)(i,j,k);
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
  /// returns elements at i,j,k
  T& operator()(int i, int j, int k)
  {
    assert (i <= j);
    assert (j <= k);
    assert( has(i,j,k) );
    if (!stored_local_)
      assert( has_local_index( trimap(i,j,k) ) );
    return store_.at( trimap(i,j,k) );
  }

  const T& operator()(int i, int j, int k) const
  {
    assert (i <= j);
    assert (j <= k);
    assert(has(i,j,k));
    if (!stored_local_)
      assert( has_local_index( trimap(i,j,k) ) );
    return store_.at(trimap(i,j,k));
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  /// ith element of local storage
  T& get_local_element(int i) { return store_.at( local_indices_.at(i) ); }
  const T& get_local_element(int i) const { return store_.at( local_indices_.at(i) ); }

  /// ith element of global storage
  T& get_global_element(int i) { return store_.at( global_indices_.at(i) ); }
  const T& get_global_element(int i) const { return store_.at( global_indices_.at(i) ); }

  T& get(const std::vector<int>& orbs) { return (*this)(orbs[0], orbs[1], orbs[2]); }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  /// returns i,j,k for ith element of global storage
  std::vector<int> unmap_global_index(int i) { return global_index_tuple_.at(i); }

  /// returns i,j,k for ith element of local storage
  const std::vector<int> unmap_local_index(int i) const { return local_index_tuple_.at(i); }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  friend std::ostream& operator<<(std::ostream& os, para_array_3d& op) { assert(false); return os; }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  para_array_3d<T>* clone() const { return new para_array_3d<T>(*this); }

  void add_local_indices(int i, int j, int k)
  {
    int index = trimap(i,j,k);
    local_indices_.push_back(index);
    local_indices_map_.at(index)= index;
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
  /**
  * make a sparse para_array_3d with specified non-zero indices
  *
  * see corresponding set_indices member fn.
  */
  void set_tuple_indices(const std::vector<std::tuple<int,int,int> >& occupied, int len, bool ut = false)
  {
    clear();
    length_ = len;
    upper_triangular_ = ut;

    // Global indices
    for (auto it = occupied.begin(); it != occupied.end(); ++it) {
      // internally implement tuple as std::vector<int>(3)
      std::vector<int> tuple = { std::get<0>(*it), std::get<1>(*it), std::get<2>(*it) };
      global_indices_.push_back( trimap( tuple[0], tuple[1], tuple[2] ) );
      global_index_tuple_.push_back( tuple );
    }

    // Global indices map
    int length_1d = tristore(len);
    global_indices_map_.resize(length_1d);
    for (int i=0; i < length_1d; ++i) {
      global_indices_map_.at(i) = -1;
    }
    for (int i=0; i < global_indices_.size(); ++i) {
      global_indices_map_.at( global_indices_.at(i) ) = global_indices_.at(i);
    }
    store_.resize(length_1d);

    // Now setup local indices
    setup_local_indices();
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

private:
  /**
   * having filled out the global indices, assign indices
   * to individual processors (i.e. local indices)
   */
  void setup_local_indices()
  {
    int length_1d = tristore(length_);
    if (stored_local_) {
      local_indices_ = global_indices_;
      local_indices_map_ = global_indices_map_;
      local_index_tuple_ = global_index_tuple_;
    }
    else {
      local_indices_map_.resize(length_1d);
      int rank = mpigetrank();
      for (int i = 0; i < length_1d; ++i) local_indices_map_.at(i) = -1;
      for (int i = 0; i < global_indices_.size(); ++i) {
        if (processorindex(global_indices_.at(i)) == rank) {
          local_indices_.push_back(global_indices_.at(i));
          local_indices_map_.at(global_indices_.at(i)) = global_indices_.at(i);
          local_index_tuple_.push_back(global_index_tuple_.at(i));
        }
      }
    }
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & stored_local_
       & upper_triangular_
       & length_
       & global_indices_
       & global_indices_map_
       & local_indices_
       & local_indices_map_
       & global_index_tuple_
       & local_index_tuple_
       & store_;
  }

  std::vector<int> global_indices_;
  std::vector<int> global_indices_map_;
  std::vector<int> local_indices_;
  std::vector<int> local_indices_map_;
  std::vector<std::vector<int> > global_index_tuple_;
  std::vector<std::vector<int> > local_index_tuple_;
  std::vector<T> store_;
  bool stored_local_;
  bool upper_triangular_;
  int length_;

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

};

#endif
