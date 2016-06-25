/*
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012
Copyright (c) 2012, Garnet K.-L. Chan

This program is integrated in Molpro with the permission of
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef FOURPDM_PARA_ARRAY_H
#define FOURPDM_PARA_ARRAY_H

#include <tuple>
#include <vector>
#include "pario.h"

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

template<class T> class para_array_4d;

// Parallel 4d array class
// -----------------------
// Mimics triang_2d array class

template<class T> class para_array_4d : public para_sparse_vector<T>
{
public:
  para_array_4d() : stored_local_(false) {}

  // This is designed for 4-index operators
  int num_indices() { return 4; }

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

   inline int tristore_4d(int p)
   {
     return p * (p+1) * (p+2) * (p+3) / 24;
   }
   
   // Returns 1d index from i,j,k,l
   // l-index has stride one
   inline int trimap_4d(int i, int j, int k, int l) const
   {
     assert ( l<=k );
     assert ( k<=j );
     assert ( j<=i );
     int base_i = i*(i+1)*(i+2)*(i+3)/24;
     int base_j = j*(j+1)*(j+2)/6;
     int base_k = k*(k+1)/2;
     return l + base_k + base_j + base_i;
   } 

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//FIXME change all .at(ix) to [ix]

  /// query whether elements are non-null
  bool has(int i, int j, int k, int l) const { return has_global_index( trimap_4d(i,j,k,l) ); }
  bool has(const std::vector<int>& orbs) const { assert(orbs.size() == 4); return has(orbs[0], orbs[1], orbs[2], orbs[3]); }

  bool has_global_index(int i, int j, int k, int l) const { return has_global_index( trimap_4d(i,j,k,l) ); }
  bool has_local_index(int i, int j, int k, int l) const { return has_local_index( trimap_4d(i,j,k,l) ); }
  bool has_global_index(int i) const { return (global_indices_map_.at(i) != -1); }
  bool has_local_index(int i) const { return (local_indices_map_.at(i) != -1); }

  const std::vector<int>& get_indices() const { return global_indices_; }
  const std::vector<int>& get_local_indices() const { return local_indices_; }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
  /// returns elements at orbs
  T& operator()(const std::vector<int>& orbs)
  {
    int i = orbs[0];
    int j = orbs[1];
    int k = orbs[2];
    int l = orbs[3];
    return (*this)(i,j,k,l);
  }

  const T& operator()(const std::vector<int>& orbs) const
  {
    int i = orbs[0];
    int j = orbs[1];
    int k = orbs[2];
    int l = orbs[3];
    return (*this)(i,j,k,l);
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
  /// returns elements at i,j,k,l
  T& operator()(int i, int j, int k, int l)
  {
    assert (l <= k);
    assert (k <= j);
    assert (j <= i);
    assert( has(i,j,k,l) );
    if (!stored_local_)
      assert( has_local_index( trimap_4d(i,j,k,l) ) );
    return store_.at( trimap_4d(i,j,k,l) );
  }

  const T& operator()(int i, int j, int k, int l) const
  {
    assert (l <= k);
    assert (k <= j);
    assert (j <= i);
    assert(has(i,j,k,l));
    if (!stored_local_)
      assert( has_local_index( trimap_4d(i,j,k,l) ) );
    return store_.at(trimap_4d(i,j,k,l));
  }

  const std::vector<T>& get_store() const  { return store_; }
  
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  /// ith element of local storage
  T& get_local_element(int i) { return store_.at( local_indices_.at(i) ); }
  const T& get_local_element(int i) const { return store_.at( local_indices_.at(i) ); }

  /// ith element of global storage
  T& get_global_element(int i) { return store_.at( global_indices_.at(i) ); }
  const T& get_global_element(int i) const { return store_.at( global_indices_.at(i) ); }

  T& get(const std::vector<int>& orbs) { return (*this)(orbs[0], orbs[1], orbs[2], orbs[3]); }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  /// returns i,j,k,l for ith element of global storage
  std::vector<int> unmap_global_index(int i) { return global_index_tuple_.at(i); }

  /// returns i,j,k,l for ith element of local storage
  const std::vector<int> unmap_local_index(int i) const { return local_index_tuple_.at(i); }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  friend std::ostream& operator<<(std::ostream& os, para_array_4d& op) { abort(); return os; }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  para_array_4d<T>* clone() const { return new para_array_4d<T>(*this); }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  void add_local_indices(int i, int j, int k, int l)
  {
    int index = trimap_4d(i,j,k,l);
    std::vector<int> tuple = { i,j,k,l };
    local_indices_.push_back(index);
    local_indices_map_.at(index) = index;
    local_index_tuple_.push_back( tuple );
    //FIXME MAW also update global, even though this wasn't done before
    global_indices_.push_back(index);
    global_indices_map_.at(index) = index;
    global_index_tuple_.push_back( tuple );
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
  /**
  * make a sparse para_array_4d with specified non-zero indices
  *
  * see corresponding set_indices member fn.
  */
  void set_tuple_indices(const std::map< std::tuple<int,int,int,int>, int >& occupied, int len)
  {
    clear();
    length_ = len;

    for (auto it = occupied.begin(); it != occupied.end(); ++it) {
      std::vector<int> tuple = { std::get<0>(it->first), std::get<1>(it->first), std::get<2>(it->first), std::get<3>(it->first) };
      int rank = it->second; 
      // Global indices
      int idx = trimap_4d( tuple[0], tuple[1], tuple[2], tuple[3] );
      global_indices_.push_back( idx );
      global_index_tuple_.push_back( tuple );
      // Local indices
      if ( stored_local_ || rank == mpigetrank() ) {
        // Assign to requested rank
        local_indices_.push_back( idx );
        local_index_tuple_.push_back( tuple );
      }
      else if ( rank == -1 ) {
        // Load balance equally over ranks using global rule
        if ( processorindex( idx ) == mpigetrank() ) {
          local_indices_.push_back( idx );
          local_index_tuple_.push_back( tuple );
        }
      }
    }
//pout << "local index tuples:\n";
//for (auto it = local_index_tuple_.begin(); it != local_index_tuple_.end(); ++it) {
//  pout << "p" << mpigetrank() << ": " << it->at(0) << "," << it->at(1) << "," << it->at(2) << endl;
//}

    // Global indices map
    int length_1d = tristore_4d(len);
    global_indices_map_.resize(length_1d);
    for (int i=0; i < length_1d; ++i) {
      global_indices_map_.at(i) = -1;
    }
    for (int i=0; i < global_indices_.size(); ++i) {
      global_indices_map_.at( global_indices_.at(i) ) = global_indices_.at(i);
    }
    store_.resize(length_1d);

    // Local indices map
    local_indices_map_.resize(length_1d);
    for (int i=0; i < length_1d; ++i) {
      local_indices_map_.at(i) = -1;
    }
    for (int i=0; i < local_indices_.size(); ++i) {
      local_indices_map_.at( local_indices_.at(i) ) = local_indices_.at(i);
    }

  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

private:

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//  /**
//   * having filled out the global indices, assign indices
//   * to individual processors (i.e. local indices)
//   */
//  void setup_local_indices( const std::map< std::tuple<int,int,int>, int >& occupied)
//  {
//    int length_1d = tristore_4d(length_);
//    if (stored_local_) {
//      abort();
////      local_indices_ = global_indices_;
////      local_indices_map_ = global_indices_map_;
////      local_index_tuple_ = global_index_tuple_;
////    }
//    else {
//      local_indices_map_.resize(length_1d);
//      for (int p = 0; p < length_1d; ++p) local_indices_map_.at(p) = -1;
//
//      for (int p = 0; p < global_indices_.size(); ++p) {
//       
//
//        if ( processorindex(global_indices_[p]) == mpigetrank() ) {
//          local_indices_.push_back(global_indices_.at(p));
//          local_indices_map_.at(global_indices_.at(p)) = global_indices_.at(p);
//          local_index_tuple_.push_back(global_index_tuple_.at(p));
//        }
//
//      }
//    }
//  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & stored_local_
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
  int length_;

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

};

#endif
