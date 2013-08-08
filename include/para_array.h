/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef PARA_ARRAY_H 
#define PARA_ARRAY_H 
#include <vector>
#include <iostream>
#include <utility>
#ifndef SERIAL
#include <boost/mpi/communicator.hpp>
#endif
#include <communicate.h>
#include <multiarray.h>     
#include <boost/serialization/serialization.hpp>

/**
 * forward declarations
 * 
 */
template<class T> struct para_sparse_vector;
template<class T> class para_array_0d;
template<class T> class para_array_1d;
template<class T> class para_array_triang_2d;

// utility functions for communication
inline int processorindex(int i)
{
#ifdef SERIAL
  return 0;
#else
  boost::mpi::communicator world;
  int size = world.size();
  return i % size;
#endif
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------

/// 1d interface for the parallel array classes
template<class T> struct para_sparse_vector
{
public:
  /// clear all elements
  virtual void clear() {};
  /// number of non-zero elements in local storage
  virtual int local_nnz() const {};

  virtual int global_nnz() const {};
  /// ith element of local storage
  virtual const T& get_local_element(int i) const {};
  /// ith element of local storage
  virtual T& get_local_element(int i) {};
  /// query nullness of element i
  virtual bool has_local_index(int i) const {};

  // virtual indexing
  virtual T& operator()(const std::vector<int>& orbs) {};
  virtual const T& operator()(const std::vector<int>& orbs) const {};

  /// expose local storage
  virtual bool is_local() const {};
  virtual bool is_distributed() const {};

  virtual bool& set_local() {};

  /// virtual constructor. caller is responsible for managing storage
  virtual para_sparse_vector<T>* clone() const {};

  //FIX ME!!
  virtual bool has_local_index(int i, int j) const {};

  virtual bool has(int i) const {};
  virtual bool has(int i, int j) const {};
  virtual bool has(int i=-1, int j=-1, int k=-1) const {};
  virtual bool has(const std::vector<int>& orbs) const {};
  virtual const std::vector<T>& get_store() const {};
  virtual const std::vector<int>& get_indices() const {};
  virtual int trimap_2d(int i, int j) const {};
  virtual T& get(const std::vector<int>& orbs)=0;
//MAW
//  virtual const std::pair<int, int> unmap_local_index(int i) const { assert(false); };
  virtual const std::vector<int> unmap_local_index(int i) const { assert(false); };
};

//------------------------------------------------------------------------------------------------------------------------------------------------------------

/// A wrapper for a single element 
/**
 * This class enables single elements to be used as if they
 * were arrays with 1 element, allowing them to have
 * the same interface as para_sparse_vector, to be used in
 * for_all loops, etc..
 * Functions are implemented in terms of the para_array_1d functions.
 */
template<class T> class para_array_0d : public para_sparse_vector<T>
{
public:

  /// implements para_sparse_vector interface by forwarding to para_array_1d
  const int num_indices() { return 0; }

  void clear() { store.clear(); }
  int local_nnz() const { return store.local_nnz(); }
  int global_nnz() const { return store.global_nnz(); }
  const T& get_local_element(int i) const { return store.get_local_element(i); }
  T& get_local_element(int i) { return store.get_local_element(i); }
  const T& get_global_element(int i) const { return store.get_global_element(i); }
  T& get_global_element(int i) { return store.get_global_element(i); }
  T& operator()(const std::vector<int>& orbs) { return store.get_local_element(orbs[0]); }
  const T& operator()(const std::vector<int>& orbs) const { return store.get_local_element(orbs[0]); }
  T& operator()(int i=-1, int j=-1, int k=-1) { return store(0); }
  const T& operator()(int i=-1, int j=-1, int k=-1) const { return store(0); }
  bool has_local_index(int i, int j=-1, int k=-1) const { return store.has_local_index(i); }
  const std::vector<T>& get_store() const { return store.get_store(); }
  void set_indices() 
  { 
    std::vector<int> i(1); i[0] = 0;
    store.set_indices(i, 1);
  }
  T& get(const std::vector<int>& orbs)
  {
    return store(0);
  }
  const T& operator()() const { return store(0); }
  T& operator()() { return store(0); }
  bool is_local() const { return store.is_local(); }
  bool is_distributed() const { return store.is_distributed(); }
  para_array_0d<T>* clone() const { return new para_array_0d<T>(*this); }
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & store;
  }
  para_array_1d<T> store;
};

//------------------------------------------------------------------------------------------------------------------------------------------------------------

/// A sparse parallel 1d array class
/**
 * 
 * The idea behind this array is that the elements have two states, one which is
 * large, one which is null. The storage therefore contains a FULL array of
 * the null objects; set_indices then sets up the array by telling it which indices
 * will contain LARGE objects. The non-null element indices are then stored in 
 * global_indices/local_indices.
 *
 * The parallelism works via local_indices being different from global_indices. all
 * processors storing the same parallel array have the same global_indices, but may
 * have different local_indices. The local_indices are partitioned equally amongst the 
 * processors in set_indices.
 *
 * The index maps store -1 if the index is invalid, otherwise they store the value of 
 * the index. e.g. if an array holds elements 2, 3, then index_map[0] = -1, but 
 * index_map[3] = 3.
 */
/// TODO: there is a lot of common ground between the different array classes: this should
/// be exploited.
template<class T> class para_array_1d : public para_sparse_vector<T>
{
public:
  para_array_1d() : stored_local(true) {}
  
  // This is designed for 1-index operators
  const int num_indices() { return 1; }

  /// clears all elements
  void clear()
  {
    global_indices.clear();
    global_indices_map.clear();
    local_indices.clear();
    local_indices_map.clear();
    store.clear();
  }
  const std::vector<int>& get_local_indices() const { return local_indices; }

  /// query whether elements are non-null, locally and globally
  bool has(int i, int j=-1, int k=-1) const
  {
    return has_global_index(i); 
  }
  bool has_local_index(int i, int j=-1, int k=-1) const 
  { 
    return local_indices_map[i] != -1; 
  }
  bool has_global_index(int i) const
  {
    return (global_indices_map[i] != -1);
  }
  bool has(const std::vector<int>& orbs) const
  {
    assert(orbs.size() == 1);
    return (global_indices_map[orbs[0]] != -1);      
  }
  
  void set_indices(const std::vector<int>& inds, int length)
  {    
    global_indices.clear();
    global_indices_map.clear();
    local_indices.clear();
    local_indices_map.clear();
      
    global_indices = inds;
    global_indices_map.resize(length);
    for (int i = 0; i < length; ++i) global_indices_map[i] = -1;
    for (int i = 0; i < inds.size(); ++i)
      global_indices_map[global_indices[i]] = global_indices[i];


//MAW
//std::cout << "1D para_array; global_indices_map: \n";
//for (int p=0; p<global_indices_map.size(); ++p) {
//  std::cout << p << "\t\t" << global_indices_map.at(p) << std::endl;
//}

    if (stored_local) {
//cout << "1D para_array stored local\n";
      local_indices = global_indices;
    }
    else {
//cout << "1D para_array spread over procs\n";
      for (int i = 0; i < global_indices.size(); ++i) {
        if (processorindex(global_indices[i]) == mpigetrank()) {
//cout << "global_indices[i], mpigetrank() = " << global_indices[i] << "; " <<  mpigetrank() << endl;
          local_indices.push_back(global_indices[i]);
        }
      }
    }

    local_indices_map.resize(length);
    for (int i = 0; i < length; ++i) local_indices_map[i] = -1;
    for (int i = 0; i < local_indices.size(); ++i) 
      local_indices_map[local_indices[i]] = local_indices[i];

    store.resize(length);    
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  void add_local_index(int i)
  {
    local_indices.push_back(i);
    local_indices_map[i] = i;
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  /// storage labels
  bool& set_local() { return stored_local; }
  bool is_local() const { return stored_local; }
  bool is_distributed() const { return !stored_local; }

  const std::vector<T>& get_store() const { return store; }  /// deprecated
  const std::vector<int>& get_indices() const { return global_indices; } /**< deprecated */

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  /// returns value at index i
  const T& operator()(int i, int j=-1, int k=-1) const
  {
    assert(has_global_index(i));
    if (is_distributed())		/**< if storage is distributed, make sure index is available on current processor */
      assert(has_local_index(i));
    return store[i];
  }    

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  /// returns value at index i
  T& operator()(int i, int j=-1, int k=-1)
  {
    assert(has_global_index(i));
    if (is_distributed())
      assert(has_local_index(i));
    return store[i];
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  /// ith element of local storage
  T& get_local_element(int i)
  {
    return store[local_indices[i]];
  }
  
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  /// ith element of local storage
  const T& get_local_element(int i) const
  {
    return store[local_indices[i]];
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  /// returns the ith element of global storage
  T& get_global_element(int i)
  {
    return store[global_indices[i]];
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  /// returns the ith element of global storage
  const T& get_global_element(int i) const
  {
    return store[global_indices[i]];
  }

  /// total non-sparse size >= global_nnz
  int max_size() { return global_indices_map.size(); }

  /// number of non-null elements in local storage
  int local_nnz() const { return local_indices.size(); }

  /// number of non-null elements in global storage
  int global_nnz() const { return global_indices.size(); }

  T& get(const std::vector<int>& orbs)
  {
    return store[orbs[0]];
  }
  
  T& operator()(const std::vector<int>& orbs)
  {
    int i = orbs[0]; 
    return (*this)(i);
  }
  const T& operator()(const std::vector<int>& orbs) const
  {
    int i = orbs[0]; 
    return (*this)(i);
  }
  

  friend std::ostream& operator<<(std::ostream& os, para_array_1d& op) 
  {
    for (int i = 0; i < op.max_size(); ++i) os << op.store[i] << std::endl;
    return os;
  }      

  para_array_1d<T>* clone() const { return new para_array_1d<T>(*this); }

private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & stored_local & global_indices & global_indices_map & local_indices & local_indices_map & store;
  }
  std::vector<int> global_indices;	/**< non-zero global indices */
  std::vector<int> global_indices_map;	
  std::vector<int> local_indices; /**< non-zero local indices */
  std::vector<int> local_indices_map; 
  std::vector<T> store;		/**< element storage (length >= global_nnz) */
  bool stored_local;		/**< are all elements stored locally, or are they distributed */
};


//===========================================================================================================================================================


inline int tristore_2d(int i)
{
  return i * (i + 1) / 2;
}

inline int trimap_2d(int i, int j, int length, bool ut = false)
{
  if (i>=j) 
    {
      //if (ut)
      //return tristore_2d(i) + j;
      
      int halflen = length/2;
      //there are three slots
      //slot1 i >= halflen and j>= halflen
      //slot2 i<halflen and j < halflen
      //slot3 i >= halflen and j<halflen
      
      //first check if our case is in slot 1
      if (i>=halflen && j >= halflen)
        return tristore_2d(length - j - 1) + length - i - 1;
      else if (i < halflen && j <halflen)
        return tristore_2d(length - halflen - 1) + length - halflen + tristore_2d(i) + j;
      else {
        int base= tristore_2d(length - halflen - 1) + length - halflen + tristore_2d(halflen);
        return base + (i-halflen)*(halflen) + (j);
      }
      //return tristore_2d(length - j - 1) + length - i - 1;
    }
  else 
    assert(false);
  return 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------


/// parallel 2d lower triangular array class
/**
 * This has a flat layout; the corresponding 1d storage can the be distributed
 * over different processors. See para_array_1d and para_array_2d for more details.
 */
template<class T> class para_array_triang_2d : public para_sparse_vector<T>
{
public:
  para_array_triang_2d() : stored_local(true), upper_triangular(false) {}

  // This is designed for 2-index operators
  const int num_indices() { return 2; }

  /// exposes storage
  const std::vector<T>& get_store() const { return store; }	/**< deprecated */

  bool is_upper() const { return upper_triangular; }

  /// clears all elements
  void clear()
  {
    global_indices.clear();
    global_indices_map.clear();
    local_indices.clear();
    local_indices_map.clear();
    global_index_pair.clear();
    local_index_pair.clear(); 
    store.clear();
    length = 0;
  }

  /// set storage flags
  bool& set_local() { return stored_local; }
  /// is storage distributed?
  bool is_distributed() const { return !stored_local; }
  /// is storage local?
  bool is_local() const { return stored_local; }

  const std::vector<int>& get_indices() const { return global_indices; }

  std::vector<int>& get_local_indices() { return local_indices; }


  /// number of non-null elements in local storage
  int local_nnz() const { return local_indices.size(); }

  /// number of non-null elements in global storage
  int global_nnz() const { return global_indices.size(); }

  /// ith element of local storage
  T& get_local_element(int i)
  {
    return store[local_indices[i]];
  }
  const T& get_local_element(int i) const
  {
    return store[local_indices[i]];
  }

  T& get_global_element(int i)
  {
    return store[global_indices[i]];
  }

  /// ith element of global storage
  const T& get_global_element(int i) const
  {
    return store[global_indices[i]];
  }


  T& get(const std::vector<int>& orbs)
  {
    return (*this)(orbs[0], orbs[1]);
  }

  /// returns elements at orbs
  T& operator()(const std::vector<int>& orbs)
  {
    int i = orbs[0]; int j = orbs[1];
    return (*this)(i, j);
  }
  const T& operator()(const std::vector<int>& orbs) const
  {
    int i = orbs[0]; int j = orbs[1];
    return (*this)(i, j);
  }
  /// returns elements at i, j
  T& operator()(int i, int j, int k=-1)
  {
    assert (i >= j);
    assert(has(i, j));
    if (!stored_local)
      assert(has_local_index(trimap_2d(i, j)));
    return store[trimap_2d(i, j)];
  }
  const T& operator()(int i, int j, int k=-1) const
  {
    assert (i >= j);
    assert(has(i, j));
    if (!stored_local)
      assert(has_local_index(trimap_2d(i, j)));
    return store[trimap_2d(i, j)];
  }

  /// query whether elements are non-null
  bool has(int i, int j, int k=-1) const { return has_global_index(trimap_2d(i, j)); }
  bool has(const std::vector<int>& orbs) const
  {
    assert(orbs.size() == 2);
    return has(orbs[0], orbs[1]);
  }
  bool has_global_index(int i, int j, int k=-1) const
  {
    return has_global_index(trimap_2d(i, j));
  }
  bool has_local_index(int i, int j, int k=-1) const
  {
    return has_local_index(trimap_2d(i, j));
  }
  bool has_global_index(int i) const
  {
    return (global_indices_map[i] != -1);
  }
  bool has_local_index(int i) const
  {
    return (local_indices_map[i] != -1);
  }
  
  friend std::ostream& operator<<(std::ostream& os, para_array_triang_2d& op)
  {
    assert(false);
    return os;
  }

  /// returns 1d index from i, j
  int trimap_2d(int i, int j) const
  {
    return ::trimap_2d(i, j, length, upper_triangular);
  }

  /// returns i j for ith element of global storage
  std::pair<int, int> unmap_global_index(int i)
  {
    return global_index_pair[i];
  }
  /// returns i j for ith element of local storage
//MAW  const std::pair<int, int> unmap_local_index(int i) const
  const std::vector<int> unmap_local_index(int i) const
  {
//MAW    return local_index_pair[i];
    std::vector<int> ret(2);
    ret[0] = local_index_pair[i].first;
    ret[1] = local_index_pair[i].second;
    return ret;
  }

  para_array_triang_2d<T>* clone() const { return new para_array_triang_2d<T>(*this); }

  void add_local_indices(int i, int j)
  {
    int index = trimap_2d(i, j);
    local_indices.push_back(index);
    local_indices_map[index]= index;
    // I am not updating local_index_pair because it seems to do nothing
    //local_index_pair.push_back(global_index_pair[index]);
  }

  /**
  * make a sparse para_array_triang2d with specified
  * non-zero indices
  *
  * see corresponding set_indices member fn.
  */
  void set_pair_indices(const std::vector<std::pair<int, int> >& occupied, 
			int len,
			bool ut = false)
  {
    clear();

    length = len;
    upper_triangular = ut;

    int length_1d = tristore_2d(len);

    /* this part is different from set_indices */
    for (std::vector<std::pair<int, int> >::const_iterator ptr = occupied.begin(); ptr != occupied.end(); ++ptr) {
      global_indices.push_back(trimap_2d(ptr->first, ptr->second));
      global_index_pair.push_back(*ptr);
    }
				 
    global_indices_map.resize(length_1d);
    for (int i = 0; i < length_1d; ++i) global_indices_map[i] = -1;
    for (int i = 0; i < global_indices.size(); ++i)
      global_indices_map[global_indices[i]] = global_indices[i];
    
    store.resize(length_1d);

    // now setup local indices
    setup_local_indices();
  }

private:
  /**
   * having filled out the global indices, assign indices
   * to individual processors (i.e. local indices)
   */
  void setup_local_indices()
  {
    int length_1d = tristore_2d(length);
    if (stored_local) {
//      cout << "2D para_array stored local\n";
      local_indices = global_indices;
      local_indices_map = global_indices_map;
      local_index_pair = global_index_pair;
    }
    else {
      local_indices_map.resize(length_1d);
      int rank = mpigetrank();
//      cout << "2D para_array spread over procs\n";

      for (int i = 0; i < length_1d; ++i) local_indices_map[i] = -1;

      for (int i = 0; i < global_indices.size(); ++i) {
        if (processorindex(global_indices[i]) == rank) {
          local_indices.push_back(global_indices[i]);
          local_indices_map[global_indices[i]] = global_indices[i];
          local_index_pair.push_back(global_index_pair[i]);
//          cout << "global_index_pair[i], mpigetrank() = " << global_index_pair[i].first << "," << global_index_pair[i].second << "; " <<  mpigetrank() << endl;
        }
      }
      
    }  
  }

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & stored_local & upper_triangular & length & global_indices & global_indices_map
       & local_indices & local_indices_map & global_index_pair & local_index_pair & store;
  }

  std::vector<int> global_indices;
  std::vector<int> global_indices_map;
  std::vector<int> local_indices;
  std::vector<int> local_indices_map;
  std::vector<std::pair<int, int> > global_index_pair; 
  std::vector<std::pair<int, int> > local_index_pair;
  std::vector<T> store;
  bool stored_local;
  bool upper_triangular;
  int length;

};

//------------------------------------------------------------------------------------------------------------------------------------------------------------

#endif
