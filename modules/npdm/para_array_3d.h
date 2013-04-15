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

template<class T> class para_array_3d;

// Parallel 3d array class
// -----------------------
// Mimics triang_2d array class

template<class T> class para_array_3d : public para_sparse_vector<T>
{
public:
  para_array_3d() : stored_local(true), upper_triangular(false) {}

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
    global_index_tuple.clear();
    local_index_tuple.clear(); 
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
      assert(has_local_index(trimap(i, j)));
    return store[trimap(i, j)];
  }
  const T& operator()(int i, int j, int k=-1) const
  {
    assert (i >= j);
    assert(has(i, j));
    if (!stored_local)
      assert(has_local_index(trimap(i, j)));
    return store[trimap(i, j)];
  }

  /// query whether elements are non-null
  bool has(int i, int j, int k=-1) const { return has_global_index(trimap(i, j)); }
  bool has(const std::vector<int>& orbs) const
  {
    assert(orbs.size() == 2);
    return has(orbs[0], orbs[1]);
  }
  bool has_global_index(int i, int j, int k=-1) const
  {
    return has_global_index(trimap(i, j));
  }
  bool has_local_index(int i, int j, int k=-1) const
  {
    return has_local_index(trimap(i, j));
  }
  bool has_global_index(int i) const
  {
    return (global_indices_map[i] != -1);
  }
  bool has_local_index(int i) const
  {
    return (local_indices_map[i] != -1);
  }
  
  friend std::ostream& operator<<(std::ostream& os, para_array_3d& op)
  {
    abort();
    return os;
  }

  /// returns 1d index from i, j
  int trimap(int i, int j) const
  {
    return ::trimap(i, j, length, upper_triangular);
  }

  /// returns i j for ith element of global storage
  std::tuple<int,int,int> unmap_global_index(int i)
  {
    return global_index_tuple[i];
  }
  /// returns i j for ith element of local storage
  //FIXME const std::tuple<int,int,int> unmap_local_index(int i) const
  const std::pair<int,int> unmap_local_index(int i) const
  {
    return local_index_tuple[i];
  }

  para_array_3d<T>* clone() const { return new para_array_3d<T>(*this); }

  void add_local_indices(int i, int j)
  {
    int index = trimap(i, j);
    local_indices.push_back(index);
    local_indices_map[index]= index;
  }

  /**
  * make a sparse para_array_triang2d with specified
  * non-zero indices
  *
  * see corresponding set_indices member fn.
  */
//MAW FIXME tuple
  void set_tuple_indices(const std::vector<std::pair<int,int> >& occupied, 
			int len,
			bool ut = false)
  {
    clear();

    length = len;
    upper_triangular = ut;

    int length_1d = tristore(len);

    /* this part is different from set_indices */
//MAW FIXME tuple
    for (std::vector<std::pair<int,int> >::const_iterator ptr = occupied.begin(); ptr != occupied.end(); ++ptr) {
      global_indices.push_back(trimap(std::get<0>(*ptr), std::get<1>(*ptr)));
      global_index_tuple.push_back(*ptr);
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
    int length_1d = tristore(length);
    if (stored_local)
      {
	local_indices = global_indices;
	local_indices_map = global_indices_map;
	local_index_tuple = global_index_tuple;
      }
    else
      {
	local_indices_map.resize(length_1d);
	int rank = mpigetrank();
	for (int i = 0; i < length_1d; ++i) local_indices_map[i] = -1;
	for (int i = 0; i < global_indices.size(); ++i)
	  if (processorindex(global_indices[i]) == rank)
	    {
	      local_indices.push_back(global_indices[i]);
	      local_indices_map[global_indices[i]] = global_indices[i];
	      local_index_tuple.push_back(global_index_tuple[i]);
	    }

      }  
  }
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & stored_local & upper_triangular & length & global_indices & global_indices_map
       & local_indices & local_indices_map & global_index_tuple & local_index_tuple & store;
  }

  std::vector<int> global_indices;
  std::vector<int> global_indices_map;
  std::vector<int> local_indices;
  std::vector<int> local_indices_map;
//MAW FIXME tuple
  std::vector<std::pair<int,int> > global_index_tuple; 
//MAW FIXME tuple
  std::vector<std::pair<int,int> > local_index_tuple;
  std::vector<T> store;
  bool stored_local;
  bool upper_triangular;
  int length;

};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

#endif
