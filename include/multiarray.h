/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is free software: you can redistribute it and/or modify         
it under the terms of the GNU General Public License as published by         
the Free Software Foundation, either version 3 of the License, or            
(at your option) any later version.                                          
                                                                             
This program is distributed in the hope that it will be useful,              
but WITHOUT ANY WARRANTY; without even the implied warranty of               
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
GNU General Public License for more details.                                 
                                                                             
You should have received a copy of the GNU General Public License            
along with this program.  If not, see <http://www.gnu.org/licenses/>.        

This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef multi_array_h
#define multi_array_h
#include <vector>
#include <string>
#include <cassert>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/serialization.hpp>

using namespace std;
template<class T> class multiarray
{
public:
  virtual T& operator() (const vector<int>& indices) = 0;
  virtual T operator() (const vector<int>& indices) const = 0;
  virtual string reflect_type() const = 0;
};

template<class T> class array_2d : public vector<T>, public multiarray<T>
{
 public:
  array_2d () : dim1_d (0), dim2_d (0), vector<T> () { }
  array_2d (const int d1, const int d2) : dim1_d (d1), dim2_d (d2), vector<T> (d1 * d2) { }
  T& operator() (const int i, const int j) 
    { 
      assert ((0 <= i) && (i < dim1_d)); assert ((0 <= j) && (j < dim2_d));
      return vector<T>::operator[] (i * dim2_d + j); 
    }
  T operator() (const int i, const int j) const
    { 
      assert ((0 <= i) && (i < dim1_d)); assert ((0 <= j) && (j < dim2_d));
      return vector<T>::operator[] (i * dim2_d + j); 
    }

  T& operator() (const vector<int>& indices)
  {
    assert(indices.size() == 2);
    return operator()(indices[0], indices[1]);
  }
  T operator() (const vector<int>& indices) const
  {
    assert(indices.size() == 2);
    return operator()(indices[0], indices[1]);
  }
  void resize (const int i, const int j) { vector<T>::resize (i * j); dim1_d = i; dim2_d = j; }
  int dim1() const { return dim1_d; }
  int dim2() const { return dim2_d; }
  string reflect_type() const { return string("array_2d"); }

 private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & dim1_d & dim2_d;
    ar & boost::serialization::base_object<std::vector<T> >(*this);
  }
  int dim1_d;
  int dim2_d;
};

template<class T> class array_3d : public vector<T>, public multiarray<T>
{
 public:
  array_3d() : dim1_d (0), dim2_d (0), dim3_d (0), dim2_times_dim3_d (0), vector<T> () { }
  array_3d(const int d1, const int d2, const int d3) : dim1_d (d1), dim2_d (d2), dim3_d (d3), dim2_times_dim3_d (d2 * d3), vector<T> (d1 * d2 * d3) { }
  T& operator() (const int i, const int j, const int k) 
    { 
      assert((0 <= i) && (i < dim1_d)); assert((0 <= j) && (j < dim2_d)); assert((0 <= k) && (k < dim3_d));
      return vector<T>::operator[](i * dim2_times_dim3_d + j * dim3_d + k); 
    }
  T operator() (const int i, const int j, const int k) const
    { 
      assert((0 <= i) && (i < dim1_d)); assert((0 <= j) && (j < dim2_d)); assert((0 <= k) && (k < dim3_d));
      return vector<T>::operator[](i * dim2_times_dim3_d + j * dim3_d + k); 
    }
  T& operator() (const vector<int>& indices)
  {
    assert(indices.size() == 3);
    return operator()(indices[0], indices[1], indices[2]);
  }
  T operator() (const vector<int>& indices) const
  {
    assert(indices.size() == 3);
    return operator()(indices[0], indices[1], indices[2]);
  }

  void resize(const int i, const int j, const int k) { vector<T>::resize(i * j * k); dim1_d = i; dim2_d = j; dim3_d = k; dim2_times_dim3_d = dim2_d * dim3_d; }
  int dim1() const { return dim1_d; }
  int dim2() const { return dim2_d; }
  int dim3() const { return dim3_d; }
  string reflect_type() const { return string("array_3d"); }

 private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & dim1_d & dim2_d & dim3_d & dim2_times_dim3_d;
    ar & boost::serialization::base_object<vector<T> >(*this);
  }
  int dim1_d;
  int dim2_d;
  int dim3_d;
  int dim2_times_dim3_d;
};

template<class T> class array_4d : public vector<T>, public multiarray<T>
{
 public:
  array_4d() : dim1_d(0), dim2_d(0), dim3_d(0), dim4_d(0), dim3_times_dim4_d(0), dim2_times_dim3_times_dim4_d (0), vector<T>() { }
  array_4d(const int d1, const int d2, const int d3, const int d4) : dim1_d (d1), dim2_d (d2), dim3_d (d3), dim4_d(d4), 
    dim3_times_dim4_d(d3 * d4), dim2_times_dim3_times_dim4_d(d2 * d3 * d4), vector<T> (d1 * d2 * d3 * d4) { }
  T& operator()(const int i, const int j, const int k, const int l) 
    { 
      assert((0 <= i) && (i < dim1_d)); assert((0 <= j) && (j < dim2_d)); assert((0 <= k) && (k < dim3_d)); assert((0 <= l) && (l < dim4_d));
      return vector<T>::operator[](i * dim2_times_dim3_times_dim4_d + j * dim3_times_dim4_d + k * dim4_d + l); 
    }
  T operator()(const int i, const int j, const int k, const int l) const
    { 
      assert((0 <= i) && (i < dim1_d)); assert((0 <= j) && (j < dim2_d)); assert((0 <= k) && (k < dim3_d)); assert((0 <= l) && (l < dim4_d));
      return vector<T>::operator[](i * dim2_times_dim3_times_dim4_d + j * dim3_times_dim4_d + k * dim4_d + l); 
    }
  array_4d<T>& operator+=(const array_4d<T>& C)
  {
    assert(dim1() == C.dim1() &&
           dim2() == C.dim2() &&
           dim3() == C.dim3() &&
           dim4() == C.dim4());
    for (int i=0; i<C.size(); ++i)
      (*this)[i]+=C[i];
    return *this;
  }
  void Clear()
  {
    for (int i=0; i<this->size(); ++i)
      (*this)[i]=0.;
  }
  T& operator() (const vector<int>& indices)
  {
    assert(indices.size() == 4);
    return operator()(indices[0], indices[1], indices[2], indices[3]);
  }
  T operator() (const vector<int>& indices) const
  {
    assert(indices.size() == 4);
    return operator()(indices[0], indices[1], indices[2], indices[3]);
  }

  void resize (const int i, const int j, const int k, const int l) 
  { 
    vector<T>::resize (i * j * k * l); dim1_d = i; dim2_d = j; dim3_d = k; dim4_d = l; 
    dim2_times_dim3_times_dim4_d = dim2_d * dim3_d * dim4_d; dim3_times_dim4_d = dim3_d * dim4_d;
  }
  string reflect_type() const { return string("array_4d"); }

  int dim1() const { return dim1_d; }
  int dim2() const { return dim2_d; }
  int dim3() const { return dim3_d; }
  int dim4() const { return dim4_d; } 
 private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & dim1_d & dim2_d & dim3_d & dim4_d & dim2_times_dim3_times_dim4_d & dim3_times_dim4_d;
    ar & boost::serialization::base_object<std::vector<T> >(*this);
  }
  int dim1_d;
  int dim2_d;
  int dim3_d;
  int dim4_d;
  int dim2_times_dim3_times_dim4_d;
  int dim3_times_dim4_d;
};

#endif
