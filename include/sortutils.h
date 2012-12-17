/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SORT_UTILS_HEADER_H
#define SORT_UTILS_HEADER_H
#include <vector>
#include <map>
template<class T> void sort_data_to_indices(const vector<T>& data, vector<int>& indices)
{
  multimap<const T, int> sorted_data;
  for (int i = 0; i < data.size(); ++i)
    sorted_data.insert(pair<const T, int>(data[i], i));
  typename multimap<const T, int>::iterator it = sorted_data.begin();
  indices.reserve(data.size());
  while (it != sorted_data.end())
    {
      indices.push_back(it->second);
      ++it;
    }
}
#endif// SORT_UTILS_HEADER_H
