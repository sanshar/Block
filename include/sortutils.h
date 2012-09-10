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
