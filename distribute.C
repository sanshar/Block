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
*/

#include "distribute.h"

namespace SpinAdapted{
#ifdef SERIAL
int receivefrom(int offsetproc)
{
  return 0;
}
#else
#include <boost/mpi/communicator.hpp>

int receivefrom(int offsetproc)
{
  boost::mpi::communicator world;
  int rank = world.rank();
  int size = world.size();
  rank -= offsetproc;

  
  if (rank < 0) rank += size; // modulo size arithmetic

  
  return ((rank - 1) / Broadcastsettings::npyramid + offsetproc) % size;
}
#endif

#ifdef SERIAL
void makesendlist(vector<int>& tolist, int offsetproc) {;}
#else
void makesendlist(vector<int>& tolist, int offsetproc)
{
  boost::mpi::communicator world;
  int rank = world.rank();
  int size = world.size();
  rank -= offsetproc; // if offsetproc == rank, then treat current proc as the root node
  if (rank < 0) rank += size;
  for (int i = 1; i <= Broadcastsettings::npyramid; ++i)
    {
      int tosend = rank * Broadcastsettings::npyramid +i;      
      if (tosend < size)
	tolist.push_back((tosend + offsetproc) % size); // if offsetproc is not zero, add it back on, and then wraparound
    }
}
#endif

}
