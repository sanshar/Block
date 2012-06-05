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


#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include <pario.h>
#include <cstdlib>

using namespace std;
int callDmrg(char*, char*);

int main(int argc, char* argv []) 
{
#ifndef SERIAL
  boost::mpi::environment env(argc, argv);
#endif
  if (argc == 1) {
    pout << "No command line argument found, expects an input file"<<endl;
    abort();
  }
  char* output = 0;
  return callDmrg(argv[1], output);
}
