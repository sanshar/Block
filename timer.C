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


#include "global.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <errno.h>
#include "malloc.h"

namespace SpinAdapted{

void mcheck(const char* message)
{
  //#ifndef NDEBUG
   struct rusage usage;
   getrusage(RUSAGE_SELF, &usage);
   double totalmem = double(usage.ru_maxrss);
   double totalmem2 = 0.;
   
   //first get the pid number
   int procid= getpid();

   char file[500];
   sprintf(file, "/proc/%d/stat",procid);
   ifstream stat_stream(file,ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss;

  long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
  double rass = rss * page_size_kb;
  pout << "\t\t\t VM: " << vsize/(1024.0*1024.0) <<" Mb = "<< vsize/(1024.0*1024.0*1024.0) <<" Gb: rss = "<<rass/(1024.0*1024.0) <<" Gb "<<endl;
   
   //#endif
}

void mdebugcheck(const char* message)
{
  if (DEBUG_MEMORY)
    mcheck(message);
}
}
