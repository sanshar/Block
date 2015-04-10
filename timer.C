/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/
#include <unistd.h>
#include "global.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <errno.h>
#if defined(__APPLE__)
   #include <sys/malloc.h>
#else
   #include <malloc.h>
#endif
#include "pario.h"

namespace SpinAdapted{

void __GetMachineName(char* machineName)
{
  char Name[150];
  int i=0;
  
  gethostname(Name, 150);
  strncpy(machineName,Name, 150);
}

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
  p3out << "\t\t\t VM: " << vsize/(1024.0*1024.0) <<" Mb = "<< vsize/(1024.0*1024.0*1024.0) <<" Gb: rss = "<<rass/(1024.0*1024.0) <<" Gb "<<endl;
   
   //#endif
}

void mdebugcheck(const char* message)
{
  if (DEBUG_MEMORY)
    mcheck(message);
}
}
