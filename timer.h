/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_TIMER_HEADER_H
#define SPIN_TIMER_HEADER_H
#include <ctime>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <execinfo.h>
#include <stdlib.h>
#include <stdio.h>

#ifndef SERIAL
#include <boost/serialization/shared_ptr.hpp>
#include <boost/mpi/timer.hpp>
#include <include/communicate.h>
#endif

using namespace std;
namespace SpinAdapted{
class cumulTimer
{
 public:

#ifndef SERIAL
 cumulTimer() : localStart(-1.0), cumulativeSum(0){t = boost::shared_ptr<boost::mpi::timer> (new boost::mpi::timer());};
#ifdef _OPENMP
  void start() {if(!omp_get_thread_num()) {localStart = t->elapsed();}}
#else
  void start() {localStart = t->elapsed();}
#endif
#else
  cumulTimer() : localStart(-1.0), cumulativeSum(0) {};
  void start() {localStart = clock();}
#endif

  void stop() 
  {
#ifdef _OPENMP
    if(!omp_get_thread_num()){
#endif

#ifndef SERIAL
      if (localStart < 0 || localStart > t->elapsed() +1 ) 
	{
	  
	  cout << "local stop called without starting first"<<endl;
	  cout << localStart<<"  "<<t->elapsed();
	  throw 20;
	  assert(1==2);
	  abort();
	}

      cumulativeSum = cumulativeSum + t->elapsed() - localStart;
#else
      if (localStart < 0 || localStart > clock() +1 ) 
	{
	  
	  cout << "local stop called without starting first"<<endl;
	  cout << localStart<<"  "<<clock();
	  throw 20;
	  assert(1==2);
	  abort();
	}

      cumulativeSum = cumulativeSum + clock() - localStart;
#endif
      localStart = 0;
#ifdef _OPENMP
    }
#endif
  }

  friend ostream& operator<<(ostream& os, const cumulTimer& t)
  {
#ifndef SERIAL
    os << t.cumulativeSum;
#else
    os << ((float)t.cumulativeSum)/CLOCKS_PER_SEC;
#endif
    return os;
  }

 private:
#ifndef SERIAL
  boost::shared_ptr<boost::mpi::timer> t;
#endif
  double localStart;
  double cumulativeSum;
};

class Timer
{
public:
  Timer(bool s) { if (s) start(); }
  Timer() { start(); }
#ifndef SERIAL
  void start() { t = boost::shared_ptr<boost::mpi::timer>(new boost::mpi::timer()); walltime = t->elapsed(); lastwalltime = walltime; wallstarttime = walltime; cputime = clock(); lastcputime = cputime; cpustarttime = cputime; }  
  double elapsedwalltime() { walltime = t->elapsed(); double elapsed = walltime - lastwalltime; lastwalltime = walltime; return elapsed; }
  double elapsedcputime() { cputime = clock(); double elapsed = double(cputime - lastcputime) / double(CLOCKS_PER_SEC); elapsed=accumclock(elapsed);lastcputime = cputime; return elapsed; }
  double totalwalltime() { return t->elapsed() - wallstarttime; }
  double totalcputime() { double total=double(clock() - cpustarttime) / double(CLOCKS_PER_SEC);total=accumclock(total);return total; }
#else
  void start() { walltime = time(NULL); lastwalltime = walltime; wallstarttime = walltime; cputime = clock(); lastcputime = cputime; cpustarttime = cputime; }  
  long elapsedwalltime() { walltime = time(NULL); time_t elapsed = walltime - lastwalltime; lastwalltime = walltime; return elapsed; }
  double elapsedcputime() { cputime = clock(); double elapsed = double(cputime - lastcputime) / double(CLOCKS_PER_SEC); lastcputime = cputime; return elapsed; }
  long totalwalltime() { return time(NULL) - wallstarttime; }
  double totalcputime() { return (double)(clock() - cpustarttime) / double(CLOCKS_PER_SEC); }
#endif
private:
#ifndef SERIAL
  boost::shared_ptr<boost::mpi::timer> t;
  double  walltime;
  double lastwalltime;
  double wallstarttime;
  double accumclock(double tcpu) 
  {
   boost::mpi::communicator world;
   double sum;
   all_reduce(world, tcpu, sum, std::plus<double>());
   return sum;
  }
#else
  time_t walltime;
  time_t lastwalltime;
  time_t wallstarttime;
#endif
  clock_t cputime;
  clock_t lastcputime;
  clock_t cpustarttime;
};

void __GetMachineName(char* machineName);
void mcheck(const char* message);
void mdebugcheck(const char* message);
}
#endif
