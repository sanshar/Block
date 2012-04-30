#ifndef SPIN_TIMER_HEADER_H
#define SPIN_TIMER_HEADER_H
#include <ctime>
#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <cassert>
#ifndef SERIAL
#include <boost/mpi/timer.hpp>
#endif

using namespace std;
namespace SpinAdapted{
class cumulTimer
{
 public:

#ifndef SERIAL
  cumulTimer() : localStart(0), cumulativeSum(0){ t = boost::mpi::timer();};
  void start() {if(!omp_get_thread_num()) {localStart = t.elapsed();}}
#else
  cumulTimer() : localStart(0), cumulativeSum(0) {};
  void start() {localStart = clock();}
#endif

  void stop() 
  {
    if(!omp_get_thread_num()){
      if (localStart == 0) 
	{
	  cout << "local stop called without starting first"<<endl;
	  throw 20;
	  assert(1==2);
	  abort();
	}
#ifndef SERIAL
      cumulativeSum = cumulativeSum + t.elapsed() - localStart;
#else
      cumulativeSum = cumulativeSum + clock() - localStart;
#endif
      localStart = 0;
    }
  }

  friend ostream& operator<<(ostream& os, const cumulTimer& t)
  {
    os << t.cumulativeSum;
    return os;
  }

 private:
#ifndef SERIAL
  boost::mpi::timer t;
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
  void start() { t = boost::mpi::timer(); walltime = t.elapsed(); lastwalltime = walltime; wallstarttime = walltime; cputime = clock(); lastcputime = cputime; cpustarttime = cputime; }  
  double elapsedwalltime() { walltime = t.elapsed(); double elapsed = walltime - lastwalltime; lastwalltime = walltime; return elapsed; }
  double elapsedcputime() { cputime = clock(); double elapsed = double(cputime - lastcputime) / double(CLOCKS_PER_SEC); lastcputime = cputime; return elapsed; }
  double totalwalltime() { return t.elapsed() - wallstarttime; }
  double totalcputime() { return (double)(clock() - cpustarttime) / double(CLOCKS_PER_SEC); }
#else
  void start() { walltime = time(NULL); lastwalltime = walltime; wallstarttime = walltime; cputime = clock(); lastcputime = cputime; cpustarttime = cputime; }  
  long elapsedwalltime() { walltime = time(NULL); time_t elapsed = walltime - lastwalltime; lastwalltime = walltime; return elapsed; }
  double elapsedcputime() { cputime = clock(); double elapsed = double(cputime - lastcputime) / double(CLOCKS_PER_SEC); lastcputime = cputime; return elapsed; }
  long totalwalltime() { return time(NULL) - wallstarttime; }
  double totalcputime() { return (double)(clock() - cpustarttime) / double(CLOCKS_PER_SEC); }
#endif
private:
#ifndef SERIAL
  boost::mpi::timer t;
  double  walltime;
  double lastwalltime;
  double wallstarttime;
#else
  time_t walltime;
  time_t lastwalltime;
  time_t wallstarttime;
#endif
  clock_t cputime;
  clock_t lastcputime;
  clock_t cpustarttime;
};
void mcheck(const char* message);
void mdebugcheck(const char* message);
}
#endif
