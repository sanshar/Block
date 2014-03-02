#ifndef _BTAS_CXX11_TIME_STAMP_H
#define _BTAS_CXX11_TIME_STAMP_H 1

#include <chrono> // C++11 clock classes

//! Clock function for timing
class time_stamp {
public:
  typedef std::chrono::time_point<std::chrono::system_clock, std::chrono::microseconds> time_point_type;
private:
  //! Starting time point
  time_point_type
    m_time_point_start;
  //! Lap time record
  time_point_type
    m_time_point_lap;
public:
  //! Default constructor
  time_stamp() { start(); }
  //! Start / Restart
  void start() {
    m_time_point_start = std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::system_clock::now());
    m_time_point_lap   = std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::system_clock::now());
  }
  //! Get elapased time in microseconds
  double elapsed(unsigned long n_periods = 1000000 /* to seconds */) const {
    time_point_type record = std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::system_clock::now());
    auto t = record - m_time_point_start;
    return static_cast<double>(t.count())/n_periods;
  }
  //! Get lap time in microseconds, and reset lap time record
  double lap(unsigned long n_periods = 1000000 /* to seconds */) {
    time_point_type record = std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::system_clock::now());
    auto t = record - m_time_point_lap;
    m_time_point_lap = record;
    return static_cast<double>(t.count())/n_periods;
  }
};

#endif // _BTAS_CXX11_TIME_STAMP_H
