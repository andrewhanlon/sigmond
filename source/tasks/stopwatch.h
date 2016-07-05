#ifndef STOPWATCH_H
#define STOPWATCH_H
#include <sys/time.h>
#include <iostream>
#include <cstdlib>
#include <string>

  // ***************************************************************
  // *                                                             *
  // *    Useful timing objects. The "start" member begins a       *
  // *    timing, and "stop" ends the timing.  "getTimeInSeconds"  *
  // *    or "getTimeInMicroSeconds" is then used to output the    *
  // *    time between the last "start" and "stop".  Calling       *
  // *    "start" again begins a new timing (not cumulative).      *
  // *                                                             *
  // *    Usage:                                                   *
  // *                                                             *
  // *      StopWatch rolex;                                       *
  // *      rolex.start(); ...statements ... rolex.stop();         *
  // *      cout << rolex.getTimeInSeconds()<<endl;                *
  // *                                                             *
  // ***************************************************************

class StopWatch 
{
 public:
   StopWatch();
   ~StopWatch();
   void reset();
   void start();
   void stop();
   double getTimeInMicroseconds();
   double getTimeInSeconds();

 private:
   long sec;
   long usec;
   bool startedP;
   bool stoppedP;
   bool state;

  struct timeval t_start;
  struct timeval t_end;
};


  // ***************************************************************
  // *                                                             *
  // *   Returns a string showing the current date (without year)  *
  // *   and time, of the form                                     *
  // *                                                             *
  // *            D:12:23-T:17:10:58                               *
  // *            D:month:day-T:hour:min:sec                       *
  // *                                                             *
  // *   This can be useful for making a unique MC observable      *
  // *   name for a temporary data entry.                          *
  // *                                                             *
  // ***************************************************************
    
std::string currDateTimeString();


  // ***************************************************************
#endif
