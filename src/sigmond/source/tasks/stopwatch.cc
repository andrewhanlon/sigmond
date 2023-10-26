#include "stopwatch.h"
#include <sstream>
#include <stdexcept>
using namespace std;


StopWatch::StopWatch() 
{
 stoppedP=false;
 startedP=false;
 state=true;
}

StopWatch::~StopWatch() {}

void StopWatch::reset() 
{
 startedP = false;
 stoppedP = false;
 state=true;
}

void StopWatch::start() 
{
 int ret_val;
 ret_val = gettimeofday(&t_start, NULL);
 if( ret_val != 0 ){
    cerr << "Gettimeofday failed in StopWatch::start()" << endl;
    state=false;}
 else 
    state=true;
 startedP = true;
 stoppedP = false;
}

void StopWatch::stop() 
{
 if( !startedP ){ 
    throw(std::runtime_error("Attempting to stop a non running stopwatch in StopWatch::stop()"));}
 int ret_val;
 ret_val = gettimeofday(&t_end, NULL);
 if( ret_val != 0 ){
    cerr << "Gettimeofday failed in StopWatch::end()" << endl;
    state=false;}
 stoppedP = true;
}

double StopWatch::getTimeInMicroseconds() 
{
 if (!state){
    cerr << "Timer in error state"<<endl;
    return -1.0;}
 long usecs=0;
 if( startedP && stoppedP ){ 
    if( t_end.tv_sec < t_start.tv_sec ){ 
      cerr << "Critical timer rollover" << endl;
      state=false;
      return -1.0;}
    else{ 
      usecs = (t_end.tv_sec - t_start.tv_sec)*1000000;
      if( t_end.tv_usec < t_start.tv_usec ){
	usecs -= 1000000;
	usecs += 1000000+t_end.tv_usec - t_start.tv_usec;}
      else{
	usecs += t_end.tv_usec - t_start.tv_usec;}}}
  else{
    throw(std::runtime_error("Either stopwatch not started, or not stopped"));}
 return (double)usecs;
}
    
double StopWatch::getTimeInSeconds()  
{
 if (!state){
    cerr << "Timer in error state"<<endl;
    return -1.0;}
 long secs=0;
 long usecs=0;
 if( startedP && stoppedP ){ 
   if( t_end.tv_sec < t_start.tv_sec ){ 
     cerr << "Critical timer rollover" << endl;
     state=false;
     return -1.0;}
   else { 
     secs = t_end.tv_sec - t_start.tv_sec;
     if( t_end.tv_usec < t_start.tv_usec ) {
       secs -= 1;
       usecs = 1000000;}
     usecs += t_end.tv_usec - t_start.tv_usec;}}
 else{
   throw(std::runtime_error("Either stopwatch not started, or not stopped"));}
 return (double)secs + ((double)usecs / 1e6);
}

// ****************************************************

string currDateTimeString()
{
 time_t tim=time(NULL);
 tm *now=localtime(&tim);
 stringstream tmp;
 tmp << "D:"<<now->tm_mon+1<<":"<<now->tm_mday
     << "-T:"<<now->tm_hour<<":"<<now->tm_min<<":"<<now->tm_sec;
 return tmp.str();
}

// ****************************************************
