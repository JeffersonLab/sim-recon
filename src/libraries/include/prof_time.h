// $Id$
//
//    File: prof_time.h
// Created: Tue Jan 11 15:13:48 EST 2011
// Creator: davidl (on Darwin eleanor.jlab.org 10.6.0 i386)
//

#ifndef _prof_time_
#define _prof_time_

#include <sys/time.h>
#include <map>
#include <string>

/// This utility class is to help in doing detailed profiling of code.
/// It is uses the unix itimer system (see man getitimer). It is assumed
/// that all 3 of the hi-res timers are running. They should be because
/// JANA starts them all at program start up.
///
/// To use this, one would want to do something like the following:
///
/// Declare a container to hold timing info outside of evnt loop
///
///  static map<string, prof_time::time_diffs> prof_times;
///
/// Then, wrap the block of code you want to time with a prof_time
/// variable and a call to its TimeDiffNow(...) method.
///
///  prof_time start_time;
///
///   <... do something ...>
///
///  start_time.TimeDiffNow(prof_times, "My code block 1");
///
/// Note that the current time is copied into the prof_time variable
/// upon it's creation. To force the recording of a time at a
/// position other than the line creating the object, call it's GetTime()
/// method.
///
/// This system is designed to check all 3 timers. This might add a little
/// more overhead than simply checking one timer, but at least this way
/// this can be used for a broad range of purposes.

class prof_time{
	public:
	
		// Class for holding all times as doubles
		class time_diffs {
			public:
				double real;
				double virt;
				double prof;
				
				time_diffs():real(0),virt(0),prof(0){}
		};
	
		
		//----------------------
		// ctor
		prof_time(){GetTime();}
		
		//----------------------
		// dtor
		virtual ~prof_time(){}

		//----------------------
		// GetTime
		void GetTime(void){
			getitimer(ITIMER_REAL, &real_tmr);
			getitimer(ITIMER_VIRTUAL, &virt_tmr);
			getitimer(ITIMER_PROF, &prof_tmr);
		}

		//----------------------
		// TimeDiffNow
		void TimeDiffNow(std::map<std::string, time_diffs> &prof_times, std::string what){ // update(add to) entry in map or create new one
			time_diffs tdiffs;
			prof_time now;
			TimeDiff(now, tdiffs.real, tdiffs.virt, tdiffs.prof);

			std::map<string, time_diffs>::iterator iter = prof_times.find(what);
			if(iter==prof_times.end()){
				prof_times[what] = tdiffs;
			}else{
				(iter->second).real += tdiffs.real;
				(iter->second).virt += tdiffs.virt;
				(iter->second).prof += tdiffs.prof;
			}
		}

		//----------------------
		// TimeDiff
		void TimeDiff(const prof_time &other_time, double &tdiff_real, double &tdiff_virt, double &tdiff_prof){
			struct itimerval tmp_real;
			struct itimerval tmp_virt;
			struct itimerval tmp_prof;
			
			// Automatically decide which to subtract from which depending
			// on which is greater. Assume if REAL timer is >, then other
			// timers follow suit.
			if(timercmp(&real_tmr.it_value, &other_time.real_tmr.it_value, >)){
				timersub(&real_tmr.it_value, &other_time.real_tmr.it_value, &tmp_real.it_value);
				timersub(&virt_tmr.it_value, &other_time.virt_tmr.it_value, &tmp_virt.it_value);
				timersub(&prof_tmr.it_value, &other_time.prof_tmr.it_value, &tmp_prof.it_value);
			}else{
				timersub(&other_time.real_tmr.it_value, &real_tmr.it_value, &tmp_real.it_value);
				timersub(&other_time.virt_tmr.it_value, &virt_tmr.it_value, &tmp_virt.it_value);
				timersub(&other_time.prof_tmr.it_value, &prof_tmr.it_value, &tmp_prof.it_value);
			}
			tdiff_real = (double)tmp_real.it_value.tv_sec + (double)tmp_real.it_value.tv_usec/1.0E6;
			tdiff_virt = (double)tmp_virt.it_value.tv_sec + (double)tmp_virt.it_value.tv_usec/1.0E6;
			tdiff_prof = (double)tmp_prof.it_value.tv_sec + (double)tmp_prof.it_value.tv_usec/1.0E6;
		}
		
		struct itimerval real_tmr;
		struct itimerval virt_tmr;
		struct itimerval prof_tmr;

};

#endif // _prof_time_

