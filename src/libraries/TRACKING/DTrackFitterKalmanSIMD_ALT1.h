// $Id$
//
//    File: DTrackFitterKalmanSIMD_ALT1.h
// Created: Tue Mar 29 09:45:14 EDT 2011
// Creator: staylor (on Linux ifarml1 2.6.18-128.el5 x86_64)
//

#ifndef _DTrackFitterKalmanSIMD_ALT1_
#define _DTrackFitterKalmanSIMD_ALT1_

#include <JANA/jerror.h>
#include <TRACKING/DTrackFitterKalmanSIMD.h>

class DTrackFitterKalmanSIMD_ALT1: public DTrackFitterKalmanSIMD{
 public:
  DTrackFitterKalmanSIMD_ALT1(JEventLoop *loop):DTrackFitterKalmanSIMD(loop){};
    //DTrackFitterKalmanSIMD_ALT1();
    virtual ~DTrackFitterKalmanSIMD_ALT1(){};
    
  jerror_t KalmanForward(double anneal,DMatrix5x1 &S,DMatrix5x5 &C,
			 double &chisq,unsigned int &numdof);
  
  // Virtual methods from TrackFitter base class
  string Name(void) const {return string("KalmanSIMD_ALT1");}
 protected:
	
  
 private:
  
};

#endif // _DTrackFitterKalmanSIMD_ALT1_

