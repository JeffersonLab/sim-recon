// $Id$
//
//    File: JEventProcessor_FDC_online.h
// Created: Fri Nov  9 11:58:09 EST 2012
// Creator: wolin (on Linux stan.jlab.org 2.6.32-279.11.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_FDC_online_
#define _JEventProcessor_FDC_online_

#include <JANA/JEventProcessor.h>


class JEventProcessor_FDC_online:public jana::JEventProcessor{
 public:
  JEventProcessor_FDC_online();
  ~JEventProcessor_FDC_online();
  const char* className(void){return "JEventProcessor_FDC_online";}


 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.

  
  unsigned int thresh;
  float wire_pitch;
  float strip_pitch_u;
  float strip_pitch_d;
  float strip_angle;
  float cell_rot_step;
  
  float ADCmax[4][6][2][192][20];
  float ADCtime[4][6][2][192][20];
  int ADCnh[4][6][2][192];
  float TDCval[4][6][96][20];
  int TDCnh[4][6][96];
  
};

#endif // _JEventProcessor_FDC_online_

