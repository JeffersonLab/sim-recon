

#ifndef _HDDVIEW_H_
#define _HDDVIEW_H_

#include <iostream>
#include <iomanip>
using namespace std;

#include <TApplication.h>
#include <TCanvas.h>

#include "DEventLoop.h"
#include "derror.h"
#include "hddm_s.h"
#include "MyProcessor.h"

extern TCanvas *maincanvas;
extern DEventLoop *eventloop;
extern MyProcessor *myproc;

derror_t hdv_getevent(void);
derror_t hdv_drawevent(void);
void RotatePoints(int how);
void Orbit(void);

#endif //_HDDVIEW_H_
