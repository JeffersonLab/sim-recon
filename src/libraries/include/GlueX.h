// $Id$
//
//    File: GlueX.h
// Created: Tue Aug 23 04:47:57 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
// Modified: yqiang, Oct 10 2012, add RICH
// Modified: jrsteven, June 22 2015, move RICH -> DIRC and remove CERE
//

#ifndef _GlueX_
#define _GlueX_

#include <string.h>

enum DetectorSystem_t{
     SYS_NULL       = 0x0000,
     SYS_CDC        = 0x0001,
     SYS_FDC        = 0x0002,
     SYS_BCAL       = 0x0004,
     SYS_TOF        = 0x0008,
     SYS_CHERENKOV  = 0x0010,
     SYS_FCAL       = 0x0020,
     SYS_UPV        = 0x0040,
     SYS_TAGM       = 0x0080,
     SYS_START      = 0x0100,
     SYS_DIRC       = 0x0200,
     SYS_CCAL       = 0x0400,
     SYS_TAGH       = 0x0800,
     SYS_RF         = 0x1000,
     SYS_PS         = 0x2000,
     SYS_PSC        = 0x4000,
     SYS_FMWPC      = 0x8000,
     SYS_TPOL       = 0x10000
};

inline const char* SystemName(DetectorSystem_t sys)
{
     switch(sys){
          case SYS_NULL:
              return "NULL_DETECTOR";
              break;
          case SYS_CDC:
              return "CDC";
              break;
          case SYS_FDC:
              return "FDC";
              break;
          case SYS_BCAL:
              return "BCAL";
              break;
          case SYS_TOF:
              return "TOF";
              break;
          case SYS_CHERENKOV:
              return "Cherenkov";
              break;
          case SYS_FCAL:
              return "FCAL";
              break;
          case SYS_UPV:
              return "UPV";
              break;
          case SYS_TAGM:
              return "TAGM";
              break;
          case SYS_TAGH:
              return "TAGH";
              break;
          case SYS_START:
              return "ST";
              break;
          case SYS_DIRC:
              return "DIRC";
              break;
          case SYS_CCAL:
              return "CCAL";
              break;
          case SYS_RF:
              return "RF";
              break;
          case SYS_PS:
              return "PS";
              break;
          case SYS_PSC:
              return "PSC";
              break;
          case SYS_FMWPC:
              return "FMWPC";
              break;
          case SYS_TPOL:
              return "TPOL";
              break;
     }
     return "UNKNOWN";
}

inline DetectorSystem_t NameToSystem(const char* locSystemName)
{
	if(strcmp(locSystemName, "CDC") == 0)
		return SYS_CDC;
	else if(strcmp(locSystemName, "FDC") == 0)
		return SYS_FDC;
	else if(strcmp(locSystemName, "BCAL") == 0)
		return SYS_BCAL;
	else if(strcmp(locSystemName, "TOF") == 0)
		return SYS_TOF;
	else if(strcmp(locSystemName, "Cherenkov") == 0)
		return SYS_CHERENKOV;
	else if(strcmp(locSystemName, "FCAL") == 0)
		return SYS_FCAL;
	else if(strcmp(locSystemName, "UPV") == 0)
		return SYS_UPV;
	else if(strcmp(locSystemName, "TAGM") == 0)
		return SYS_TAGM;
	else if(strcmp(locSystemName, "TAGH") == 0)
		return SYS_TAGH;
	else if(strcmp(locSystemName, "ST") == 0)
		return SYS_START;
	else if(strcmp(locSystemName, "SC") == 0)
		return SYS_START;
	else if(strcmp(locSystemName, "START") == 0)
		return SYS_START;
	else if(strcmp(locSystemName, "DIRC") == 0)
		return SYS_DIRC;
	else if(strcmp(locSystemName, "CCAL") == 0)
		return SYS_CCAL;
	else if(strcmp(locSystemName, "RF") == 0)
		return SYS_RF;
	else if(strcmp(locSystemName, "PS") == 0)
		return SYS_PS;
	else if(strcmp(locSystemName, "PSC") == 0)
		return SYS_PSC;
	else if(strcmp(locSystemName, "FMWPC") == 0)
		return SYS_FMWPC;
	else if(strcmp(locSystemName, "TPOL") == 0)
		return SYS_TPOL;
	else
		return SYS_NULL;
}

#endif // _GlueX_
