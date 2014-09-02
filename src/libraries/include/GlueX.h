// $Id$
//
//    File: GlueX.h
// Created: Tue Aug 23 04:47:57 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
// Modified: yqiang, Oct 10 2012, add RICH
//

#ifndef _GlueX_
#define _GlueX_

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
     SYS_RICH       = 0x0200,
     SYS_CCAL       = 0x0400,
     SYS_TAGH       = 0x0800
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
          case SYS_RICH:
              return "RICH";
              break;
          case SYS_CCAL:
              return "CCAL";
              break;
     }
     return "UNKNOWN";
}

#endif // _GlueX_
