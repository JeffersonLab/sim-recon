/*
 * crc.h
 *
 *
 * CRC computation routines derived from Chuck Forsberg's Zmodem
 * implementation for UNIX. It is public domain, isn't it?
 *
 * -rev 04-16-87
 *  This file contains Unix specific code for setting terminal modes,
 *  very little is specific to ZMODEM or YMODEM per se (that code is in
 *  sz.c and rz.c).  The CRC-16 routines used by XMODEM, YMODEM, and ZMODEM
 *  are also in this file, a fast table driven macro version
 *
 */

#ifndef _CRCH_INCLUDED
#define _CRCH_INCLUDED

#define CRC16INIT 0
#define CRC32INIT (0xFFFFFFFFL)

#define CRC16MAGIC (0xFFFF)
#define CRC32MAGIC (0xDEBB20E3L)

#include <ntypes.h>

uint32 data_crc(const void*data,int length);

#endif
/* End of crc.h */
