// ----------------------------------------------------------------------------
// Projekt       : Analytica AnaGate API
// File          : AnaGateDllDIO.h
// Author        : Axel Schmidt
// Copyright     : (C) 2005 by Analytica GmbH
// ----------------------------------------------------------------------------
// $Id: AnaGateDllDIO.h,v 1.1.1.1 2010/07/16 18:12:17 ctarbert Exp $
//
// $Log: AnaGateDllDIO.h,v $
// Revision 1.1.1.1  2010/07/16 18:12:17  ctarbert
//
//
// Revision 1.4  2009/04/30 13:24:56  stefanwelisch
// Referenzen als Pointer implementiert
//
// Revision 1.3  2007/05/21 11:43:39  axelschmidt
// standard timeout 500 -> 1000 ms
//
// Revision 1.2  2005/08/04 11:02:21  axelschmidt
// reads also output register
//
// Revision 1.1  2005/08/04 09:52:05  axelschmidt
// intial
//
// ----------------------------------------------------------------------------

#ifndef ANAGATE_DLL_DIO_H
#define ANAGATE_DLL_DIO_H

// Defines --------------------------------------------------------------------
#ifndef LUA_LIB
   #ifdef WIN32
      #ifdef ANADIODLL_EXPORTS
         #define ANADIODLL_API __declspec(dllexport)
      #else
         #define ANADIODLL_API __declspec(dllimport)
      #endif
   #else
      #define ANADIODLL_API
   #endif
#else
   #define ANADIODLL_API
#endif

#include "AnaGateDLL.h"

// Prototyping ----------------------------------------------------------------
extern "C"
{
   /** Opens an AnaGate DIO device.
      @param pHandle Pointer to an integer, in which the device handle is stored, if device is
                     opened successfully.
      @param pcIPAddress Tcp/ip address of the AnaGate device (Port is always 5000).
      @param nTimeout Standard tcp/ip timeout in millseconds.
      @return If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANADIODLL_API int  DIOOpenDevice( int *pHandle, const char * pcIPAddress, int nTimeout = 1000 );

   /** Closes an open AnaGate DIO device.
       @param hHandle Device handle (from a successfull #OpenDevice call).
       @return If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANADIODLL_API int  DIOCloseDevice( int hHandle );

   /** Reads data bits from the DIO partner.
      @param hHandle  Device handle (from a successfull #OpenDevice call).
      @param nInBits  Pointer to the variable that receives the digtial input register read.
      @param nOutBits Pointer to the variable that receives the digtial output register read.
      @return. If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANADIODLL_API int  DIORead ( int hHandle, unsigned long * pnInBits, unsigned long * pnOutBits );

   /** Writes data bits to the DIO partner.
      @param hHandle Device handle (from a successfull #OpenDevice call).
      @param nOutBits Holds the digital output register bits to write.
      @return. If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANADIODLL_API int  DIOWrite ( int hHandle, unsigned long nOutBits );

   /** Retrieves a textual error description of the supplied return code.
      @param nRC Return code.
      @param pcMessage Pointer to a c-style character buffer, in which the retrieved error string is stored.
      @param nMessageLen Length of the supplied charcter buffer. If the error string does not fit into the
             buffer, the string is shortened.
      @return The byte of the returned error string.
   */
   ANADIODLL_API int  DIOErrorMessage( int nRC, char *pcMessage, int nMessageLen );
}

#endif // ANAGATE_DLL_DIO_H
