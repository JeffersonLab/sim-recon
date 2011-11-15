// ----------------------------------------------------------------------------
// Projekt       : Analytica AnaGate Renesas DLL
// File          : AnaGateDllRenesas.h
// Author        : Axel Schmidt
// Copyright     : (C) 2007 by Analytica GmbH
// ----------------------------------------------------------------------------
// $Id: AnaGateDllRenesas.h,v 1.1.1.1 2010/07/16 18:12:17 ctarbert Exp $
//
// $Log: AnaGateDllRenesas.h,v $
// Revision 1.1.1.1  2010/07/16 18:12:17  ctarbert
//
//
// Revision 1.3  2009/04/30 13:24:56  stefanwelisch
// Referenzen als Pointer implementiert
//
// Revision 1.2  2007/10/25 07:02:35  axelschmidt
// Block-Erase added
//
// Revision 1.1  2007/05/21 11:41:00  axelschmidt
// initial
//
// ----------------------------------------------------------------------------

#ifndef ANAGATE_DLL_RENESAS_H
#define ANAGATE_DLL_RENESAS_H

// Defines --------------------------------------------------------------------
// Includes--------------------------------------------------------------------

#include "AnaGateDLL.h"

// Prototyping ----------------------------------------------------------------
extern "C"
{
   /** Opens an AnaGate Renesas device.
      @param pHandle Pointer to an integer, in which the device handle is stored, if device is
                     opened successfully.
      @param pcIPAddress Tcp/ip address of the AnaGate device.
      @param nTimeout Standard tcp/ip timeout in millseconds.
      @return If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int  RenesasOpenDevice( int *pHandle, const char * pcIPAddress, int nTimeout = 1000 );

   /** Closes an open AnaGate Renesas device.
       @param hHandle Device handle (from a successfull #OpenDevice call).
       @return If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int  RenesasCloseDevice( int hHandle );

   ANAGATEDLL_API int  RenesasInit( int hHandle );
   ANAGATEDLL_API int  RenesasCheckID( int hHandle, const char *pcIDBuffer, int nBufLen );
   ANAGATEDLL_API int  RenesasComplete( int hHandle );
   ANAGATEDLL_API int  RenesasEraseAll( int hHandle );
   ANAGATEDLL_API int  RenesasBlockErase( int hHandle, unsigned short nPageAddress );
   ANAGATEDLL_API int  RenesasGetVersion( int hHandle, char *pcBuffer, int nBufLen );
   ANAGATEDLL_API int  RenesasSetBaudrate( int hHandle, int nBaudrate );
   ANAGATEDLL_API int  RenesasReadPage( int hHandle, unsigned short nPageAddress, char *pcBuffer, int nNumberOfBytesToRead, int * pnNumberOfBytesRead  );
   ANAGATEDLL_API int  RenesasWritePage ( int hHandle, unsigned short nPageAddress, const char *pcBuffer, int nNumberOfBytesToWrite );
   ANAGATEDLL_API int  RenesasVerifyCheck ( int hHandle, unsigned short nStartAddress, unsigned short nEndAddress, unsigned short * pnChecksum);

   /** Retrieves a textual error description of the supplied return code.
      @param nRC Return code.
      @param pcMessage Pointer to a c-style character buffer, in which the retrieved error string is stored.
      @param nMessageLen Length of the supplied charcter buffer. If the error string does not fit into the
             buffer, the string is shortened.
      @return The byte of the returned error string.
   */
   ANAGATEDLL_API int  RenesasErrorMessage( int nRC, char *pcMessage, int nMessageLen );

   /** Reads data from digital io register of the partner.
      @param hHandle Device handle (from a successfull #OpenDevice call).
      @param nInputBits   Pointer to the variable that receives the digtial input register read.
      @param nOutputBits  Pointer to the variable that receives the digtial output register read.
      @return. If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int RenesasReadDigital( int hHandle, unsigned long * pnInputBits, unsigned long * pnOutputBits );

   /** Writes data to digital io register of the partner.
      @param hHandle Device handle (from a successfull #OpenDevice call).
      @param nOutputBits  Refererence to the variable that hold the digtial IO output register bits to write.
      @return. If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int RenesasWriteDigital( int hHandle, unsigned long nOutputBits );

}

#endif // ANAGATE_DLL_RENESAS_H
