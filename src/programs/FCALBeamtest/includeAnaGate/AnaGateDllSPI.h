// ----------------------------------------------------------------------------
// Projekt       : Analytica AnaGate API
// File          : AnaGateDllSPI.h
// Author        : Axel Schmidt
// Copyright     : (C) 2005 by Analytica GmbH
// ----------------------------------------------------------------------------
// $Id: AnaGateDllSPI.h,v 1.1.1.1 2010/07/16 18:12:17 ctarbert Exp $
//
// $Log: AnaGateDllSPI.h,v $
// Revision 1.1.1.1  2010/07/16 18:12:17  ctarbert
//
//
// Revision 1.4  2009/04/30 13:24:56  stefanwelisch
// Referenzen als Pointer implementiert
//
// Revision 1.3  2007/05/21 11:43:39  axelschmidt
// standard timeout 500 -> 1000 ms
//
// Revision 1.2  2005/12/27 08:26:37  AxelSchmidt
// base functions of SPI interface
//
// Revision 1.1  2005/11/28 17:05:48  axelschmidt
// SPI integrated, initial version
//
// ----------------------------------------------------------------------------

#ifndef ANAGATE_DLL_SPI_H
#define ANAGATE_DLL_SPI_H

// Defines --------------------------------------------------------------------
#include "AnaGateDLL.h"

// Prototyping ----------------------------------------------------------------
extern "C"
{
   /** Opens an AnaGate SPI device.
      @param pHandle Pointer to an integer, in which the device handle is stored, if device is
                     opened successfully.
      @param pcIPAddress Tcp/ip address of the AnaGate device (Port is always 5000).
      @param nTimeout Standard tcp/ip timeout in millseconds.
      @return If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int  SPIOpenDevice( int *pHandle, const char * pcIPAddress, int nTimeout = 1000 );

   /** Closes an open AnaGate SPI device.
       @param hHandle Device handle (from a successfull #OpenDevice call).
       @return If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int  SPICloseDevice( int hHandle );

   /** Writes data to the SPI partner The data is routed to the partner directly without any transformation,
      so the user of the interface has to build a correct data buffer.
      @param hHandle Device handle (from a successfull #OpenDevice call).
      @param pcBufWrite    Pointer to the buffer containing the data to be written.
      @param nBufWriteLen  Number of bytes to be written.
      @param pcBufRead     Pointer to the buffer that receives the data read.
      @param nBufReadLen   Number of bytes to be read (should be same as number of bytes written).
      @return. If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int  SPIDataReq( int hHandle, const char * pcBufWrite, int nBufWriteLen,
                                                char * pcBufRead, int nBufReadLen );

   /** Reads data from the specified offset to a data buffer.
      @remarks This command is a specialized version of the #Read method and is optimized for SPI EEProms.
               <p>This version implicit performs a write request to set the offset where to read from and the
               following read request in a single tcpip telegramm.</p>
               <p>Use this method only for EEProm partners</p>
      @param hHandle Device handle (from a successfull #OpenDevice call).
      @param nChipAddress  Chip enable address (bit 1 - 3 of the slave address).
      @param nOffset     Offset position where data is to be read (EEProm partner).
      @param pcBuffer    Pointer to the buffer that receives the data read.
      @param nBufferLen  Number of bytes to be read.
      @param nOffsetFormat Offset format. To create a I2C command, the method must specify the offset in
                           the correct format. At this time two formats are supported (8-bit and 16bit addresses).
      @return. If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int  SPIErrorMessage( int nRC, char *pcMessage, int nMessageLen );

   /** Sets the global settings of an AnaGate I2C device.
       @param hHandle Device handle (from a successfull #OpenDevice call).
       @param nBaudrate Baud rate.
       @param nSigLevel Signal level.
       @param nAuxVoltage Auxiliary voltage supply.
       @param nClockMode The phase and the polarity of the clock line.
       @return If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int  SPISetGlobals( int hHandle, int nBaudrate, unsigned char nSigLevel,
                                                   unsigned char nAuxVoltage, unsigned char nClockMode );
   /** Gets the global settings of an open AnaGate I2C device.
       @param hHandle Device handle (from a successfull #OpenDevice call).
       @param nBaudrate Pointer to the variable that receives Baud rate.
       @param nSigLevel Pointer to the variable that receives Signal level.
       @param nAuxVoltage Pointer to the variable that receives Auxiliary voltage supply.
       @param nClockMode Pointer to the variable that receives The phase and the polarity of the clock line.
       @return If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int  SPIGetGlobals( int hHandle, int * pnBaudrate, unsigned char * pnSigLevel,
                                                   unsigned char * pnAuxVoltage, unsigned char * pnClockMode );

   /** Reads data from digital io register of the partner.
      @param hHandle Device handle (from a successfull #OpenDevice call).
      @param nInputBits   Pointer to the variable that receives the digtial input register read.
      @param nOutputBits  Pointer to the variable that receives the digtial output register read.
      @return. If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int SPIReadDigital( int hHandle, unsigned long * pnInputBits, unsigned long * pnOutputBits );

   /** Writes data to digital io register of the partner.
      @param hHandle Device handle (from a successfull #OpenDevice call).
      @param nOutputBits  Refererence to the variable that hold the digtial IO output register bits to write.
      @return. If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int SPIWriteDigital( int hHandle, unsigned long nOutputBits );
}

#endif // ANAGATE_DLL_I2C_H
