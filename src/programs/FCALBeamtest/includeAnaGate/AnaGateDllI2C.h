// ----------------------------------------------------------------------------
// Projekt       : Analytica AnaGate API
// File          : AnaGateDllI2C.h
// Author        : Axel Schmidt
// Copyright     : (C) 2004 by Analytica GmbH
// ----------------------------------------------------------------------------
// $Id: AnaGateDllI2C.h,v 1.1.1.1 2010/07/16 18:12:17 ctarbert Exp $
//
// $Log: AnaGateDllI2C.h,v $
// Revision 1.1.1.1  2010/07/16 18:12:17  ctarbert
//
//
// Revision 1.12  2009/04/30 13:24:56  stefanwelisch
// Referenzen als Pointer implementiert
//
// Revision 1.11  2007/05/21 11:43:38  axelschmidt
// standard timeout 500 -> 1000 ms
//
// Revision 1.10  2006/03/29 10:12:04  axelschmidt
// comment
//
// Revision 1.9  2005/11/30 16:16:04  axelschmidt
// Sequence command implemented
//
// Revision 1.8  2005/11/28 17:07:47  axelschmidt
// errorbyte parameter for I2CWrite function
//
// Revision 1.7  2005/09/16 09:18:48  axelschmidt
// Digital IO support for I2C and CAN
//
// Revision 1.6  2004/08/13 12:52:33  AxelSchmidt
// import/export-defines changed
//
// Revision 1.5  2004/08/09 15:48:17  AxelSchmidt
// div. Parameter von signed auf unsigned umgesetzt
//
// Revision 1.4  2004/08/02 07:07:15  AxelSchmidt
// Schnittstelle zur Klasse CAnaGateClientI2C geändert
//
// Revision 1.3  2004/07/23 15:53:49  AxelSchmidt
// new method for device Versions
//
// Revision 1.2  2004/07/21 16:13:33  AxelSchmidt
// new I2COpenDevice and I2CCloseDevice
//
// Revision 1.1  2004/07/21 06:53:41  AxelSchmidt
// initial version
//
// ----------------------------------------------------------------------------

#ifndef ANAGATE_DLL_I2C_H
#define ANAGATE_DLL_I2C_H

// Defines --------------------------------------------------------------------
#include "AnaGateDLL.h"

// Prototyping ----------------------------------------------------------------
extern "C"
{
   /** Opens an AnaGate I2C device.
      @param pHandle Pointer to an integer, in which the device handle is stored, if device is
                     opened successfully.
      @param nBaudrate Baudrate (two values are supported: 100 and 400)
      @param pcIPAddress Tcp/ip address of the AnaGate device (Port is always 5000).
      @param nTimeout Standard tcp/ip timeout in millseconds.
      @return If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int  I2COpenDevice( int *pHandle, unsigned int nBaudrate,
                                     const char * pcIPAddress, int nTimeout = 1000 );

   /** Closes an open AnaGate I2C device.
       @param hHandle Device handle (from a successfull #OpenDevice call).
       @return If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int  I2CCloseDevice( int hHandle );

   /** Writes data to the I2C partner The data is routed to the partner directly without any transformation,
      so the user of the interface has to build a correct data buffer.
      @param hHandle Device handle (from a successfull #OpenDevice call).
      @param nSlaveAddress Slave address of the I2C partner.
      @param pcBuffer    Pointer to the buffer containing the data to be written.
      @param nBufferLen  Number of bytes to be written.
      @param pnErrorByte Pointer to an integer value which holds the byte number which produces the last error.
                         This is only valid if return value is not equal 0.
      @return. If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int  I2CWrite( int hHandle, unsigned short nSlaveAddress, const char * pcBuffer, int nBufferLen, int * pnErrorByte );

   /** Reads data from the I2C partner. The data offset has to be set by a seperate write command.
      @param hHandle Device handle (from a successfull #OpenDevice call).
      @param nSlaveAddress Slave address of the I2C partner.
      @param pcBuffer    Pointer to the buffer that receives the data read.
      @param nBufferLen  Number of bytes to be read.
      @return. If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int  I2CRead ( int hHandle, unsigned short nSlaveAddress, char * pcBuffer, int nBufferLen );

   /** Reads data from the specified offset to a data buffer.
      @remarks This command is a specialized version of the #Read method and is optimized for I2C EEProms.
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
   ANAGATEDLL_API int  I2CReadEEProm ( int hHandle, unsigned short nChipAddress, unsigned int nOffset,
                                      char * pcBuffer, int nBufferLen, unsigned int nOffsetFormat = 16 );

   /** Writes data to the specified offset. The user of this method has to take care that the offset is at
       page size bound. If the offset of the write command is not at page size bound, all is possible ...
      @remarks This command is a specialized version of the #Write method and is optimized for I2C EEProms.
               <p>Use this method only for EEProm partners</p>
      @param hHandle Device handle (from a successfull #OpenDevice call).
      @param nChipAddress  Chip enable address (bit 1 - 3 of the slave address).
      @param nOffset     Offset position where data is to be written (EEProm partner).
      @param pcBuffer    Pointer to the buffer containing the data to be written.
      @param nBufferLen  Number of bytes to be written.
      @param nOffsetFormat Offset format. To create a I2C command, the method must specify the offset in
                           the correct format. At this time two formats are supported (8-bit and 16bit addresses).
      @return. If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int  I2CWriteEEProm( int hHandle, unsigned short nChipAddress, unsigned int nOffset,
                                       const char * pcBuffer, int nBufferLen, unsigned int nOffsetFormat = 16 );

   /** Retrieves a textual error description of the supplied return code.
      @param nRC Return code.
      @param pcMessage Pointer to a c-style character buffer, in which the retrieved error string is stored.
      @param nMessageLen Length of the supplied charcter buffer. If the error string does not fit into the
             buffer, the string is shortened.
      @return The byte of the returned error string.
   */
   ANAGATEDLL_API int  I2CErrorMessage( int nRC, char *pcMessage, int nMessageLen );

   /** Reads data from digital io register of the partner.
      @param hHandle Device handle (from a successfull #OpenDevice call).
      @param nInputBits   Pointer to the variable that receives the digtial input register read.
      @param nOutputBits  Pointer to the variable that receives the digtial output register read.
      @return. If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int I2CReadDigital( int hHandle, unsigned long * pnInputBits, unsigned long * pnOutputBits );

   /** Writes data to digital io register of the partner.
      @param hHandle Device handle (from a successfull #OpenDevice call).
      @param nOutputBits  Refererence to the variable that hold the digtial IO output register bits to write.
      @return. If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int I2CWriteDigital( int hHandle, unsigned long nOutputBits );

   /** Writes an I2C command sequence to the AnaGate I2C. The sequence already converted to a character buffer.
      @param hHandle Device handle (from a successfull #OpenDevice call).
      @param pcWriteBuffer          Pointer to the buffer containing the data to be written.
      @param nNumberOfBytesToWrite  Number of bytes to be written.
      @param pcReadBuffer           Pointer to the buffer that receives the data read.
      @param nNumberOfBytesToRead   Number of bytes to be read.
      @param nNumberOfBytesRead     Pointer to the variable that receives the number of bytes read.
      @param nByteNumberLastError   If an error occurs then this is set to the byte number which produces the
                                    error. If slave address produces the error it is set to zero.
      @return. If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int I2CSequence ( int          hHandle,
                                    const char * pcWriteBuffer,
                                    int          nNumberOfBytesToWrite,
                                    char *       pcReadBuffer,
                                    int          nNumberOfBytesToRead,
                                    int *        pnNumberOfBytesRead,
                                    int *        pnByteNumberLastError );

   /** Resets the connected I2C bus the AnaGate device.
      @param hHandle Device handle (from a successfull #OpenDevice call).
      @return. If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int  I2CReset( int hHandle );
}

#endif // ANAGATE_DLL_I2C_H
