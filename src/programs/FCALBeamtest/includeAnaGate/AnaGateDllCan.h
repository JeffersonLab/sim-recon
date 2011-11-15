// ----------------------------------------------------------------------------
// Projekt       : Analytica AnaGate API
// File          : AnaGateDllCAN.h
// Author        : Axel Schmidt
// Copyright     : (C) 2004 by Analytica GmbH
// ----------------------------------------------------------------------------
// $Id: AnaGateDllCan.h,v 1.1.1.1 2010/07/16 18:12:17 ctarbert Exp $
//
// $Log: AnaGateDllCan.h,v $
// Revision 1.1.1.1  2010/07/16 18:12:17  ctarbert
//
//
// Revision 1.17  2009/04/30 13:24:56  stefanwelisch
// Referenzen als Pointer implementiert
//
// Revision 1.16  2009/03/30 11:17:18  stefanwelisch
// Timestamp hinzugefügt
//
// Revision 1.15  2008/05/13 12:19:33  stefanwelisch
// DevicePort bei CANOpenDevice hinzugefügt
//
// Revision 1.14  2008/04/18 10:23:39  stefanwelisch
// neue Parameter bei CAN-Kommunikation (Termination, HighSpeed, ConsiderPrio)
//
// Revision 1.13  2008/01/14 17:00:53  axelschmidt
// zusätzlicher Parameter für CAN callback
//
// Revision 1.12  2007/11/21 09:19:22  axelschmidt
// real c interface for DLL
//
// Revision 1.11  2007/11/08 11:11:43  StefanWelisch
// Defaultwert bei SetGlobals() angegeben
//
// Revision 1.10  2007/11/07 12:09:14  StefanWelisch
// In Get/SetGlobals() die Terminierung mit aufgenommen
//
// Revision 1.9  2007/10/29 08:53:11  axelschmidt
// new define for VB6 support (special DLL)
//
// Revision 1.8  2007/05/21 11:43:39  axelschmidt
// standard timeout 500 -> 1000 ms
//
// Revision 1.7  2005/09/16 09:18:48  axelschmidt
// Digital IO support for I2C and CAN
//
// Revision 1.6  2005/08/15 14:54:49  axelschmidt
// TCP-Protocol changed for global settings
//
// Revision 1.5  2005/08/03 16:28:39  axelschmidt
// SetFiter/GetFilter implementation
//
// Revision 1.4  2005/07/08 12:04:07  AxelSchmidt
// Baurate is int32
//
// Revision 1.3  2005/06/09 16:11:43  axelschmidt
// new CANWrite + CANSetCallback
//
// Revision 1.2  2005/02/09 12:09:08  axelschmidt
// function to set baudrate
//
// Revision 1.1  2005/01/31 16:24:31  axelschmidt
// initial
//
// ----------------------------------------------------------------------------

#ifndef _ANAGATE_DLL_CAN_H
#define _ANAGATE_DLL_CAN_H

// Defines --------------------------------------------------------------------
#include "AnaGateDLL.h"

// Prototyping ----------------------------------------------------------------
#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

   typedef void (WINAPI * CAN_PF_CALLBACK)    ( int nID, const char *pcBuf, int nLen, int nFlags, int hHandle );
   typedef void (WINAPI * CAN_PF_CALLBACK_EX) ( int nID, const char *pcBuf, int nLen, int nFlags, int hHandle, long nSeconds, long nMicroseconds );

   /** Opens an AnaGate CAN device.
      @param pHandle Pointer to an integer, in which the device handle is stored, if device is
                     opened successfully.
      @param bDataCnf     Telegrams are to be confirmed by client and AnaGate CAN.
      @param bMonitor     AnaGate CAN sents all CAN messages to client anyway.
      @param nDevicePort  Portnumber of the Angate (0 - DEVICE_PORT_MAX).
      @param pcIPAddress  Tcp/ip address of the AnaGate device (Port is always 5001).
      @param nTimeout     Standard tcp/ip timeout in millseconds.
      @return If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int  ANAGATEDLL_CALLING CANOpenDevice( int *pHandle, BOOL bDataCnf, BOOL bMonitor, int nDevicePort,
                                     const char * pcIPAddress, int nTimeout = 1000 );

   /** Closes an open AnaGate CAN device.
       @param hHandle Device handle (from a successfull #OpenDevice call).
       @return If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int ANAGATEDLL_CALLING  CANCloseDevice( int hHandle );


   /** Restarts an AnaGate CAN device.
   @param pcIPAddress  Tcp/ip address of the AnaGate device (Port is always 5001).
   @param nTimeout     Standard tcp/ip timeout in millseconds.
   @return If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
    */
   ANAGATEDLL_API int  ANAGATEDLL_CALLING CANRestart( const char * pcIPAddress, int nTimeout = 1000 );


   /** Gets the defined filter set of the current connection.
       @param hHandle Device handle (from a successfull #OpenDevice call).
       @param pnFilter Pointer to 8 software filter definitions. A single filter definiton contains
                        of two 32 bit values, so pnFilter should be a pointer array of 16 integer values.
       @return If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int ANAGATEDLL_CALLING CANGetFilter( int hHandle, int * pnFilter );

   /** Sets the filter set of the current connection.
       @param hHandle Device handle (from a successfull #OpenDevice call).
       @param pnFilter Pointer to 8 software filter definitions. A single filter definiton contains
                        of two 32 bit values, so pnFilter should be a pointer array of 16 integer values.
       @return If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int ANAGATEDLL_CALLING CANSetFilter( int hHandle, const int * pnFilter );

   /** Sets the global settings of an AnaGate CAN device.
       @param hHandle        Device handle (from a successfull #OpenDevice call).
       @param nBaudrate      Baud rate.
       @param nBaudrate      Operating mode.
       @param bTermination   Termination  on/off.
       @param bHighSpeed     HighSpeed    on/off.
       @param bTimestamp     Timestamp in DataIndication on/off.
       @return If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int ANAGATEDLL_CALLING CANSetGlobals( int hHandle, int nBaudrate, unsigned char nOperatingMode, BOOL bTermination, BOOL bHighSpeed, BOOL bTimestampOn );

   /** Gets the global settings of an open AnaGate CAN device.
       @param hHandle        Device handle (from a successfull #OpenDevice call).
       @param nBaudrate      Baud rate.
       @param nBaudrate      Operating mode.
       @param bTermination   Termination  on/off.
       @param bHighSpeed     HighSpeed    on/off.
       @param bTimestamp     Timestamp in DataIndication on/off.
       @return If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int ANAGATEDLL_CALLING CANGetGlobals( int hHandle, int * pnBaudrate, unsigned char * pnOperatingMode, BOOL * pbTermination, BOOL * pbHighSpeed, BOOL * pbTimestampOn);

   /** Sends telegramm to the CAN bus.
       @param hHandle Device handle (from a successfull #OpenDevice call).
       @param nIdentifier            Id of the sendet message.
       @param pcBuffer               Pointer to the buffer containing the CAN data.
       @param nBufferLen             Length of the data buffer.
       @param nFlags                 format flags (bit 0 = extended CAN id, bit 1 = remote telegram).
       @return. If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int ANAGATEDLL_CALLING CANWrite ( int hHandle, int nIdentifier, const char *pcBuffer, int nBufferLen, int nFlags );

   /** Sends telegramm to the CAN bus.
       @param hHandle Device handle (from a successfull #OpenDevice call).
       @param nIdentifier            Id of the sendet message.
       @param pcBuffer               Pointer to the buffer containing the CAN data.
       @param nBufferLen             Length of the data buffer.
       @param nFlags                 format flags (bit 0 = extended CAN id, bit 1 = remote telegram).
       @param pnSeconds              Pointer to long in which the seconds of the timeval will be written.
       @param pnMicroseconds         Pointer to long in which the microseconds of the timeval will be written.
       @return. If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int ANAGATEDLL_CALLING CANWriteEx ( int hHandle, int nIdentifier, const char *pcBuffer, int nBufferLen, int nFlags, long * pnSeconds, long * pnMicroseconds );


   /** Sets the callback function which is called every time an incoming CAN message arrives.
       @param hHandle Device handle (from a successfull #OpenDevice call).
       @param pFunciotn Pointer to the function which is to be called (NULL for reset).
       @return. If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int ANAGATEDLL_CALLING CANSetCallback( int hHandle, CAN_PF_CALLBACK pFunction );

   /** Sets the callbackex function which is called every time an incoming CAN message arrives.
       @param hHandle Device handle (from a successfull #OpenDevice call).
       @param pFunciotn Pointer to the function which is to be called (NULL for reset).
       @return. If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int ANAGATEDLL_CALLING CANSetCallbackEx( int hHandle, CAN_PF_CALLBACK_EX pFunction );

   /** Retrieves a textual error description of the supplied return code.
      @param nRC Return code.
      @param pcMessage Pointer to a c-style character buffer, in which the retrieved error string is stored.
      @param nMessageLen Length of the supplied charcter buffer. If the error string does not fit into the
             buffer, the string is shortened.
      @return The byte of the returned error string.
   */
   ANAGATEDLL_API int ANAGATEDLL_CALLING CANErrorMessage( int nRC, char *pcMessage, int nMessageLen );

   /** Reads data from digital io register of the partner.
      @param hHandle Device handle (from a successfull #OpenDevice call).
      @param nInputBits   Pointer to the variable that receives the digtial input register read.
      @param nOutputBits  Pointer to the variable that receives the digtial output register read.
      @return. If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int ANAGATEDLL_CALLING CANReadDigital( int hHandle, unsigned long * pnInputBits, unsigned long * pnOutputBits );

   /** Writes data to digital io register of the partner.
      @param hHandle Device handle (from a successfull #OpenDevice call).
      @param nOutputBits  Variable that hold the digtial IO output register bits to write.
      @return. If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int ANAGATEDLL_CALLING CANWriteDigital( int hHandle, unsigned long nOutputBits );


   /** Writes the time to the CAN partner.
      @param nSeconds        seconds elapsed since 01.01.1970 (attribut tv_sec of struct timeval).
      @param nMicroseconds   microseconds since the last full second (attribut tv_usec of struct timeval).
      @return. If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int ANAGATEDLL_CALLING CANSetTime( int hHandle, unsigned long nSeconds, unsigned long nMicroseconds );


#ifdef __cplusplus
}
#endif // __cplusplus

#endif // _ANAGATE_DLL_CAN_H
