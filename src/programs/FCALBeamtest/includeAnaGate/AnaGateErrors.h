// ----------------------------------------------------------------------------
// Projekt       : Analytica AnaGate API
// File          : AnaGateErrors.h
// Author        : Axel Schmidt
// Copyright     : (C) 2004 by Analytica GmbH
// ----------------------------------------------------------------------------
// $Id: AnaGateErrors.h,v 1.1.1.1 2010/07/16 18:12:17 ctarbert Exp $
//
// $Log: AnaGateErrors.h,v $
// Revision 1.1.1.1  2010/07/16 18:12:17  ctarbert
//
//
// Revision 1.12  2009/03/30 11:17:18  stefanwelisch
// Timestamp hinzugefügt
//
// Revision 1.11  2008/12/12 12:23:42  StefanWelisch
// neuen Errorcode für CAN Setbaudrate hinzugefügt
//
// Revision 1.10  2007/11/12 10:50:20  axelschmidt
// Byte-Alignment
//
// Revision 1.9  2007/11/08 09:07:16  axelschmidt
// additional CAN error removed
//
// Revision 1.8  2007/11/07 12:11:45  StefanWelisch
// Neue Fehlermeldunen für Anagate CAN
//
// Revision 1.7  2007/06/05 14:01:38  axelschmidt
// new function AnaGateDeviceName
//
// Revision 1.6  2007/05/31 10:15:37  stefanwelisch
// Fehlermeldungen für Anagate-Renesas hinzugefügt
//
// Revision 1.5  2005/07/28 16:24:21  axelschmidt
// error codes angepasst
//
// Revision 1.4  2005/07/22 07:11:59  AxelSchmidt
// new global error
//
// Revision 1.3  2004/08/13 12:51:49  AxelSchmidt
// defines I2C-Error changed
//
// Revision 1.2  2004/07/08 16:02:25  AxelSchmidt
// new I2C errors
//
// Revision 1.1  2004/07/01 16:28:58  AxelSchmidt
// initial version
//
// ----------------------------------------------------------------------------

#ifndef _ANAGATE_ERRORS_H
#define _ANAGATE_ERRORS_H

// Includes -------------------------------------------------------------------
#include <string>

// Defines --------------------------------------------------------------------

/** Possible return codes from the AnaGate API methods. */
enum EAnaGateRetcode
{
   ERR_NONE                   =  0,             ///< No error = success.
   ERR_TCPIP_SOCKET           = 0x020000,       ///< TCP/IP socket error.
   ERR_TCPIP_NOTCONNECTED     = 0x030000,       ///< Not connected to TCP/IP partner.
   ERR_TCPIP_TIMEOUT          = 0x040000,       ///< TCP/IP timeout occured.
   ERR_TCPIP_CALLNOTALLOWED   = 0x050000,       ///< TCP/IP request not allowed now, stack is intializing.
   ERR_TCPIP_NOT_INITIALIZED  = 0x060000,       ///< TCP/IP stack is not initialized.

   ERR_INVALID_CRC            = 0x0A0000,       ///< AnaGate TCP/IP telegram has invalid CRC.
   ERR_INVALID_CONF           = 0x0B0000,       ///< AnaGate TCP/IP telegram is not confirmed by partner.
   ERR_INVALID_CONF_DATA      = 0x0C0000,       ///< AnaGate TCP/IP telegram confirmation is invalid.

   ERR_INVALID_DEVICE_HANDLE  = 0x900000,       ///< Device handle is invalid.
   ERR_INVALID_DEVICE_TYPE    = 0x910000,       ///< Device type invalid (handle is connected device type, which does not support called function).
   ERR_INVALID_PARAM          = 0x920000,       ///< Invalid parameters for called function.

   ERR_OPEN_MAX_CONN          = 0x000001,       ///< Open failed, max. connections
   ERR_CONFIG_MEMORY          = 0x000002,       ///< read error eeprom memory
   
   ERR_OP_CMD_FAILED          = 0x0000FF,       ///< command failed, unkown reason

   // device specific error codes 
   ERR_I2C_NACK                   = 0x000120,       ///< I2C NACK
   ERR_I2C_TIMEOUT                = 0x000121,       ///< I2C time out

   ERR_CAN_NACK                   = 0x000220,         ///< CAN NACK   (Obsolete)
   ERR_CAN_TX_ERROR               = 0x000221,         ///< CAN Transmit errror
   ERR_CAN_SEND                   = ERR_CAN_TX_ERROR, ///< CAN Send error (compatibility)  
   ERR_CAN_TX_BUF_OVERLOW         = 0x000222,         ///< CAN buffer overflow
   ERR_CAN_TX_MLOA                = 0x000223,         ///< CAN Lost Arbitration
   ERR_CAN_NO_VALID_BAUDRATE      = 0x000224,         ///< CAN Setting no valid Baudrate

   ERR_RENESAS_TIMEOUT            = 0x000920,         ///< Renesas timeout
   ERR_RENESAS_INVALID_ID         = 0x000921,         ///< Renesas Invalid ID
   ERR_RENESAS_FLASH_ERASE_FAILED = 0x000922,         ///< Renesas failed erase the flash
   ERR_RENESAS_PAGE_PROG_FAILED   = 0x000923,         ///< Renesas failed prog the page
   
   ERR_DISPLAY_WRONG_BITMAP_SIZE_X    = 0x000A20,     ///< Display wrong bitmap size x
   ERR_DISPLAY_WRONG_BITMAP_SIZE_Y    = 0x000A21,     ///< Display wrong bitmap size y
   ERR_DISPLAY_WRONG_BITMAP_SIZE_DATA = 0x000A22      ///< Display wrong bitmap size data
};

/** Returns a description of the given return code.
   @param nRetCode Return code.
   @return Description of return code.
*/
std::string AnaGateErrorMsg( int nRetCode );

/** Returns the name of the AnaGate device.
   @param nDeviceType Device type.
   @return Name of Anagate device.
*/
std::string AnaGateDeviceName( int nDeviceType );

#endif // #ifndef _ANAGATE_ERRORS_H
