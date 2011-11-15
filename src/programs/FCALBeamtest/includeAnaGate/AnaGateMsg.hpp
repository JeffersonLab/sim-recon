// ----------------------------------------------------------------------------
// Projekt       : Analytica AnaGate API
// File          : AnaGateMsg.hpp
// Author        : Axel Schmidt
// Copyright     : (C) 2004 by Analytica GmbH
// ----------------------------------------------------------------------------
// $Id: AnaGateMsg.hpp,v 1.1.1.1 2010/07/16 18:12:17 ctarbert Exp $
//
// $Log: AnaGateMsg.hpp,v $
// Revision 1.1.1.1  2010/07/16 18:12:17  ctarbert
//
//
// Revision 1.4  2009/04/30 13:24:56  stefanwelisch
// Referenzen als Pointer implementiert
//
// Revision 1.3  2005/06/14 14:46:27  axelschmidt
// new static FromChar members
//
// Revision 1.2  2005/06/07 16:37:38  axelschmidt
// new func to check if message is a data confirmation
//
// Revision 1.1  2004/07/01 16:28:58  AxelSchmidt
// initial version
//
// ----------------------------------------------------------------------------

#ifndef _ANAGATE_MSG_HPP
#define _ANAGATE_MSG_HPP

// Includes -------------------------------------------------------------------
#include <AnaCommon.h>

// packing mode 1 Byte
//#pragma pack(1)

#ifndef USE_DEFAULT_PACKING
   // packing mode 1 Byte
#pragma pack(push,1)
#endif


/** @short Manages the tcp/ip messages in the AnaGate telegram format.

    The tcp/ip messages for all devices of AnaGate series have the following
    data format (byte-aligned):

     <code><pre>
    // 0x00: 2 Byte : length of message (little endian))
    // 0x02: 2 Byte : command code
    // 0x04: 2 Byte : command id
    // 0x06: x Byte : data (depending on command code)
    //       1 Byte : CRC
     </pre></code>

     @author Axel Schmidt
*/
class CAnaGateMsg
{
public:
   /** Standard Constructor */
   CAnaGateMsg( void );

   /** Constructor
      @param nCmdCode command code
      @param nMsgID   message id
      @param pcData   data bytes
      @param nDataLen number of data bytes
   */
   CAnaGateMsg( unsigned short nCmdCode, short nMsgID, const char * pcData, short nDataLen );

   /** Constructor (normally used when parsing incoming telegram)
      @param pcMsgData   data bytes of the complete message.
      @param nMsgDataLen number of message data bytes.
   */
   CAnaGateMsg( const char * pcMsgData, short nMsgDataLen );

   /** Destructor */
   virtual ~CAnaGateMsg ( void );

   /** Get command code.
      @return command code
   */
   unsigned short GetCmdCode( void ) const;

   /** Checks, if message is a confirmation of a data request.
       @return true, if confirmation.
    */
   bool IsConfirmation( void ) const;

   /** Get message id.
      @return message id
   */
   short GetMsgID( void ) const;

   /** Get data bytes.
      @return data bytes
   */
   const std::string & GetData( void ) const;

   /** Get complete message data bytes.
      @return data bytes
   */
   const std::string GetMessageData( void ) const;

   /** Validates the format of the message object (do not validate the CRC).
      @return true, if message format is okay.
   */
   virtual bool IsValid( void ) const;

   /** Validates the message CRC.
      @return true, if the CRC is correct.
   */
   bool IsCRCValid( void ) const;

   /** Get the current crc byte.
      @return crc
   */
   char GetCRC( void ) const;

   /** Computes the CRC for the current message.
      @return CRC
   */
   virtual char ComputeCRC( void ) const;

   /** Converts a short value to its little endian aequivalent.
      @return always an string with 2 bytes.
      */
   static const std::string ToChar( short nDataByte );

   /** Converts a unsigned short value to its little endian aequivalent.
      @return always an string with 2 bytes.
      */
   static const std::string ToChar( unsigned short nDataByte );

   /** Converts a long value to its little endian aequivalent.
      @return always an string with 4 bytes.
      */
   static const std::string ToChar( long nDataByte );

   /** Converts an character array (2 bytes) to an short integer.
      @param nNumber result of the function
      @param pcData  input for the conversion
      */
   static void FromChar( short & nNumber, const char *pcData );

   /** Converts an character array (4 bytes) to an long integer.
      @param nNumber result of the function
      @param pcData  input for the conversion
      */
   static void FromChar( long & nNumber, const char *pcData );
   static void FromChar( unsigned long & nNumber, const char *pcData );

   /** Converts an character array (2 bytes) to an unsigned short integer.
      @param nNumber result of the function
      @param pcData  input for the conversion
      */
   static void FromChar( unsigned short & nNumber, const char *pcData );

   /** Initialize the message.
      @param nCmdCode command code
      @param nMsgID   message id
      @param pcData   data bytes
      @param nDataLen number of data bytes
   */
   void SetMessage( unsigned short nCmdCode, short nMsgID, const char * pcData, short nDataLen );

   /** Initializes the message (with raw incoming message).
      @param pcData   message data bytes
      @param nDataLen number of message data bytes
      @return true, if message data is valid.
   */
   bool SetMessage( const char * pcMsgData, short nMsgDataLen );

protected:
   short m_nDataLen;             ///< length of message
   unsigned short m_nCmdCode;    ///< command code
   short m_nMsgID;               ///< message id
   char  m_nCRC;                 ///< crc
   std::string m_acData;         ///< data bytes
};

// packing mode default
//#pragma pack()

#ifndef USE_DEFAULT_PACKING
   // packing mode default
#pragma pack(pop)
#endif

#endif // #ifndef _ANAGATE_MSG_HPP










