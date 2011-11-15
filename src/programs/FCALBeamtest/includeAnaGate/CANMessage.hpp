// ----------------------------------------------------------------------------
// Projekt       : Analytica AnaGate API
// File          : CANMessage.hpp
// Author        : Axel Schmidt
// Copyright     : (C) 2005 by Analytica GmbH
// ----------------------------------------------------------------------------
// $Id: CANMessage.hpp,v 1.1.1.1 2010/07/16 18:12:17 ctarbert Exp $
//
// $Log: CANMessage.hpp,v $
// Revision 1.1.1.1  2010/07/16 18:12:17  ctarbert
//
//
// Revision 1.7  2009/04/30 13:25:03  stefanwelisch
// Referenzen als Pointer implementiert
//
// Revision 1.6  2009/03/30 11:17:18  stefanwelisch
// Timestamp hinzugefügt
//
// Revision 1.5  2007/11/12 10:50:20  axelschmidt
// Byte-Alignment
//
// Revision 1.4  2005/08/19 13:21:34  axelschmidt
// CAN-Message mit Zeitstempel
//
// Revision 1.3  2005/08/03 16:29:47  axelschmidt
// new message format: remote telegram
//
// Revision 1.2  2005/07/22 16:14:46  AxelSchmidt
// baudrate values changed
//
// Revision 1.1  2005/07/11 16:18:43  axelschmidt
// initial
//
// ----------------------------------------------------------------------------

#ifndef _CAN_MESSAGE_HPP
#define _CAN_MESSAGE_HPP

// Includes -------------------------------------------------------------------
#include <sys/types.h>
#include <iostream>
#include <AnaLocalTime.hpp>

// packing mode 1 Byte
//#pragma pack(1)

#ifndef USE_DEFAULT_PACKING
   // packing mode 1 Byte
#pragma pack(push,1)
#endif

/** Class implements a CAN message.
    @author  Axel Schmidt
    @version 1.0
*/

class CCanMessage
{
   public:
      enum { MAX_CAN_LENGTH  = 8 };

      /** Supported baudrates for CAN bus */
      enum EBaudrate
      {
         KB_20    =   20000,  ///< 20kBit/s
         KB_50    =   50000,  ///< 50kBit/s
         KB_100   =  100000,  ///< 100kBit/s
         KB_250   =  250000,  ///< 250kBit/s
         KB_500   =  500000,  ///< 500kBit/s
         KB_1000  = 1000000,  ///< 1MBit/s
         MB_1     = KB_1000   ///< 1MBit/s
      };

      /** Format flags */
      enum EFormatFlags
      {
         FMT_EXTENDED  = 0x01,    ///< CAN identifier is estended format (29 bit)
         FMT_REMOTE    = 0x02,    ///< telegram is a remote telegram (no data)
         FMT_TIME      = 0x04,    ///< telegram includes timestamp
         FMT_TCP_MSG   = 0x08     ///< telegram was sending from another TCP client to anagate
      };

      /** Construktor. Creates a telegram with id 0 and an empty data buffer. */
      CCanMessage();

      /** Constructor.
          @param nIdentifier    Telegram id.
          @param pcData         Pointer to data buffer.
          @param nDataLength    Length og data buffer.
          @param nFormatFlags   Format flags.
          @param bShowTimeDiff  Show time diff by using stream operator.
      */
      CCanMessage( const unsigned int nIdentifier,
                   const unsigned char * pcData,
                   const int nDataLength,
                   const int nFormatFlags   = 0,
                   const bool bIncoming     = true,
                   const bool bShowTimeDiff = false );

      /** Destructor */
      virtual ~CCanMessage();

      /** Sets the CAN telegram buffer. A CAN telegram has a maximum of 8 data bytes.
          @param  pcData      Pointer to message buffer.
          @param  nDataLength Message length.
      */
      void SetCANData(const unsigned char * pcData, int nDataLength);

      /** Sets the time of CAN telegram.
          @param  sTime       struct of timeval.
      */
      void SetTime( timeval sTime );

      /** Calculates the timediff to CAN telegram.
         @param  void.
      */
      void CalcDiffTime( CAnaLocalTime & oTimeLastMessage); //, unsigned int & nMicrosecondsLastMessage );

      /** Sets the timediff to previous CAN telegram.
          @param  void.
      */
      void SetDiffTime( void );

      /** Gets the number of data bytes in the CAN telegram.
         @return Buffer length.
      */
      int GetBufferLength( void ) const;

      /** Gets a pointer to the data buffer of teh CAN telegram.
         @return Identifier.
      */
      const unsigned char * GetBuffer( void ) const;

      /** Gets the identifier of the CAN telegram.
         @return Identifier.
      */
      unsigned int GetIdentifier( void ) const;

      /** Returns true, if the CAN identifier is in extended bit format.
         @return True or false.
      */
      bool  GetIsExtended( void ) const;

      /** Returns true, if the CAN telegram is a remote telegram.
         @return True or false.
      */
      bool  GetIsRemote( void ) const;

      /** Returns true, if the CAN telegram has timestamp from anagate.
         @return True or false.
      */
      bool GetHasTimestamp( void ) const;

      /** Returns true, if the CAN telegram is a message from an other tcp client.
         @return True or false.
      */
      bool GetIsTCPMsg( void ) const;

      /** Returns the format flags of the telegram.
         @return format flags..
      */
      int  GetFormatFlags( void ) const;

      /** Returns true, if the CAN telegram is a incoming message over TCP  / false for generated messages.
         @return True or false.
      */
      bool GetIsIncoming( void ) const;

      /** Gets the time stamp of a received CAN telegram.
         @return Receive time.
      */
      CAnaLocalTime GetReceiveTime( void ) const;

      
      /** Gets the timediff to previous message
         @return timediff.
      */
      CAnaTimeDiff GetTimeDiff( void ) const;

      /** Gets the microseconds from time of CAN telegram.
      @return microseconds.
       */
      //unsigned int GetMicroseconds( void ) const;

      /** Ausgabestreamoperator */
      friend std::ostream & operator<<( std::ostream & out, const CCanMessage & oCanMsg );

      /** Virtuale Darstellungsmethode für Ausgabeströme
          @param  out Referenz auf den Ausgabestream
          @return Referenz auf den Ausgabestream
      */
      virtual std::ostream & Put( std::ostream & out ) const;

protected:
      int m_nDataLength;                                 ///< Message length
      unsigned char  m_acDataBuffer[MAX_CAN_LENGTH];     ///< Data buffer
      unsigned int   m_nIdentifier;                      ///< Message id
      CAnaLocalTime  m_oReceiveTime;                     ///< Time of message receipt

      bool           m_bShowTimeDiff;                    ///< TRUE show diff time / FALSE show absolute time
      CAnaTimeDiff   m_oTimeDiff;                        ///< Timediff to previous message
      
      bool           m_bIncoming;                        ///< TRUE message received from anagate / FALSE message is going to anagate

private:
      int m_nFormatFlags;    ///< Format flags

};

// packing mode default
//#pragma pack()

#ifndef USE_DEFAULT_PACKING
   // packing mode default
#pragma pack(pop)
#endif

//==============================================================================
inline std::ostream & operator<<( std::ostream & out, const CCanMessage & oCanMsg )
{
   return oCanMsg.Put( out );
}

#endif // #ifndef _CAN_MESSAGE_HPP










