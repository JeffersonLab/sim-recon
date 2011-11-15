#ifndef _CANMSG_
#define _CANMSG_

#ifdef ANAGATE
#include <AnaGateDllCan.h>
#include <linux/types.h>
#else
#include <pcan.h>
#include <libpcan.h>
#endif

using namespace std;

typedef __u32 u32_t;
typedef __u16 u16_t;
typedef __u8 u8_t;

class Message
{

 public:
  Message( int id, int *data, int ndata ):
    m_messageID( id ),
    m_nData( ndata ){
    for( int i=0; i<m_nData; i++ )m_messageData[i] = data[i];
  };

  void setID( int id ){ m_messageID = id; };
  void setNData( int ndata ){ m_nData = ndata; };
  void setData( int* data ){ 
    for( int i=0; i<m_nData; i++ )m_messageData[i] = data[i];
  };
  int getID(){ return m_messageID; };
  int* getData(){ return &m_messageData[0]; };
  int getNData(){ return m_nData; };

 private:
  int m_messageID;
  int m_messageData[8];
  int m_nData;

};

#ifdef ANAGATE

class Device
{
 public:
  Device( int handle ):
    m_handle( handle ){};
    
  void setHandle( int handle ){ m_handle = handle; };
  int getHandle(){ return m_handle; };

 private:
  int m_handle;

};

#else 

class Device
{
 public:
  Device( HANDLE handle ):
    m_handle( handle ){};
    
  void setHandle( HANDLE handle ){ m_handle = handle; };
  HANDLE getHandle(){ return m_handle; };

 private:
  HANDLE m_handle;

};

#endif

// Global variables for control protocol
const int PD_INIT = 0xAA;
const int PU_INIT = 0xAB;
const int ADC_INIT = 0xBB;
const int DA_INIT = 0xCC;
const int IAP_INIT = 0xDD;
const int PULSER_INIT = 0xEE;
const int RC_INIT = 0xFF;
const int HV_INIT = 0x99;
const int ID_INIT = 0x88;
const int RLED_INIT = 0x77;
const int VER_INIT = 0x66;

// LED colours
const int LEDGREEN = 0x01;
const int LEDRED = 0x00;
const int LEDOFF = 0x02;

// ADC control variables 
const int ADC_MVB = 0x00;            // Medium Voltage (bottom)
const int ADC_MVT = 0x01;            // Medium Voltage (top)
const int ADC_DYN = 0x02;            // 1st Dynode
const int ADC_CAT = 0x03;            // Photocathode
const int ADC_DAC = 0x04;            // DAC
const int ADC_TEM = 0x05;            // Temperature Monitor
const int ADC_CUR = 0x06;            // Current Monitor

// Test pulser control variables
const int TP_FIRE = 0x00;          // fire test pulser
const int TP_ENBSYNC = 0x01;       // enable external sync of test pulser
const int TP_DISSYNC = 0x02;       // disable external sync of test pulser

//Bootloader Sync
const int SYNC = 0x7F;


void WriteMessage( Message &mess, Device &dev );
Message* ReadMessage( Device &dev ); 
Device* OpenCANBus( char* port, char* address );
#ifdef ANAGATE
void WINAPI ReadCANMessage( int, const char*, int, int, int );
#endif

void LEDSwitch( Device &dev, int id, int color );
void TestPulser( Device &dev, int id, int command );
void getID( Device &dev, int id, int nBases );
void HVControl( Device &dev, int id, int b1 , int b2, bool HVread, int nReplies );
void SetVoltage( Device &dev, int id, float Voltage );
void PowerDown( Device &dev, int id );
void PowerUp( Device &dev, int id );
void ReadHV( Device &dev, int id, int nReplies );
void ADC( Device &dev, int id, int command, int nReplies );
void ReadVersion( Device &dev, int id, int nReplies );
void BootloaderAct( Device &dev, int id);

#endif
