#include <iostream>
#include <iomanip>
#include <sstream>
#include <fcntl.h>
#include <assert.h>
#include "gluex_newCanlib.h"
#include <sys/time.h>

#ifdef ANAGATE
#include <AnaGateDllCan.h>
#else
#include <libpcan.h>
#include <pcan.h>
#endif


//---------------------------------------------------
// Global variables needed to handle messages from AnaGate
// device.
unsigned char wAddr = 0;
unsigned char rAddr = 0;

int m_MessageBuffer[266][11]; //0 is ID, 1 is Message Length, 2-9 is 8(max) integers of data.
//---------------------------------------------------

void WriteMessage( Message& mess, Device& dev ){
  /*
    Send a message over the bus.
    mess = message to send
    dev  = handle for device
  */

#ifdef ANAGATE

  int nFlag = 0x01;  // bit 0 => extended ID, bit 1 => remote telegram
  int id = mess.getID();
  const int nData = (const int)mess.getNData();
  int *data = mess.getData();
  //char dataChar[nData];
  char *dataChar = (char*)malloc(nData*sizeof(char));
  for(int i=0; i<nData; i++ ){
    dataChar[i] = (char)data[i];
  }

  CANWrite( dev.getHandle(), id, dataChar, nData, nFlag );
  
  free(dataChar);

#else 

  int *data = mess.getData();
  int nData = mess.getNData();
  int id = mess.getID();
  TPCANMsg LMsg;      
  LMsg.ID = (u32_t)id;
  LMsg.LEN = (u8_t)nData;
  LMsg.MSGTYPE = MSGTYPE_EXTENDED;
  for( int i=0; i<nData; i++ ){
    LMsg.DATA[i] = data[i];
  }

  CAN_Write( dev.getHandle(), &LMsg );  

#endif

}

//---------------------------------------------------

Message* ReadMessage( Device& dev ){ 
  /*
    Reads a CAN message from device dev and returns
    an Message object.
  */

  Message* recMess = new Message( 0x0, NULL, 0 ); 

#ifdef ANAGATE

  // Set up a function to read CAN message.  When a CAN message 
  // is received, ReadCANMessage will be called...
  int nRC = CANSetCallback( dev.getHandle(), ReadCANMessage ); 
  int data[8];
  int toFlag = 0; //timeoutflag
  
  struct timeval start, current;

  gettimeofday(&start, NULL);
  while(wAddr==rAddr)
    {
      gettimeofday(&current, NULL);
      if((current.tv_sec - start.tv_sec) > 1)
      {
	  cout<<"TIMEOUT";
	  toFlag = 1;
	  wAddr++;
	  break;
       }
    }
  
  for( int i=0; i< m_MessageBuffer[rAddr][1]; i++ ){
      data[i] = (int)(m_MessageBuffer[rAddr][i+2] & 0x000000FF );
      if (toFlag==1)
	{
	data[i]=0;
	m_MessageBuffer[rAddr][0] = 0x00000000;
	}
    }
  
  recMess->setID( m_MessageBuffer[rAddr][0] );  
  recMess->setNData( m_MessageBuffer[rAddr][1] );
  recMess->setData( &data[0] );
  
  // cout<< m_MessageBuffer[rAddr][0] <<" | "<<m_MessageBuffer[rAddr][1]<<" | "<<data[0]<<" | "<<data[1]<<" | "<<data[2]<<" | r: "<<(int)rAddr<<" | w: "<<(int)wAddr<<"\n";
 
  rAddr++;

#else

    TPCANMsg readMsg;
    if (!(CAN_Read( dev.getHandle(), &readMsg ) & CAN_ERR_QRCVEMPTY)){
      
      int data[8];
      for ( int i=0; i<(int)readMsg.LEN; i++ ){
	data[i] = (int)( readMsg.DATA[i] & 0x000000FF );
      }
      recMess->setID( (int)readMsg.ID );
      recMess->setNData( (int)( readMsg.LEN & 0x000000FF ));
      recMess->setData( &data[0] );
      
    }
        
#endif
    return recMess;

}

//-----------------------------------------------
#ifdef ANAGATE

void WINAPI ReadCANMessage(int nIdentifier, const char * pcBuffer, 
			   int nBufferLen, int nFlags, int hHandle )
{
  /* 
     Call back function needed to read messages from the AnaGate
     device.  When AnaGate device receives a message, this function
     should be called.
  */

  for( int i=0; i<nBufferLen; i++ )
    m_MessageBuffer[wAddr][i+2] = 0;

  for( int i = 0; i < nBufferLen; i++ ){
    m_MessageBuffer[wAddr][i+2] = pcBuffer[i];
  }
  
  m_MessageBuffer[wAddr][1] = nBufferLen;
  m_MessageBuffer[wAddr][0] = nIdentifier;
  
  wAddr++;

  return;
}

#endif

//----------------------------------------------------------

Device* OpenCANBus( char* port, char* address ){
  /*
    Function to connect to CAN bus via either PCAN-USB dongle
    or AnaGate TCP-IP gateway, and configure gateway device
    with appropriate baudarate etc.
    Function arguments are only used if opening AnaGate 
    device, in which case:

    char *port = "A", "B", "C" or "D"
    char *address = IP address of gateway
  */


#ifdef ANAGATE

  Device *dev;

  int nPort = 0;  // Anagate port
  if( strcmp( port,"A" ) == 0 )nPort = 0;
  if( strcmp( port,"B" ) == 0 )nPort = 1;
  if( strcmp( port,"C" ) == 0 )nPort = 2;
  if( strcmp( port,"D" ) == 0 )nPort = 3;

  cout << "Trying to open CAN/TCP-IP gateway port " << port << "... " << endl;
  int h = 0; 
  int nRC = 0;
  nRC = CANOpenDevice( &h, TRUE, TRUE, nPort, address, 1000 );
  
  if ( nRC != 0 ){
    cout << "Couldn't open device..... try again!" << endl;
    nRC = CANCloseDevice( h );
    return NULL;
  }
  
  cout << "Gateway is open... " << endl;

  // Set up CAN communication
  // baudrate = 50000, standard modus, No Termination, no Highspeed, no TimeStamp
  nRC = CANSetGlobals( h, 50000, 0, FALSE, FALSE, FALSE );

  dev = new Device( h );
  return dev;

#else

  HANDLE h;
  unsigned int wBTR0BTR1 = 0x472F;  // bit rate 50kB/s
  char *szDevNode = "/dev/pcan32";  // device mount point
  h = LINUX_CAN_Open( szDevNode, O_RDWR );  
  if( h == NULL ){
    cout << "Couldn't open USB device..." << endl;
    cout << "\tFirst check the dongle is connected, then check " ;
    cout << "the device permissions." << endl;
    return NULL;
  }
    
  // A cout line is needed after opening the device
  // otherwise program hangs.
  cout << "Opened USB-Linux device... " << endl;

  CAN_Status( h );
  // initialise comm... bit rate 50kB/s, extended message id
  int error_no = CAN_Init( h, wBTR0BTR1, CAN_INIT_TYPE_EX );
  if( error_no ) {  
    u32_t status = CAN_Status( h );
    cerr << "Status: " << hex << status;
    return NULL;
  }


  Device *dev = new Device( h );
  return dev;

#endif

}

//----------------------------------------------------------

void LEDSwitch(Device &dev,
	       int id, int color)
{
  /*
    Send a CAN message to change the LED.
    id    = base ID to address
    color = 0x00 => red
          = 0x01 => green
	  = 0x02 => off
  */

  // set up data field:
  int data[] = { RLED_INIT, color };

  Message mess( id, &data[0], 2 );
  WriteMessage( mess, dev );

}

//----------------------------------------------------------

void getID(Device &dev, int id, int nReplies){
  /*
    Sends a CAN message to all bases on the bus
    asking them to send back their ID.
    
    id     = base ID to address
    nBases = number of bases on the bus i.e. the 
    number of replies that the program
    will wait for.
  */

  int read1 = 0;

  // set up data field:
  int data[] = { ID_INIT };
  int ndata = 1;

  Message mess( id, &data[0], ndata );
  WriteMessage( mess, dev );

  while( read1 < nReplies ){

    Message *recMess = ReadMessage( dev );
    cout << "board id 0x" << hex << recMess->getID() << endl;
    cout << dec;
    delete recMess;
    read1++;
  
  
  }
  
  
  return;
  
}

//----------------------------------------------------------------------

void TestPulser(Device &dev, int id, int command)
{
  /*
    Send a message to configure or start the test pulser
    id      = base ID to address
    command = TP_FIRE       start test pulser
	      TP_ENBSYNC    enable external sync of pulser
 	      TP_DISSYNC    disable external sync of pulser

  */

  // set up data field:
  int data[] = { PULSER_INIT, command };

  Message mess( id, &data[0], 2 );
  WriteMessage( mess, dev );
}

//----------------------------------------------------------

void HVControl(Device &dev, int id, int cmd1, int cmd2, bool HVread, int nReplies ){
  /* 
     Sends CAN message which switches the HV on/off by
     enabling/disabling  MVB, MVT, HVB, HVT pins.

     id = base ID to address
     HVread = 0   no message returned by base
            = 1   base will return HV status
     nReplies = no of replies to wait for

     1st data byte enables/disables MVB, MVT, HVB, HVT pins

     2nd data byte enables/disables JAM pins

     To enable HV:
     cmd1 = 0xFF: enable HV
     cmd2 = 0x00

     To disable HV:
     cmd1 = 0xF0: disable HV
     cmd2 = 0x00
     
     Enabling of JAM pins should only be done by experts!

  */

  // set up data field:
  int data[] = { HV_INIT, cmd1, cmd2 };

  Message mess( id, &data[0], 3 );
  WriteMessage( mess, dev );

  if( cmd1 == 0xFF && cmd2 == 0x00 )
    cout << "Switching HV on by enabling MVB, MVT, HVB, HVT" << endl;
  if( cmd1 == 0xF0 && cmd2 == 0x00 )
    cout << "Switching HV off by disabling MVB, MVT, HVB, HVT" << endl;

  int read1 = 0;
  if( HVread == true ){ 
 
    while( read1 < nReplies ){
      
      Message *recMess = ReadMessage( dev );
      cout << hex;
      cout << recMess->getID() << endl;
      cout << endl << "Board id: 0x" << hex << recMess->getID() << endl;
      cout << dec;
      int *recData = recMess->getData();
      int HV_stat;
      HV_stat = recData[1];
      cout << "HV status: " << HV_stat << endl;
      if((HV_stat & 0x40)==0x40) cout<<"ENBMVB = ON\n";  // ENBMVB STATUS
      else cout<<"ENBMVB = OFF\n";
    
      if((HV_stat & 0x80)==0x80) cout<<"ENBMVT = ON\n";  // ENBMVT STATUS
      else cout<<"ENBMVT = OFF\n";
    
      if((HV_stat & 0x10)==0x10) cout<<"ENBHVB = ON\n";	// ENBHVB STATS
      else cout<<"ENBHVB = OFF\n";
    
      if((HV_stat & 0x04)==0x04) cout<<"ENBHVT = ON\n";	// ENBHVT STATUS
      else cout<<"ENBHVT = OFF\n";
    
      if((HV_stat & 0x20)==0x20) cout<<"JAMHVB = ON\n";	// JAMHVB STATUS
      else cout<<"JAMHVB = OFF\n";
    
      if((HV_stat & 0x08)==0x08) cout<<"JAMHVT = ON\n";  // JAMHVT STATUS
      else cout<<"JAMHVT = OFF\n";

      delete recMess;
      read1++;
    
    }

  }
}

//----------------------------------------------------------

void SetVoltage(Device &dev, int id, float Voltage){
  /*
    This sends a CAN message instructing the MCU to write
    a new value to the 12-bit DAC.

    id      = base id to address
    Voltage = new voltage to set

    Voltage must be a float between 0 and 2047.5
    12-bit DAC => takes input values 0-4096

    1. Multiply Voltage by 2 and cast as integer
    => integer from 0-4096 i.e. a 12-bit number
    2. Split this value into 8 most significant bits
    and 4 least significant bits, so that they 
    can fit in the 8-bit CAN message data fields.
    3. Send CAN message.
  */

  int DAchannel;
  int value = (int)(Voltage*2.);
  DAchannel = value & 0x0000FFFF;

  // set up data field:
  int data[] = { DA_INIT, ((DAchannel>>4) & 0x0000FF), ((DAchannel<<4) & 0x0000F0) };

  Message mess( id, &data[0], 3 );
  WriteMessage( mess, dev );

  // for debugging... print bytes to screen
  /*
  int dummy1 = (DAchannel>>4) & 0x000000FF;
  int dummy2 = (DAchannel<<4) & 0x000000F0;
  cout << "1st CAN byte: " << DA_INIT << endl;
  cout << "2nd CAN byte: " << dummy1 << endl;
  cout << "3rd CAN byte: " << dummy2 << endl;
  */

}

//----------------------------------------------------------

void PowerDown(Device &dev, int id){
  /*
    Function sends CAN message putting control board
    in low power mode

    id = id of base to address
  */

  // set up data field:
  int data[] = { PD_INIT };

  Message mess( id, &data[0], 1 );
  WriteMessage( mess, dev );


}

//----------------------------------------------------------

void PowerUp(Device &dev, int id){
  /*
    Send a CAN message to MCU - triggers an interrupt in
    firmware which will wake it up from halt mode.
  */

  // set up data field:
  int data[] = { PU_INIT };

  Message mess( id, &data[0], 1 );
  WriteMessage( mess, dev );
 
    
}

//----------------------------------------------------------

void ADC( Device &dev, int id, int command, int nReplies ){
  /*
    Read one of the ADCs.

    id       = id of base to address
    command  = ADC_MVB   medium voltage (bottom)
	       ADC_MVT   medium voltage (top)
	       ADC_DYN   first dynode voltage
	       ADC_CAT   photocathod voltage
	       ADC_DAC   DAC voltage
	       ADC_TEM   temperature monitor
	       ADC_CUR   current monitor
     nReplies = number of replies to wait for
	       
     The message sent back by the base contains the ADC value in 
     channels, this function also calibrates that value into volts
     C or A.
  */



  int i;

  // set up data field:
  int data[] = { ADC_INIT, command };

  Message mess( id, &data[0], 2 );
  WriteMessage( mess, dev );

  
  int value0;
  int value1;
  int read1 = 0;
  while ( read1 < nReplies ){
    Message *recMess = ReadMessage( dev );
    read1++;
    
    int *recData = recMess->getData();

    value0 = recData[1];
    value1 = recData[2];
    int ADCchannels =  (value0 << 8) | value1 ;
    float ADCvoltage = (float)ADCchannels;

    cout << "board id 0x" << hex << recMess->getID() << endl;
    cout << dec;

    if( command == ADC_CAT){
      ADCvoltage = (float)(ADCchannels/1023.) * 3.07 * 1000.;
      cout << ADCvoltage << " V" << endl << endl;
      /*
      cout<<"Photocathode voltage is " << ADCvoltage << " volts.\n";
      cout<<"The photocathode voltage in ADC channels (0 - 1023) is:  " << ADCchannels << endl;
      */
    }
    else if(command == ADC_DYN){
      ADCvoltage = (float)(ADCchannels/1023.) * 3.07 * 1000.;
      cout<<"1st dynode voltage is " << ADCvoltage << " volts.\n";
      cout<<"The 1st dynode voltage in ADC channels (0 - 1023) is:  " << ADCchannels << endl;
    }
    else if(command == ADC_MVT){
      ADCvoltage = (float)(ADCchannels/1023.)*3.07;
      cout<<"Medium voltage (top) is " << ADCvoltage << " volts.\n";
      cout<<"The medium voltage (top) in ADC channels (0 - 1023) is:  " << ADCchannels << endl;
    }
    else if(command == ADC_MVB){
      ADCvoltage = (float)(ADCchannels/1023.)*3.07;
      cout<<"Medium voltage (bottom) is " << ADCvoltage << " volts.\n";
      cout<<"The medium voltage (bottom) in ADC channels (0 - 1023) is:  " << ADCchannels << endl;
    }
    else if(command == ADC_TEM){
      ADCvoltage = (float)(ADCchannels/1023.)*3.07*100;
      cout<<"Temperature is " << ADCvoltage << " C.\n";
      cout<<"The temperature monitor in ADC channels (0 - 1023) is:  " << ADCchannels << endl;
    }
    else if(command == ADC_CUR){
      ADCvoltage = (float)(ADCchannels/1023.)*3.07/200.;
      cout<<"Current is " << ADCvoltage << " A.\n";
      cout<<"The current monitor in ADC channels (0 - 1023) is:  " << ADCchannels << endl;
    }
    else if(command == ADC_DAC){
      ADCvoltage = (float)(ADCchannels/1023.)*3.07;
      cout<<"DAC voltage is " << ADCvoltage << " volts.\n";
      cout<<"The DAC voltage in ADC channels (0 - 1023) is:  " << ADCchannels << endl;
    }
    else{
      cout<<"The voltage in ADC channels (0 - 1023) is:  " << ADCchannels << endl;
    }
    
    delete recMess;
  }

}

//----------------------------------------------------------

void ReadVersion( Device &dev, int id, int nReplies ){

  /*
    Sends a CAN message out with identifier id.

    id = id of base to address; 0x00 to address all bases on bus
    nBases = number of bases on the bus i.e. the number of replies 
    that the program will wait for.
    

  */

  // set up data field:
  int data[] = { VER_INIT };

  Message mess( id, &data[0], 1 );
  WriteMessage( mess, dev );

  int read1 = 0;
  
  while( read1 < nReplies ){

    Message *recMess = ReadMessage( dev );
    int read_id = recMess->getID();
    int *recData = recMess->getData();

    cout << "board id 0x" << hex << read_id;
    cout << dec;
    cout << " is running firmware version " << recData[1];
    cout << "." << recData[2] << endl;
    delete recMess;
    read1++;

  }
  
  return;
  

}

//----------------------------------------------------------

void BootloaderAct( Device &dev, int id ){

  int data[] = { IAP_INIT };
  Message mess( id, &data[0], 1);
  WriteMessage( mess, dev);
  cout<<"IAP_Init sent"<< endl;


  int data2[] = { SYNC };
  Message  mess2( id, &data2[0], 1);
  int data3[] = { 0x00, 0xFF };
  Message  mess3( id, &data3[0], 2);

  for(int i = 0; i < 100; i++)
  {
  usleep(500000);
  cout<<"...Waiting 500ms"<< endl;

  WriteMessage( mess2, dev);
  cout<<"Sending SYNC byte"<< endl;

  usleep(100000);

  WriteMessage( mess3, dev);
  cout<<"Sending GET Byte"<<endl;
  }

}
