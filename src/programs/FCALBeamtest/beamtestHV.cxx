/***************************************************************
 *
 * 2011/11/04 Kei Moriya
 *
 * This program is set up to control the HV of each individual
 * base that is connected to the ANAGATE device for the beamtest.
 *
 * We will need to set the HV on each base individually at
 * different voltages. To do this, we need to set up 2
 * CAN ports on the ANAGATE.
 *
 * The program is a modification of newBaseControl in
 * /home/gxfcal/newBaseControl2/newBaseControl.cxx
 * I have tried to keep the arguments the same for the
 * commands that are the same as newBaseControl.
 *
 ***************************************************************/

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <assert.h>
#include <fcntl.h>
#include <unistd.h>
#include <cmath>
#include "gluex_newCanlib.h"

#ifdef ANAGATE
#include <AnaGateDllCan.h>
#else 
#include <libpcan.h>
#endif

using namespace std;

// This function will return the V0 for the trigger module
// from a fit result to a graph of V0 against threshold.
double getV0_trigger(int channel, double threshold);
double getV0_detector(int channel, double threshold);

void ADC_menu( Device &dev , int id, int nReplies );

// To be able to access these variables within the
// functions getV0_detector, we need to put them
// in global scope.
double C0_trigger[10];
double C1_trigger[10];
double C2_trigger[10];
double C0_detector[25];
double C1_detector[25];
double C2_detector[25];

// These are for the power function
// These values were copied from the txt files
// fit_coeffs_trigger_for_wiki_power.txt etc
// The original values came from fits to the PMT tester data,
// where a fit was done to V0 as a function of threshold
// (threshold = 50,75,100,125,150 mV).
// The fit is of the form V0 = A * threshold^(1/B)
const double A_trigger[10] = {1077.2, 1069.88, 1101.5, 1077.88, 1067.42,
			      1072.25, 1091.35, 1022.64, 1070.95, 1058.38};
const double B_trigger[10] = {12.9054, 12.5077, 12.3997, 12.5126, 11.6277,
			      11.8904, 11.6944, 11.7823, 11.865, 11.9711};
const double A_detector[25] = {1013.65, 1029.88, 1004.5, 1070.68, 1013.29, 1064.8, 1091.22, 1048.86, 1087.06, 1064.46,
			       1072.96, 1000, 1044.56, 1041.17, 1061.99, 1052.34, 1017.46, 1020.5, 1079.39, 1039.75,
			       1058.06, 1064.54, 1079.63, 1043.96, 1044.79};
const double B_detector[25] = {10.9025, 11.6602, 10.62, 12.5188, 12.3737, 12.1175, 12.765, 13.0452, 12.7062, 12.2504,
			       12.4227, 12.6281, 12.0193, 13.6443, 11.1453, 12.5427, 12.7099, 12.3328, 11.8016, 11.57,
			       13.2598, 13.2966, 12.781, 12.3665, 12.3983};

//----------------------------------------------------------

int main(int argc, char ** argv){

  char ipaddress[] = "129.79.157.190";
  
  cout << "We assume the IP address of the ANAGATE device is:" << endl
       << ipaddress << endl;
  
  cout << "The program assumes that " << endl
       << "port A on the ANAGATE is the trigger box (10 bases)," << endl
       << "and" << endl
       << "port B on the ANAGATE is the detector (25 bases)" << endl;

  Device *dev_trigger = OpenCANBus("A",ipaddress);
  while(!dev_trigger){
    cout << "was not able to open the CAN device" << endl;
    cout << "Retrying..." << endl;
    sleep(2);
    dev_trigger = OpenCANBus("A",ipaddress);
  }

  Device *dev_detector = OpenCANBus("B",ipaddress);
  while(!dev_detector){
    cout << "was not able to open the CAN device" << endl;
    cout << "Retrying..." << endl;
    sleep(2);
    dev_detector = OpenCANBus("B",ipaddress);
  }

  int Exit_Prog = 1;
  int cmd = 0;

  int nReplies_trigger  = 10;
  int nReplies_detector = 25;

  // If this ext_id is set to 0, then commands will
  // be sent to all bases. To communicate with a
  // specific base, we need to set this to the base
  // we want to talk to.
  int ext_id = 0x00000000;

  cout << "The base IDs are as follows." << endl;
  cout << "(if you are getting timeout errors, then there is a problem with" << endl
       << "the connection)" << endl;

  cout << "--------------------   IDs for trigger   --------------------" << endl;
  getID(*dev_trigger, ext_id, nReplies_trigger);
  cout << "--------------------   IDs for detector  --------------------" << endl;
  getID(*dev_detector, ext_id, nReplies_detector);

  // Insert the base IDs that we want to control.
  // These should be in the order that we have positioned them
  // in the setup.

  // For the triggers, these will be in the order
  // H1, H2, H3, H4, H5,
  // V1, V2, V3, V4, V5,
  const int ext_id_trigger[] = {0x732c79, 0x72c932, 0x72360e, 0x738545, 0x733d10,  // this row all horizontal bars
				0x733e06, 0x72dacc, 0x72bb3f, 0x72579b, 0x738516}; // this row all vertical   bars

  // For the detectors, we start from the bottom row, row 1,
  // and go from left to right in each row, as seen from the back.
  const int ext_id_detector[] = {0x72f75e, 0x72a530, 0x72db49, 0x7256d9, 0x72a5c2,  // row 1 (bottom)
				 0x72361d, 0x7256ae, 0x72355b, 0x7349b1, 0x72bb50,  // row 2
				 0x72c9e2, 0x738553, 0x7281f8, 0x7282ef, 0x72daa4,  // row 3
				 0x72a60a, 0x72c8e8, 0x733d39, 0x72bade, 0x7256f4,  // row 4
				 0x722bc0, 0x7256d8, 0x726dcc, 0x72d915, 0x7256cb}; // row 5 (top)

  // Corresponding HV values for the triggers above, measured with cosmics
  const float fixed_turnon_hv_trigger[] = {1521.96, 1513.45, 1546.51, 1593.05, 1583.08,
					   1551.23, 1577.14, 1465.49, 1542.73, 1536.49};

  // Corresponding HV values for the triggers above, measured with PMTtester2
  // These 
  const float fixed_turnon_hv_detector[] = {1552.22, 1541.39, 1553.12, 1563.52, 1492.01,
					    1570.59, 1583.08, 1516.38, 1581.55, 1567.04,
					    1570.9 , 1461.84, 1547.28, 1487.68, 1610.7 ,
					    1537.7 , 1483.87, 1502.57, 1605.11, 1561.12,
					    1521.46, 1529.19, 1567.64, 1533.96, 1532.8 };


  // To be able to set the thresholds to any arbitrary value, we use
  // the files fit_coeffs_detector_for_wiki.txt and fit_coeffs_trigger_for_wiki.txt
  // These are the fit results to the PMT tester V0 value as a function of threshold.
  // We store the values of the fit coefficients of a 2nd order polynomial in the
  // following arrays.

  /*
  ifstream IN_trigger("fit_coeffs_trigger_for_wiki_poly.txt");
  for(int i=0;i<10;i++)
    IN_trigger >> C0_trigger[i] >> C1_trigger[i] >> C2_trigger[i];

  ifstream IN_detector("fit_coeffs_detector_for_wiki_poly.txt");
  for(int i=0;i<25;i++)
    IN_detector >> C0_detector[i] >> C1_detector[i] >> C2_detector[i];
  */

  int color;
  int command[2];
  int Valid_Response=1;
  float Voltage=0.;

  cout << "New session?" << endl;
  cout << "1. yes (THIS WILL RESET ALL CURRENT HV VALUES TO 0)" << endl;
  cout << "2. no  (do nothing)" << endl;
  int restart;
  cin >> restart;

  if(restart==1){
    // This is command "7", set to 0 V at beginning.
    SetVoltage(*dev_trigger, ext_id, 0);
    SetVoltage(*dev_detector, ext_id, 0);
    
    // This is command "5", turn on at beginning
    command[0]=0xFF;
    command[1]=0x00;
    HVControl(*dev_trigger, ext_id, command[0], command[1], false, nReplies_trigger);
    HVControl(*dev_detector, ext_id, command[0], command[1], false, nReplies_detector);
  }

  while(Exit_Prog == 1){
    cout << "\t 1.    Get Base IDs"<<endl;
    cout << "\t 5.    HV on"<<endl;
    cout << "\t 6.    HV off"<<endl;
    cout << "\t 7.    HV set"<<endl;
    cout << "\t 13.   Read an ADC"<<endl;
    cout << "\t"<<endl;
    cout << "\t 101.  Set single base HV"<<endl;
    cout << "\t 102.  Set single base HV based on PMT thresholds"<<endl;
    cout << "\t 301.  Set voltage individually to fixed values for different bases on TRIGGER"<<endl;
    cout << "\t 302.  Set voltage individually to fixed values for different bases on DETECTOR"<<endl;
    cout << "\t 303.  Set voltage individually to fixed values for different bases on ALL"<<endl;
    cout << "\t 300.  Test what HV will be when setting a specific module threshold" << endl;
    cout << "\t 311.  Set voltage individually based on PMT thresholds to fixed values for different bases on TRIGGER"<<endl;
    cout << "\t 312.  Set voltage individually based on PMT thresholds to fixed values for different bases on DETECTOR"<<endl;
    cout << "\t 313.  Set voltage individually based on PMT thresholds to fixed values for different bases on ALL"<<endl;
    cout << "\t 401.  Set voltage individually to fixed values for different bases on TRIGGER,  inside arrays only"<<endl;
    cout << "\t 402.  Set voltage individually to fixed values for different bases on DETECTOR, inside arrays only"<<endl;
    cout << "\t 403.  Set voltage individually to fixed values for different bases on ALL,      inside arrays only"<<endl;
    cout << "To exit program enter 0"<<endl;
    cout << "> ";
    cin >> cmd;

    if((cmd > -1) && (cmd <19) || cmd==101 || cmd==102 ||
       cmd==300 || cmd==301 || cmd==302 || cmd==303 ||
       cmd==311 || cmd==312 || cmd==313 ||
       cmd==401 || cmd==402 || cmd==403){
      switch(cmd){
      case(0):
	// Quit program
	Exit_Prog=0;
	break;
      case(1):
	// Get ID
	cout << "--------------------   IDs for trigger   --------------------" << endl;
	getID(*dev_trigger, ext_id, nReplies_trigger);
	cout << "--------------------   IDs for detector  --------------------" << endl;
	getID(*dev_detector, ext_id, nReplies_detector);
	break;
      case(5):
	// HV ON
	command[0]=0xFF;
	command[1]=0x00;
	HVControl(*dev_trigger, ext_id, command[0], command[1], false, nReplies_trigger);
	HVControl(*dev_detector, ext_id, command[0], command[1], false, nReplies_detector);
	break;
      case(6):
	// Switch HV OFF
	command[0]=0xF0;
	command[1]=0x00;
	HVControl(*dev_trigger, ext_id, command[0], command[1], false, nReplies_trigger);	  
	HVControl(*dev_detector, ext_id, command[0], command[1], false, nReplies_detector);	  
	break;
      case(7):
	// Set HV
	cout << "Which HV to set?" << endl;
	cout << "1. trigger" << endl;
	cout << "2. detector" << endl;
	int type;
	cin >> type;
	if(!(type== 1 || type ==2)){
	  cout << "Enter either 1 or 2" << endl;
	  cout << "Going back to main menu." << endl;
	  break;
	}
	Valid_Response=1;
	while(Valid_Response==1){
	  cout<<"Enter a value from 0 - 2047.5 (in steps of 0.5):\n";
	  cout<<"> ";
	  cin>>Voltage;
	  if(Voltage>=0. && Voltage<2048.){
	    Valid_Response=0;
	    if(type==1) SetVoltage(*dev_trigger, ext_id, Voltage);
	    if(type==2) SetVoltage(*dev_detector, ext_id, Voltage);
	  } else{
	    cout<<"Invalid Response.\n";
	  }
	}
	break;
      case(13):
	// Read out HV
	{
	  cout << "----------------------------   reading out trigger HV values   ----------------------------" << endl;
	  ADC_menu( *dev_trigger, ext_id, nReplies_trigger );
	  sleep(1);

	  cout << "----------------------------   reading out detector HV values  ----------------------------" << endl;
	  ADC_menu( *dev_detector, ext_id, nReplies_detector );
	  sleep(1);
	}
	break;
      case(101):
	// Set HV on 1 base
	{

	  // First specify which type (trigger/detector)
	  cout << "Which HV to set?" << endl;
	  cout << "1. trigger" << endl;
	  cout << "2. detector" << endl;
	  int type;
	  cin >> type;
	  if(!(type== 1 || type ==2)){
	    cout << "Enter either 1 or 2" << endl;
	    cout << "Going back to main menu." << endl;
	    break;
	  }

	  // The ID of the ith position base is given by my_ext_id[i-1]
	  int position;
	  cout << "Which position base do you want to change?" << endl;
	  cout << "(Type 0 to exit this option without quitting the program)" << endl;
	  cin >> position;
	  if(position==0) break;

	  if(type==1 && position >=10){
	    cout << "Trigger has only 10 bases on." << endl;
	    cout << "Going back to main menu." << endl;
	    break;
	  }
	  if(type==2 && position >=25){
	    cout << "Detector has only 25 bases on." << endl;
	    cout << "Going back to main menu." << endl;
	    break;
	  }

	  // Which HV to set to
	  float hv;
	  cout << "HV value?" << endl;
	  cout << "(Type -999 to exit this option without quitting the program)" << endl;
	  cin >> hv;
	  if(hv==-999) break;
	  if(hv < 0 || 2100 < hv){
	    cout << "Specified HV out of range..." << endl;
	    break;
	  }
	  
	  // Set HV on 1 base
	  if(type==1) SetVoltage(*dev_trigger, ext_id_trigger[position-1], hv);
	  if(type==2) SetVoltage(*dev_detector, ext_id_detector[position-1], hv);

	}
	break;
      case(102):
	// Set HV on 1 base
	{

	  // First specify which type (trigger/detector)
	  cout << "Which HV to set?" << endl;
	  cout << "1. trigger" << endl;
	  cout << "2. detector" << endl;
	  int type;
	  cin >> type;
	  if(!(type== 1 || type ==2)){
	    cout << "Enter either 1 or 2" << endl;
	    cout << "Going back to main menu." << endl;
	    break;
	  }

	  // The ID of the ith position base is given by my_ext_id[i-1]
	  int position;
	  cout << "Which position base do you want to change?" << endl;
	  cout << "(Type 0 to exit this option without quitting the program)" << endl;
	  cin >> position;
	  if(position==0) break;

	  if(type==1 && (position<1 || 10 < position)){
	    cout << "Trigger has only 10 bases on." << endl;
	    cout << "Specify 1-10" << endl;
	    cout << "Going back to main menu." << endl;
	    break;
	  }
	  if(type==2 && (position<1 || 25 < position)){
	    cout << "Detector has only 25 bases on." << endl;
	    cout << "Specify 1-25" << endl;
	    cout << "Going back to main menu." << endl;
	    break;
	  }

	  // Whiat threshold is being set
	  double threshold;
	  cout << "Threshold value?" << endl;
	  cout << "(Type -999 to exit this option without quitting the program)" << endl;
	  cin >> threshold;
	  if(threshold==-999) break;
	  if(threshold < 0 || 1000 < threshold){
	    cout << "Specified HV out of range..." << endl;
	    break;
	  }
	  
	  // Set HV on 1 base
	  if(type==1){
	    cout << "Setting HV on trigger, position " << position << " to be " << getV0_trigger(position,threshold) << " [V]" << endl;
	    SetVoltage(*dev_trigger, ext_id_trigger[position-1], getV0_trigger(position,threshold));
	  }
	  if(type==2){
	    cout << "Setting HV on detector, position " << position << " to be " << getV0_detector(position,threshold) << " [V]" << endl;
	    SetVoltage(*dev_detector, ext_id_detector[position-1], getV0_detector(position,threshold));
	  }

	}
	break;
      case(300):
	// Set to individual values of HV on each base on TRIGGER
	{

	  // First specify which type (trigger/detector)
	  cout << "Which HV to set?" << endl;
	  cout << "1. trigger" << endl;
	  cout << "2. detector" << endl;
	  int type;
	  cin >> type;
	  if(!(type== 1 || type ==2)){
	    cout << "Enter either 1 or 2" << endl;
	    cout << "Going back to main menu." << endl;
	    break;
	  }

	  // The ID of the ith position base is given by my_ext_id[i-1]
	  int position;
	  cout << "Which position base do you want to change?" << endl;
	  cout << "(Type 0 to exit this option without quitting the program)" << endl;
	  cin >> position;
	  if(position==0) break;

	  if(type==1 && position >=10){
	    cout << "Trigger has only 10 bases on." << endl;
	    cout << "Going back to main menu." << endl;
	    break;
	  }
	  if(type==2 && position >=25){
	    cout << "Detector has only 25 bases on." << endl;
	    cout << "Going back to main menu." << endl;
	    break;
	  }
	  
	  // Which threshold to set to
	  cout << "Which threshold to set to?" << endl;
	  double threshold;
	  cin >> threshold;
	  double hv = type==1 ? getV0_trigger(position,threshold) : getV0_detector(position,threshold);
	  if(type==1) cout << "The HV on the trigger, position " << position << " willl be " << hv << " [V]" << endl;
	  if(type==2) cout << "The HV on the detector, position " << position << " willl be " << hv << " [V]" << endl;
	}
	break;
      case(301):
	// Set to individual values of HV on each base on TRIGGER
	{
	  cout << "The values for the trigger will be +100 V compared to" << endl
	       << "the cosmic ray tests." << endl
	       << "Add additional offset of how many volts?" << endl;
	  float offset;
	  cin >> offset;
	  for(int position=0;position<nReplies_trigger;position++){
	    // Set HV on 1 base
	    SetVoltage(*dev_trigger, ext_id_trigger[position], fixed_turnon_hv_trigger[position]+offset);
	  }
	}
	break;
      case(302):
	// Set to individual values of HV on each base on DETECTOR
	{
	  cout << "The values for the trigger will be +100 V compared to" << endl
	       << "the PMT tests." << endl
	       << "Add additional offset of how many volts?" << endl;
	  float offset;
	  cin >> offset;
	  for(int position=0;position<nReplies_detector;position++){
	    // Set HV on 1 base
	    SetVoltage(*dev_detector, ext_id_detector[position], fixed_turnon_hv_detector[position]+offset);
	  }
	}
	break;
      case(303):
	// Set to individual values of HV on each base on BOTH TRIGGER AND DETECTOR
	{
	  cout << "Setting HV for all triggers..." << endl;
	  cout << "The values for the trigger will be +100 V compared to" << endl
	       << "the cosmic ray tests." << endl
	       << "Add additional offset of how many volts?" << endl;
	  float offset;
	  cin >> offset;
	  for(int position=0;position<nReplies_trigger;position++){
	    // Set HV on 1 base
	    SetVoltage(*dev_trigger, ext_id_trigger[position], fixed_turnon_hv_trigger[position]+offset);
	  }
	  cout << "Setting HV for all detectors..." << endl;
	  cout << "The values for the trigger will be +100 V compared to" << endl
	       << "the PMT tests." << endl
	       << "Add additional offset of how many volts?" << endl;
	  cin >> offset;
	  for(int position=0;position<nReplies_detector;position++){
	    // Set HV on 1 base
	    SetVoltage(*dev_detector, ext_id_detector[position], fixed_turnon_hv_detector[position]+offset);
	  }
	}
	break;


      case(311):
	// Set to individual values of HV on each base on TRIGGER, based on fit to V0 vs threshold
	{
	  cout << "What is the common threshold for the trigger?" << endl;
	  double threshold;
	  cin >> threshold;
	  for(int position=1;position<=nReplies_trigger;position++){
	    // Set HV on 1 base
	    cout << "Setting HV on trigger, position " << position << " to be " << getV0_trigger(position,threshold) << " [V]" << endl;
	    SetVoltage(*dev_trigger, ext_id_trigger[position-1], getV0_trigger(position,threshold));
	  }
	}
	break;
      case(312):
	// Set to individual values of HV on each base on DETECTOR, based on fit to V0 vs threshold
	{
	  cout << "What is the common threshold for the detector?" << endl;
	  double threshold;
	  cin >> threshold;
	  cout << "The values for the trigger will be +100 V compared to" << endl
	       << "the PMT tests." << endl
	       << "Add additional offset of how many volts?" << endl;
	  float offset;
	  cin >> offset;
	  for(int position=1;position<=nReplies_detector;position++){
	    // Set HV on 1 base
	    cout << "Setting HV on detector, position " << position << " to be " << getV0_detector(position,threshold) << " [V]" << endl;
	    SetVoltage(*dev_detector, ext_id_detector[position-1], getV0_detector(position,threshold));
	  }
	}
	break;
      case(313):
	// Set to individual values of HV on each base on BOTH TRIGGER AND DETECTOR, based on fit to V0 vs threshold
	{
	  cout << "Setting HV for all triggers..." << endl;
	  cout << "What is the common threshold for the trigger and detector?" << endl;
	  double threshold;
	  cin >> threshold;
	  for(int position=1;position<=nReplies_trigger;position++){
	    // Set HV on 1 base
	    cout << "Setting HV on trigger, position " << position << " to be " << getV0_trigger(position,threshold) << " [V]" << endl;
	    SetVoltage(*dev_trigger, ext_id_trigger[position-1], getV0_trigger(position,threshold));
	  }
	  cout << "Setting HV for all detectors..." << endl;
	  for(int position=1;position<=nReplies_detector;position++){
	    // Set HV on 1 base
	    cout << "Setting HV on detector, position " << position << " to be " << getV0_detector(position,threshold) << " [V]" << endl;
	    SetVoltage(*dev_detector, ext_id_detector[position-1], getV0_detector(position,threshold));
	  }
	}
	break;




      case(401):
	// Set to individual values of HV on each base on TRIGGER, inside arrays only
	{
	  cout << "The values for the trigger will be +100 V compared to" << endl
	       << "the cosmic ray tests." << endl
	       << "Add additional offset of how many volts?" << endl;
	  float offset;
	  cin >> offset;
	  for(int position=0;position<nReplies_trigger;position++){
	    // Set HV on 1 base, skip positions 0,4,5,9
	    if(!(position%5==0 || position%5==4))
	      SetVoltage(*dev_trigger, ext_id_trigger[position], fixed_turnon_hv_trigger[position]+offset);
	  }
	}
	break;
      case(402):
	// Set to individual values of HV on each base on DETECTOR, inside arrays only
	{
	  cout << "The values for the trigger will be +100 V compared to" << endl
	       << "the PMT tests." << endl
	       << "Add additional offset of how many volts?" << endl;
	  float offset;
	  cin >> offset;
	  for(int position=0;position<nReplies_detector;position++){
	    // Set HV on 1 base, skip positions 0,4,5,9
	    if(!(position%5==0 || position%5==4))
	      SetVoltage(*dev_detector, ext_id_detector[position], fixed_turnon_hv_detector[position]+offset);
	  }
	}
	break;
      case(403):
	// Set to individual values of HV on each base on BOTH TRIGGER AND DETECTOR
	{
	  cout << "Setting HV for inside triggers..." << endl;
	  cout << "The values for the trigger will be +100 V compared to" << endl
	       << "the cosmic ray tests." << endl
	       << "Add additional offset of how many volts?" << endl;
	  float offset;
	  cin >> offset;
	  for(int position=0;position<nReplies_trigger;position++){
	    // Set HV on 1 base
	    if(!(position%5==0 || position%5==4))
	      SetVoltage(*dev_trigger, ext_id_trigger[position], fixed_turnon_hv_trigger[position]+offset);
	  }
	  cout << "Setting HV for inside detectors..." << endl;
	  cout << "The values for the trigger will be +100 V compared to" << endl
	       << "the PMT tests." << endl
	       << "Add additional offset of how many volts?" << endl;
	  cin >> offset;
	  for(int position=0;position<nReplies_detector;position++){
	    // Set HV on 1 base
	    if(!(position%5==0 || position%5==4))
	      SetVoltage(*dev_detector, ext_id_detector[position], fixed_turnon_hv_detector[position]+offset);
	  }
	}
	break;
      } // End of switch over cmd
    } // End of condition 0 < cmp < 18
  } // End of Exit_Prog not 1
    
  return 0;
}

//----------------------------------------------------------

//----------------------------------------------------------

// Used by 13
void ADC_menu( Device &dev , int ext_id, int nReplies )
{
	int position;
	int command = 0xFF;
	int value[2];
	float ADCvoltage;
	int ADCchannels;
	cout<<"ADC Menu.\n";
	cout<<"\t 1. Read Medium Voltage Bottom\n";
	cout<<"\t 2. Read Medium Voltage Top\n";
	cout<<"\t 3. Read 1st Dynode\n";
	cout<<"\t 4. Read Photocathode\n";
	cout<<"\t 5. Read DAC voltage\n";
	cout<<"\t 6. Read temperature monitor\n";
	cout<<"\t 7. Read current monitor\n";
	cout<<"> ";
	cin>>position;
	switch(position){
	case(1):
	  command = ADC_MVB;
	  break;
	case(2):
	  command = ADC_MVT;
	  break;
	case(3):
	  command = ADC_DYN;
	  break;
	case(4):
	  command = ADC_CAT;
	  break;
	case(5):
	  command = ADC_DAC;
	  break;
	case(6):
	  command = ADC_TEM;
	  break;
	case(7):
	  command = ADC_CUR;
	  break;
	}
	ADC( dev, ext_id, command, nReplies );
}
//----------------------------------------------------------

double getV0_trigger(int channel, double threshold){
  if(channel<1 || 10<channel){
    cout << "Trigger channel has to be specified between 1-10" << endl;
    cout << "This function will just return 0 for now" << endl;
    return 0;
  }

  // Fit coefficients of V0 against threshold are given as
  // C0_trigger[i] etc for ith trigger.
  // double value = C0_trigger[channel-1] + C1_trigger[channel-1]*threshold + C2_trigger[channel-1]*pow(threshold,2.);

  // Use the power function
  double value = A_trigger[channel-1] * pow(threshold,1./B_trigger[channel-1]);

  // Make sure the value doesn't exceed some crazy value
  if(value < 0 || 1900. < value){
    cout << "The voltage you tried to apply is " << value << " [V], which exceeds" << endl
	 << "the current limit of 1900 [V]. The voltage will be set to 1900 [V]." << endl
	 << "To change this setting, change the function getV0_trigger and recompile." << endl;
    return 1900;
  }
  return value;
}

double getV0_detector(int channel, double threshold){
  if(channel<1 || 25<channel){
    cout << "Detector channel has to be specified between 1-25" << endl;
    cout << "This function will just return 0 for now" << endl;
    return 0;
  }

  // Fit coefficients of V0 against threshold are given as
  // C0_detector[i] etc for ith detector.
  // double value = C0_detector[channel-1] + C1_detector[channel-1]*threshold + C2_detector[channel-1]*pow(threshold,2.);

  // Use the power function
  double value = A_detector[channel-1] * pow(threshold,1./B_detector[channel-1]);

  // Make sure the value doesn't exceed some crazy value
  if(value < 0 || 1900. < value){
    cout << "The voltage you tried to apply is " << value << " [V], which exceeds" << endl
	 << "the current limit of 1900 [V]. The voltage will be set to 1900 [V]." << endl
	 << "To change this setting, change the function getV0_detector and recompile." << endl;
    return 1900;
  }
  return value;
}
