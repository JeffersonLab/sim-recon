/***************************************************************
 *
 * 2011/11/14 Kei Moriya
 *
 * This program will print out the HV values for each base on
 * the trigger and detector bases.
 *
 * The program uses the same code as beamtestHV.cxx, but
 * this is a ROOT macro that does not need any libraries or
 * compiling.
 * 
 * The base IDs are placed in an array of std::strings, and the
 * turnon curves measured with the PMTtester setup are
 * parametrized by the function
 * HV = A * threshold^(1/B)
 * with A and B parameters. These parameters are also
 * placed in an array of floats.
 *
 * Usage:
 * > root
 * > .x printHVvalues.cc
 *
 ***************************************************************/

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

using namespace std;


// To be able to access these variables within the
// functions getV0_detector, we need to put them
// in global scope.

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

double getV0_trigger(int channel, double threshold);
double getV0_detector(int channel, double threshold);
//----------------------------------------------------------

const int colors[] = {kRed, kGreen, kBlue, kMagenta, kBlack};

int printHVvalues(){

  const int nReplies_trigger  = 10;
  const int nReplies_detector = 25;

  // Insert the base IDs that we want to control.
  // These should be in the order that we have positioned them
  // in the setup.

  // For the triggers, these will be in the order
  // H1, H2, H3, H4, H5,
  // V1, V2, V3, V4, V5
  const int ext_id_trigger[] = {0x732c79, 0x72c932, 0x72360e, 0x738545, 0x733d10,  // this row all horizontal bars
				0x733e06, 0x72dacc, 0x72bb3f, 0x72579b, 0x738516}; // this row all vertical   bars

  // For the detectors, we start from the bottom row, row 1,
  // and go from left to right in each row, as seen from the back.
  const int ext_id_detector[] = {0x72f75e, 0x72a530, 0x72db49, 0x7256d9, 0x72a5c2,  // row 1 (bottom)
				 0x72361d, 0x7256ae, 0x72355b, 0x7349b1, 0x72bb50,  // row 2
				 0x72c9e2, 0x738553, 0x7281f8, 0x7282ef, 0x72daa4,  // row 3
				 0x72a60a, 0x72c8e8, 0x733d39, 0x72bade, 0x7256f4,  // row 4
				 0x722bc0, 0x7256d8, 0x726dcc, 0x72d915, 0x7256cb}; // row 5 (top)
  //___________________________________________________________________________

  ////////////////////////////////////////
  //   Start of main part of program   ///
  ////////////////////////////////////////
  int input = 0;
  int exit  = 0;

  int createdCanvas = 0;

  while(exit==0){
    cout << "--------------------------------------------------------------------------" << endl;
    cout << "Options:" << endl;
    cout << "   1. Print out base IDs" << endl;
    cout << "   2. Print out HV values based on threshold for trigger  (10 values)" << endl;
    cout << "   3. Print out HV values based on threshold for detector (25 values)" << endl;
    cout << "   4. Create graphs showing HV values as a function of threshold" << endl;
    cout << "-999. Quit" << endl;
    cout << "--------------------------------------------------------------------------" << endl;
    cout << ">";

    cin >> input;

    switch(input){
      case(-999):
	exit = 1;
	break;
    case(1):
      cout << "------------------------- Base IDs for trigger  -------------------------------" << endl;
      cout << "\thorizontal bars:" << endl;
      for(int i=0;i<5;i++)
	cout << "\t\t0x" << hex << ext_id_trigger[i] << endl;
      cout << "\tvertical bars:" << endl;
      for(int i=0;i<5;i++)
	cout << "\t\t0x" << hex << ext_id_trigger[i+5] << endl;

      cout << "------------------------- Base IDs for detector -------------------------------" << endl;
      cout << "\t1st (bottom) row:" << endl;
      for(int i=0;i<5;i++)
	cout << "\t\t0x" << hex << ext_id_detector[i] << endl;
      cout << "\t2nd row:" << endl;
      for(int i=0;i<5;i++)
	cout << "\t\t0x" << hex << ext_id_detector[i+5] << endl;
      cout << "\t3rd row:" << endl;
      for(int i=0;i<5;i++)
	cout << "\t\t0x" << hex << ext_id_detector[i+10] << endl;
      cout << "\t4th row:" << endl;
      for(int i=0;i<5;i++)
	cout << "\t\t0x" << hex << ext_id_detector[i+15] << endl;
      cout << "\t5th (top) row:" << endl;
      for(int i=0;i<5;i++)
	cout << "\t\t0x" << hex << ext_id_detector[i+20] << endl;
      
      break;
    case(2):
      cout << "Threshold value for triggers [mV]?" << endl;
      double threshold;
      cin >> threshold;
      cout << "Specify an additional offset that will be added to the fitted function" << endl;
      double offset;
      cin >> offset;
      cout << "------------------------- threshold values for trigger  -------------------------------" << endl;
      cout << "\thorizontal bars:" << endl;
      for(int i=0;i<5;i++)
	cout << "\t\t0x" << hex << ext_id_trigger[i] << "\t" << getV0_trigger(i+1,threshold)+offset << endl;
      cout << "\tvertical bars:" << endl;
      for(int i=0;i<5;i++)
	cout << "\t\t0x" << hex << ext_id_trigger[i] << "\t" << getV0_trigger(i+1+5,threshold)+offset << endl;
      break;
    case(3):
      cout << "Threshold value for detector [mV]?" << endl;
      double threshold;
      cin >> threshold;
      cout << "Specify an additional offset that will be added to the fitted function" << endl;
      double offset;
      cin >> offset;
      cout << "------------------------- threshold values for trigger  -------------------------------" << endl;
      cout << "(starting from bottom row)" << endl;
      for(int i=0;i<25;i++)
	cout << "\t\t0x" << hex << ext_id_detector[i] << "\t" << getV0_detector(i+1,threshold)+offset << endl;
      break;
    case(4):
      createdCanvas++;
      char canvasname[50];
      char canvastitle[50];

      /////////////////////////////
      ///  canvas for triggers  ///
      /////////////////////////////
      sprintf(canvasname,"c_trigger_%d",createdCanvas);
      sprintf(canvastitle,"trigger %d",createdCanvas);
      TCanvas *c_trigger = new TCanvas(canvasname,canvastitle);
      TF1 *f_trigger[10];
      TLegend *legend_trigger = new TLegend(0.20,0.30,0.60,0.93);
      legend_trigger->SetNColumns(2);
      legend_trigger->SetHeader("V_{0} for trigger bases");
      for(int i=0;i<10;i++){
	char funcname[50];
	sprintf(funcname,"f_trigger_%d___%d",i,createdCanvas);
	f_trigger[i] = new TF1(funcname,"[0]*pow(x,1./[1])",0,2000);
	f_trigger[i]->SetParameter(0,A_trigger[i]);
	f_trigger[i]->SetParameter(1,B_trigger[i]);
	f_trigger[i]->SetTitle("");
	f_trigger[i]->GetXaxis()->SetTitle("threshold [mV]");
	f_trigger[i]->GetXaxis()->CenterTitle();
	f_trigger[i]->GetYaxis()->SetTitle("V_{0} [V]");
	f_trigger[i]->GetYaxis()->CenterTitle();
	f_trigger[i]->SetLineStyle(i/5);
	f_trigger[i]->SetLineWidth(1);
	f_trigger[i]->SetLineColor(colors[i%5]);
	f_trigger[i]->SetMarkerStyle(20+i/5);
	f_trigger[i]->SetMarkerSize(1.5);
	f_trigger[i]->SetMarkerColor(colors[i%5]);
	f_trigger[i]->SetMinimum(1000);
	f_trigger[i]->SetMaximum(2500);
	if(i==0) f_trigger[i]->Draw();
	else     f_trigger[i]->Draw("same");
	char name[30];
	sprintf(name,"0x%x",ext_id_trigger[i]);
	legend_trigger->AddEntry(f_trigger[i],name,"L");
      }
      legend_trigger->Draw("same");
      c_trigger->SaveAs("V0_vs_threshold_trigger.ps");

      /////////////////////////////
      ///  canvas for detector  ///
      /////////////////////////////
      sprintf(canvasname,"c_detector_%d",createdCanvas);
      sprintf(canvastitle,"detector %d",createdCanvas);
      TCanvas *c_detector = new TCanvas(canvasname,canvastitle);
      TF1 *f_detector[25];
      TLegend *legend_detector = new TLegend(0.20,0.30,0.60,0.93);
      legend_detector->SetNColumns(2);
      legend_detector->SetHeader("V_{0} for detector bases");
      for(int i=0;i<25;i++){
	char funcname[50];
	sprintf(funcname,"f_detector_%d___%d",i,createdCanvas);
	f_detector[i] = new TF1(funcname,"[0]*pow(x,1./[1])",0,2000);
	f_detector[i]->SetParameter(0,A_detector[i]);
	f_detector[i]->SetParameter(1,B_detector[i]);
	f_detector[i]->SetTitle("");
	f_detector[i]->GetXaxis()->SetTitle("threshold [mV]");
	f_detector[i]->GetXaxis()->CenterTitle();
	f_detector[i]->GetYaxis()->SetTitle("V_{0} [V]");
	f_detector[i]->GetYaxis()->CenterTitle();
	f_detector[i]->SetLineStyle(i/5);
	f_detector[i]->SetLineWidth(1);
	f_detector[i]->SetLineColor(colors[i%5]);
	f_detector[i]->SetMarkerStyle(20+i/5);
	f_detector[i]->SetMarkerSize(1.5);
	f_detector[i]->SetMarkerColor(colors[i%5]);
	f_detector[i]->SetMinimum(1000);
	f_detector[i]->SetMaximum(2500);
	if(i==0) f_detector[i]->Draw();
	else     f_detector[i]->Draw("same");
	char name[30];
	sprintf(name,"0x%x",ext_id_detector[i]);
	legend_detector->AddEntry(f_detector[i],name,"L");
      }
      legend_detector->Draw("same");
      c_detector->SaveAs("V0_vs_threshold_detector.ps");

      break;
    default:
      break;
    }
    // end of switch over input
    
  }
  // end of while loop over input
    
  return 0;
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

  // Will not have this restriction on HV for this program.
  // This will be applied in the beamtestHV program.
  /*
  // Make sure the value doesn't exceed some crazy value
  if(value < 0 || 1900. < value){
    cout << "The voltage you tried to apply is " << value << " [V], which exceeds" << endl
	 << "the current limit of 1900 [V]. The voltage will be set to 1900 [V]." << endl
	 << "To change this setting, change the function getV0_trigger and recompile." << endl;
    return 1900;
    }
  */

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

   // Will not have this restriction on HV for this program.
  // This will be applied in the beamtestHV program.
  /*
  // Make sure the value doesn't exceed some crazy value
  if(value < 0 || 1900. < value){
    cout << "The voltage you tried to apply is " << value << " [V], which exceeds" << endl
	 << "the current limit of 1900 [V]. The voltage will be set to 1900 [V]." << endl
	 << "To change this setting, change the function getV0_detector and recompile." << endl;
    return 1900;
  }
  */
  
  return value;
}
