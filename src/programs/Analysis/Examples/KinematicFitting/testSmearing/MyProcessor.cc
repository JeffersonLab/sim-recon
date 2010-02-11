// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.cc
//

#include <iostream>
using namespace std;

#include <stdio.h>
#include <unistd.h>

#include "MyProcessor.h"
#include "PID/DPhoton.h"
#include "PID/DKinFit.h"
#include "TRACKING/DMCThrown.h"
#include "TRACKING/DTrack.h"


int PAUSE_BETWEEN_EVENTS = 1;
int SKIP_BORING_EVENTS = 0;
int PRINT_ALL=0;
float MAX_EVENTS = 1e20;
int COUNT = 0;
//TFile *fout = NULL;
//TH1F *hpi0[4];

vector<string> toprint;

#define ansi_escape		((char)0x1b)
#define ansi_bold 		ansi_escape<<"[1m"
#define ansi_black		ansi_escape<<"[30m"
#define ansi_red			ansi_escape<<"[31m"
#define ansi_green		ansi_escape<<"[32m"
#define ansi_blue			ansi_escape<<"[34m"
#define ansi_normal		ansi_escape<<"[0m"
#define ansi_up(A)		ansi_escape<<"["<<(A)<<"A"
#define ansi_down(A)		ansi_escape<<"["<<(A)<<"B"
#define ansi_forward(A)	ansi_escape<<"["<<(A)<<"C"
#define ansi_back(A)		ansi_escape<<"["<<(A)<<"D"

//------------------------------------------------------------------
// brun
//------------------------------------------------------------------
jerror_t MyProcessor::init(void)
{
  char name[256];
  ///< Called once at program start.
  fout = new TFile(OUTNAME, "RECREATE","Output file");
  cout << "Opened ROOT file " << OUTNAME <<endl;
  for(int i=0;i<2;i++)
  {
    for(int j=0;j<4;j++)
    {
      for(int k=0;k<5;k++)
      {
        sprintf(name,"hmm%d_%d_%d",i,j,k);
        hmm[i][j][k] = new TH1F(name,name,100,-0.01,0.01);

        sprintf(name,"hcl%d_%d_%d",i,j,k);
        hcl[i][j][k] = new TH1F(name,name,100,0.0,1.0);

        for(int p=0;p<4;p++)
        {
          sprintf(name,"hinvmass%d_%d_%d_%d",i,j,k,p);
          hinvmass[i][j][k][p] = new TH1F(name,name,50,0.0,1.0);
        }

        for(int p=0;p<24;p++)
        {
          sprintf(name,"hpulls%d_%d_%d_%d",i,j,k,p);
          hpulls[i][j][k][p] = new TH1F(name,name,50,-5.0,5.0);
        }
      }
    }
  }

  return NOERROR;
}

//------------------------------------------------------------------
// brun
//------------------------------------------------------------------
jerror_t MyProcessor::brun(JEventLoop *eventLoop, int runnumber)
{
  //vector<string> factory_names = eventLoop->GetFactoryNames();

  usleep(100000); //this just gives the Main thread a chance to finish printing the "Launching threads" message
  cout<<endl;

  // If int PRINT_ALL is set then add EVERYTHING.
  if(PRINT_ALL){
    //toprint = factory_names;
    SKIP_BORING_EVENTS = 0; // with PRINT_ALL, nothing is boring!
  }else{
    // make sure factories exist for all requested data types
    // If a factory isn't found, but one with a "D" prefixed
    // is, go ahead and correct the name.
    vector<string> really_toprint;
    for(unsigned int i=0; i<toprint.size();i++)
    {
      int found = 0;
      int dfound = 0;
      if(found)
        really_toprint.push_back(toprint[i]);
      else if(dfound)
        really_toprint.push_back("D" + toprint[i]);
      else
        cout<<ansi_red<<"WARNING:"<<ansi_normal
          <<" Couldn't find factory for \""
          <<ansi_bold<<toprint[i]<<ansi_normal
          <<"\"!"<<endl;
    }

    toprint = really_toprint;
  }

  // At this point, toprint should contain a list of all factories
  // in dataClassName:tag format, that both exist and were requested.
  // Seperate the tag from the name and fill the fac_info vector.
  fac_info.clear();
  for(unsigned int i=0;i<toprint.size();i++){
    string name = toprint[i];
    string tag = "";
    unsigned int pos = name.rfind(":",name.size()-1);
    if(pos != (unsigned int)string::npos){
      tag = name.substr(pos+1,name.size());
      name.erase(pos);
    }
    factory_info_t f;
    f.dataClassName = name;
    f.tag = tag;
    fac_info.push_back(f);
  }

  cout<<endl;

  return NOERROR;
}

//------------------------------------------------------------------
// evnt
//------------------------------------------------------------------
jerror_t MyProcessor::evnt(JEventLoop *eventLoop, int eventnumber)
{

  if(VERBOSE) cerr << "First thing in evnt....." << endl;

  vector<const DMCThrown*> mcthrown;
  vector<const DTrack*> tracks;

  if(eventLoop->Get(mcthrown))
  { 
    for(int i=0;i<(int)mcthrown.size();i++)
    {
      if(VERBOSE)
      {
        cerr << i << "\t" << mcthrown[i]->momentum().Mag() << "\t";
        cerr << mcthrown[i]->momentum().Theta()*180.0/3.14159 << " ";
        cerr << mcthrown[i]->momentum().Phi() << endl;
      }
    }

  }

  if(eventLoop->Get(tracks, "THROWN"))
  { 
    for(int i=0;i<(int)tracks.size();i++)
    {
      if(VERBOSE)
      {
        cerr << i << "\t" << tracks[i]->momentum().Mag() << "\t";
        cerr << tracks[i]->momentum().Theta()*180.0/3.14159 << " ";
        cerr << tracks[i]->momentum().Phi() << endl;
      }
    }

  }

  COUNT++;
  if(COUNT>=MAX_EVENTS) 
  {
    fini();
    exit(1);
  }
  return NOERROR;
}

//------------------------------------------------------------------
// fini
//------------------------------------------------------------------
jerror_t MyProcessor::fini(void)
{
  fout->cd();
  fout->Write();
  fout->cd();
  fout->Close();
  //  delete fout;
  cout<<endl<<"Closed ROOT file"<<endl;

  ///< Called after last event of last event source has been processed.
  return NOERROR;
}				
