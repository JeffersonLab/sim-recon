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
  for(int i=0;i<4;i++)
  {

    sprintf(name,"hclGood%d",i);
    hclGood[i] = new TH1F(name,name,10,-0.5,9.5);

    for(int j=0;j<4;j++)
    {
      for(int k=0;k<5;k++)
      {

        sprintf(name,"hpi0%d_%d_%d",i,j,k);
        hpi0[i][j][k] = new TH1F(name,name,160,0.0,0.8);

        sprintf(name,"hpi0fit%d_%d_%d",i,j,k);
        hpi0fit[i][j][k] = new TH1F(name,name,160,0.0,0.8);

        sprintf(name,"hcl%d_%d_%d",i,j,k);
        hcl[i][j][k] = new TH1F(name,name,100,0.0,1.0);

        for(int p=0;p<6;p++)
        {

          sprintf(name,"hpulls%d_%d_%d_%d",i,j,k,p);
          hpulls[i][j][k][p] = new TH1F(name,name,100,-5.0,5.0);

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
      //for(unsigned int j=0;j<factory_names.size();j++)
      {
        //if(factory_names[j] == toprint[i])found = 1;
        //if(factory_names[j] == "D" + toprint[i])dfound = 1;
      }
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

  vector<const DPhoton*> photons;
  vector<DKinematicData> kindata;
  vector<const DMCThrown*> mcthrown;

  DKinFit *kfit = new DKinFit();
  kfit->SetVerbose(VERBOSE);

  float numclGood = 0;
  float pulls[6];
  float clCut[5] = {-1.0, 0.01, 0.1, 0.2, 0.5};

  if(VERBOSE>0)
  {
    if(eventLoop->Get(mcthrown))
    { 
      for(int i=0;i<(int)mcthrown.size();i++)
      {
        cerr << COUNT << ": type/pdgtype/mass: " << mcthrown[i]->type << " " << mcthrown[i]->pdgtype << " " <<  mcthrown[i]->mass << endl;
        cerr << "E/p/theta/phi: " << mcthrown[i]->E << " "<< mcthrown[i]->p << " "<< mcthrown[i]->theta << " "<< mcthrown[i]->phi << endl;
        cerr << "E/p/theta/phi: " << mcthrown[i]->lorentzMomentum().E() << " " << mcthrown[i]->lorentzMomentum().Rho() << " " << mcthrown[i]->lorentzMomentum().Theta() << " " << mcthrown[i]->lorentzMomentum().Phi() << endl;
        cerr << "photon0 errmat: ";
        cerr << mcthrown[i]->errorMatrix()(0,0) << " ";
        cerr << mcthrown[i]->errorMatrix()(1,1) << " ";
        cerr << mcthrown[i]->errorMatrix()(2,2) << " ";
        cerr << mcthrown[i]->errorMatrix()(3,3)  << endl;
      }
    }
  }
  // If in FCAL tag == 0
  // If in BCAL tag == 1
  if(eventLoop->Get(photons))
  { 
    if(COUNT%100==0) cerr << COUNT << " " << MAX_EVENTS << endl;
    if(COUNT>=MAX_EVENTS) 
    {
      fini();
      exit(1);
    }
    int tagCount = 1;
    int numPhotons = (int) photons.size();
    int numPhotonsCount = numPhotons - 1;
    if(VERBOSE>0)
    {
      cerr << "Before anything..." << endl;
      cerr << "tagCount: " << tagCount << endl;
      cerr << "numPhotons: " << numPhotons << endl;
      cerr << "numPhotonsCount: " << numPhotonsCount << endl;
    }
    if(numPhotons>=2 && numPhotons<=4)
    {
      if(VERBOSE)
      {
        cerr << "photon0 errmat: ";
        cerr << photons[0]->errorMatrix()(0,0) << " ";
        cerr << photons[0]->errorMatrix()(1,1) << " ";
        cerr << photons[0]->errorMatrix()(2,2) << " ";
        cerr << photons[0]->errorMatrix()(3,3)  << endl;
      }
      for(int i=0;i<numPhotons;i++)
      {
        if(VERBOSE>0)
        {
          cerr << COUNT << " photons: " << photons[i]->mass() << " ";
          cerr << photons[i]->energy() << " "<< photons[i]->px() << " "<< photons[i]->py() << " "<< photons[i]->pz() << endl;
        }
        tagCount = 1;
        kindata.clear();
        //kindata.push_back((DKinematicData*)photons[i]);
        kindata.push_back(*photons[i]);
        tagCount += photons[i]->getTag();

        for(int j=i+1;j<numPhotons;j++)
        {
          //kindata.push_back((DKinematicData*)photons[j]);
          kindata.push_back(*photons[j]);
          tagCount += photons[j]->getTag();

          if(kindata.size()==2)
          {
            if(VERBOSE>0)
            {
              cerr << "tagCount: " << tagCount << endl;
              cerr << "numPhotons: " << numPhotons << endl;
              cerr << "numPhotonsCount: " << numPhotonsCount << endl;
            }
            float pi0mass = (kindata[0].lorentzMomentum() + kindata[1].lorentzMomentum() ).M();

            if(fabs(pi0mass-MASS)<0.15)
            {
              hpi0[0][0][0]->Fill( pi0mass );
              hpi0[0][tagCount][0]->Fill( pi0mass );
              hpi0[numPhotonsCount][0][0]->Fill( pi0mass );
              hpi0[numPhotonsCount][tagCount][0]->Fill( pi0mass );

              kfit->SetFinal(kindata);
              kfit->FitTwoGammas(MASS, ERRMATRIXWEIGHT);

              if(VERBOSE) cerr << "Successful fit........" << endl;

              float prob = kfit->Prob();
              float fitpi0mass = (kfit->GetFinal_out()[0].lorentzMomentum() + kfit->GetFinal_out()[1].lorentzMomentum() ).M();

              hcl[0][0][0]->Fill(prob);
              hcl[0][tagCount][0]->Fill(prob);
              hcl[numPhotonsCount][0][0]->Fill(prob);
              hcl[numPhotonsCount][tagCount][0]->Fill(prob);

              hpi0fit[0][0][0]->Fill( fitpi0mass );
              hpi0fit[numPhotonsCount][0][0]->Fill( fitpi0mass );
              hpi0fit[0][tagCount][0]->Fill( fitpi0mass );
              hpi0fit[numPhotonsCount][tagCount][0]->Fill( fitpi0mass );

              for(int k=0;k<6;k++)
              {
                pulls[k] = kfit->GetPull(k);
                hpulls[0][0][0][k]->Fill(pulls[k]);
                hpulls[0][tagCount][0][k]->Fill(pulls[k]);
                hpulls[0][numPhotonsCount][0][k]->Fill(pulls[k]);
                hpulls[numPhotonsCount][tagCount][0][k]->Fill(pulls[k]);
              }

              for(int clCount = 1;clCount<5;clCount++)
              {
                if(prob > clCut[clCount] )
                {

                  hpi0[0][0][clCount]->Fill( pi0mass );
                  hpi0[0][tagCount][clCount]->Fill( pi0mass );
                  hpi0[numPhotonsCount][0][clCount]->Fill( pi0mass );
                  hpi0[numPhotonsCount][tagCount][clCount]->Fill( pi0mass );

                  hcl[0][0][clCount]->Fill(prob);
                  hcl[0][tagCount][clCount]->Fill(prob);
                  hcl[numPhotonsCount][0][clCount]->Fill(prob);
                  hcl[numPhotonsCount][tagCount][clCount]->Fill(prob);

                  hpi0fit[0][0][clCount]->Fill( fitpi0mass );
                  hpi0fit[numPhotonsCount][0][clCount]->Fill( fitpi0mass );
                  hpi0fit[0][tagCount][clCount]->Fill( fitpi0mass );
                  hpi0fit[numPhotonsCount][tagCount][clCount]->Fill( fitpi0mass );

                  for(int k=0;k<6;k++)
                  {
                    hpulls[0][0][clCount][k]->Fill(pulls[k]);
                    hpulls[0][tagCount][clCount][k]->Fill(pulls[k]);
                    hpulls[0][numPhotonsCount][clCount][k]->Fill(pulls[k]);
                    hpulls[numPhotonsCount][tagCount][clCount][k]->Fill(pulls[k]);
                  }

                  if(clCount == 1) numclGood++;
                }
              }
              if(VERBOSE)
              {
                cerr << "input: " << kfit->GetFinal_in()[0].energy() << " ";
                cerr << kfit->GetFinal_in()[0].px() << " ";
                cerr << kfit->GetFinal_in()[0].py() << " ";
                cerr << kfit->GetFinal_in()[0].pz() << " ";
                cerr << kfit->GetFinal_in()[1].energy() << " ";
                cerr << kfit->GetFinal_in()[1].px() << " ";
                cerr << kfit->GetFinal_in()[1].py() << " ";
                cerr << kfit->GetFinal_in()[1].pz() << endl;
                cerr << "output: " << kfit->GetFinal_out()[0].energy() << " ";
                cerr << kfit->GetFinal_out()[0].px() << " ";
                cerr << kfit->GetFinal_out()[0].py() << " ";
                cerr << kfit->GetFinal_out()[0].pz() << " ";
                cerr << kfit->GetFinal_out()[1].energy() << " ";
                cerr << kfit->GetFinal_out()[1].px() << " ";
                cerr << kfit->GetFinal_out()[1].py() << " ";
                cerr << kfit->GetFinal_out()[1].pz() << endl;
                cerr << "prob after fit: " << kfit->Prob() << endl;
                cerr << "Chi2 after fit: " << kfit->Chi2() << endl;
              }
            }
          }
          else
          {
            cerr << "OOPS!!! Not the right size kindata object!" << endl;
          }

          // Remove the last element
          kindata.pop_back();
          tagCount -= photons[j]->getTag();

        }

      }
      hclGood[0]->Fill(numclGood);
      hclGood[numPhotonsCount]->Fill(numclGood);
    }
    COUNT++;
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
