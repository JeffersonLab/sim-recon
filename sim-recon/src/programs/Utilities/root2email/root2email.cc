// $Id$
//
// Author: David Lawrence  April 7, 2006
//
//
// root2email.cc
//

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>
#include <stdlib.h>
using namespace std;

#include <TCanvas.h>
#include <TFrame.h>
#include <TFile.h>
#include <TH2F.h>
#include <TROOT.h>

const char *filename = "mctrk_ana.root";

typedef struct{
	string message;
	string gif_file;
}plot_t;

vector<string> email_address;
vector<string> histnames;
vector<plot_t> plots;

void ParseCommandLineArguments(int narg, char *argv[]);
void Usage(void);
void AddPlot(TFile &f, string what);

int colors[] = {kRed, kBlue, kGreen};
int Ncolors = 3;

//-----------
// main
//-----------
int main(int narg, char *argv[])
{

	ParseCommandLineArguments(narg, argv);

	// open ROOT file for reading
	TFile f(filename);
	if(!f.IsOpen()){
		cout<<"Can't open file \""<<filename<<"\" !"<<endl;
		return 0;
	}
	cout<<"Opened ROOT file \""<<f.GetName()<<"\""<<endl;
	
	// Generate plots
	for(unsigned int i=0; i<histnames.size(); i++){
		AddPlot(f, histnames[i]);
	}
	
	// Create command to send e-mail (excluding the e-mail address)
	stringstream cmd;
	cmd<<"(echo \"Auto-generated from "<<filename<<":\"; ";
	cmd<<"echo \" \"; echo \" \"; ";
	for(unsigned int i=0;i<plots.size();i++){
		cmd<<"echo \""<<plots[i].message<<"\"; ";
		cmd<<"uuencode "<<plots[i].gif_file<<" "<<plots[i].gif_file<<"; ";
	}
	if(plots.size()==0)cmd<<"echo \"No histograms specified!!\"; ";
	cmd<<") | mail -s \"ROOT plots\" ";
	
	// Loop over e-mail addresses, sending the message to each one
	for(unsigned int i=0; i<email_address.size(); i++){
		stringstream my_cmd;
		my_cmd<<cmd.str()<<email_address[i];
		cout<<my_cmd.str()<<endl;
		system(my_cmd.str().c_str());
	}
	
	
	return 0;
}

//------------------------------
// AddPlot
//------------------------------
void AddPlot(TFile &f, string histname)
{
	plot_t plot;
	stringstream mess;

	TH1F *hist=(TH1F*)gROOT->FindObject(histname.c_str());
	//f.GetObject(histname.c_str(), hist);
	if(!hist){
		mess<<"Unable to read in histogram '"<<histname<<"' !";
		plot.gif_file = string("");
	}else{
		// This would be soooo much better if one of the output
		// binary formats(gif, jpg, png ..) worked right here.
		// As it is, the output only looks correct for vector-oriented
		// outputs (eps, pdf ,...) so we have to convert it externally.
		
		// Generate EPS file using ROOT
		TCanvas c1("c1",histname.c_str(),200,10,700,500);
   	c1.SetGrid();
		c1.SetFillColor(19);
		hist->SetLineColor(colors[plots.size()%Ncolors]);
		hist->Draw();
		c1.Print("tmp.eps");
		
		// Convert EPS to PPM using ghostscript(gs)
		system("echo quit | gs -sDEVICE=ppm -r72x72 -g565x405 -sOutputFile=tmp.ppm -dNOPAUSE tmp.eps > /dev/null");

		// Generate unique name of GIF file
		stringstream outfile;
		outfile<<"plot"<<plots.size()<<".gif";
		plot.gif_file = outfile.str();

		// Convert PPM to GIF using ppmtogif
		stringstream cmd;
		cmd<<"ppmtogif tmp.ppm > "<<outfile.str();
		system(cmd.str().c_str());
		
		// Add histogram name to message
		mess<<histname;
		
		// Delete temporary files
		system("rm -f tmp.ppm tmp.eps");
	}
	
	// Add plot to list. We do this even if we failed to generate a plot
	// so the e-mail can have an error message
	plot.message = mess.str();
	plots.push_back(plot);
}

//------------------------------
// ParseCommandLineArguments
//------------------------------
void ParseCommandLineArguments(int narg, char *argv[])
{
	if(narg<2){
		Usage();
	}

	for(int i=1; i<narg; i++){
		if(string(argv[i]) == "-cf"){
		}else if(string(argv[i]) == "-r"){
			if(i<narg-1){
				email_address.push_back(argv[i+1]);
				i++;
			}
		}else if(string(argv[i]) == "-H"){
			if(i<narg-1){
				histnames.push_back(argv[i+1]);
				i++;
			}
		}else if(string(argv[i]) == "-if"){
			if(i<narg-1){
				filename = argv[i+1];
				i++;
			}
		}else if(string(argv[i]) == "-h"){
			Usage();
		}else{
		
		}
	}
	
	if(email_address.size() == 0){
		cout<<endl;
		cout<<"You MUST provide at least one e-mail address!!"<<endl<<endl;
		Usage();
	}
}

//------------------------------
// Usage
//------------------------------
void Usage(void)
{
	cout<<endl;
	cout<<"Usage:"<<endl;
	cout<<"      root2email [options]"<<endl;
	cout<<endl;
	cout<<"   -r email-address     recipient e-mail address"<<endl;
	cout<<"                        (can be given multiple times)"<<endl;
	cout<<"   -H histogram         name of histogram. If histogram"<<endl;
	cout<<"                        is in a TDirectory, include the"<<endl;
	cout<<"                        path. (can be given multiple times)"<<endl;
	cout<<"   -if rootfile         name of ROOT file to use"<<endl;
	cout<<"   -h                   print this help message."<<endl;
	
	cout<<endl;
	exit(0);
}

