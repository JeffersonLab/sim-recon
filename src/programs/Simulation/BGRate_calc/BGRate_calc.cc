#include <iostream>
#include <sstream>
using namespace std;
#include <string>
#include <stdio.h>
#include <cstdlib>
#include <ctype.h>

#include <typeinfo>
#include "AMPTOOLS_MCGEN/CobremsGeneration.hh"
#include "TFile.h"
#include "TH1D.h"

// /w/halld-scifs1a/home/scole/gluex_top/sim-recon/master/Linux_CentOS7-x86_64-gcc4.8.5/bin is location
// of the executable

static void show_usage(string argStr)
{
	cout<<endl<<endl;
	cout<<"Usage, print, and write arguemtns:"<<endl;
	cout<<"\t-h or --help"<<endl;
	cout<<"\t\tDefault: false"<<endl;
	cout<<"\t-w or --write"<<endl;
        cout<<"\t\tDefault: false"<<endl;
        cout<<"\t\tNote: Only use for saving the histogram to a root file for"<<endl;
        cout<<"\t\tstudying the generation beyond the integration value."<<endl;
	cout<<endl;

	cout<<"Necessary arguments to MCWrapper, otherwise optional:"<<endl;
	cout<<endl;
	cout<<"\t--coherent_peak <coherent peak>"<<endl;
	cout<<"\t\tDefault: 9.0 GeV"<<endl;
	cout<<"\t\tNote: If set to 0.0, will treat it as an amorphous run."<<endl;
	cout<<"\t--beam_on_current <beam on current>"<<endl;
	cout<<"\t\tDefault: 0.01 uA"<<endl;
	cout<<"\t--beam_energy <beam energy>"<<endl;
	cout<<"\t\tDefault: 12.0 GeV"<<endl;
	cout<<"\t--collimator_diameter <collimator diameter>"<<endl;
	cout<<"\t\tDefautl: 0.005 m"<<endl;
	cout<<"\t--radiator_thickness <radiator thickness>"<<endl;
	cout<<"\t\tDefault: 20 um"<<endl;
	cout<<"\t--endpoint_energy_low <endpoint energy low>"<<endl;
	cout<<"\t\tDefault: 3.0 GeV"<<endl;
	cout<<"\t\tNote: Low end of the endpoint energy cannot go to 0 due to"<<endl; 
	cout<<"\t\tvertical asymptote in photon eneryg spectrum at 0. Set to"<<endl;
	cout<<"\t\tlow end of the tagger."<<endl;
	cout<<"\t--endpoint_energy_high <endpoint energy high>"<<endl;
	cout<<"\t\tDefault: 12.0 GeV"<<endl;

	cout<<endl;
	cout<<"Optional arguments, subject to change for MCWrapper:"<<endl;
	cout<<endl;
	cout<<"\t--runNo <run number>"<<endl;
	cout<<"\t\tDefault: 0"<<endl;
	cout<<"\t--beam_energy_rms <beam energy rms>"<<endl;
	cout<<"\t\tDefault: 6.0e-4 GeV"<<endl;
	cout<<"\t--beam_emittance <beam emittance>"<<endl;
	cout<<"\t\t Default: 10.0e-9"<<endl;
	cout<<"\t--collimator_distance <collimator distance>"<<endl;
	cout<<"\t\tDefault: 76.0 m"<<endl;
	cout<<"\t--nbins <nbins>"<<endl;
	cout<<"\t\tDefault: 2000"<<endl;
	cout<<"\t\tNote: Will be removed in the future when binning selection decided."<<endl;
	cout<<endl;
}

double getTargetRadiationLength_Schiff(double Z, double N, double a)
{
	const double me = 0.510998910e-3;
	const double alpha = 7.2973525698e-3;
	const double hbarc = 0.1973269718e-15;

	double zeta = log(1440 * pow(Z, -2/3.)) / log(183 * pow(Z, -1/3.));
	double s = 4 * N * pow(alpha, 3) * pow(hbarc/(a*me), 2) / a * Z * (Z + zeta) * log(183 * pow(Z, -1/3.));
	return 1/s;
}

int main(const int argc, char* argv[])
{
	if (argc == 1)
	{
		show_usage(argv[0]);
		return 1;
	}

	int runNo = 0;
	bool write = false;
	double coherent_peak = 9.0;
        double beam_on_current = 0.01;
	double beam_energy = 12.0;
	double beam_energy_rms = 6.0e-4;
	double beam_emittance = 10.0e-09;
	double collimator_distance = 76.0;
	double collimator_diameter = 0.005;
	double radiator_thickness = 20e-6;
	int nbins = 2000;
	double photon_energy_min = 0.0;
	double endpoint_energy_low = 3.0;
	double endpoint_energy_high = 12.0;

	for(int i = 1; i < argc; i++)
	{
		string arg = argv[i];
		if ((arg == "-h") || (arg == "--help"))
		{
			show_usage(argv[0]);
			return 0;
		}

		else if ((arg == "-w") || (arg == "--write"))
                {
                        write = true;
                        continue;
                }

		i++;
		
		if (arg == "--runNo")
		{
			runNo = atoi(argv[i]); 
		}

		else if (arg == "--coherent_peak")
		{
			coherent_peak = atof(argv[i]);
		}

		else if (arg == "--beam_energy")
		{
			beam_energy = atof(argv[i]);
		}

		else if (arg == "--beam_on_current")
		{
			beam_on_current = atof(argv[i]);
		}

		else if (arg == "--beam_energy_rms")
		{
			beam_energy_rms = atof(argv[i]);
		}
		
		else if (arg == "--beam_emittance")
		{
			beam_emittance = atof(argv[i]);
		}

		else if (arg == "--collimator_distance")
		{
			collimator_distance = atof(argv[i]);
		}

		else if (arg == "--collimator_diameter")
		{
			collimator_diameter = atof(argv[i]);
		}

		else if (arg == "--radiator_thickness")
		{
			radiator_thickness = atof(argv[i]);
		}

		else if (arg == "--nbins")
		{
			nbins = atoi(argv[i]);
		}

		else if (arg == "--photon_energy_min")
		{
			photon_energy_min = atof(argv[i]);
		}

		else if (arg == "--endpoint_energy_low")
		{
			endpoint_energy_low = atof(argv[i]);
		}

		else if (arg == "--endpoint_energy_high")
		{
			endpoint_energy_high = atof(argv[i]);
		}
	}

	cout<<"Building off given information below"<<endl;
	cout<<"Run number: "<<runNo<<endl;
        cout<<"Coherent peak: "<<coherent_peak<<endl;
	cout<<"Beam energy: "<<beam_energy<<endl;
	cout<<"Beam energy rms: "<<beam_energy_rms<<endl;
	cout<<"Beam curret: "<<beam_on_current<<endl;
	cout<<"Beam emittance: "<<beam_emittance<<endl;
	cout<<"Collimator diameter: "<<collimator_diameter<<endl;
	cout<<"Collimator distance: "<<collimator_distance<<endl;
	cout<<"Radiator thickness: "<<radiator_thickness<<endl;
	cout<<"Minimum photon energy: "<<photon_energy_min<<endl;
	cout<<"Endpoint energy low: "<<endpoint_energy_low<<endl;
	cout<<"Endpoint energy high: "<<endpoint_energy_high<<endl;
	cout<<"Number of bins: "<<nbins<<endl;
	
	CobremsGeneration* cobrems = new CobremsGeneration(beam_energy,coherent_peak);
	cobrems->setBeamEnergy(beam_energy);
	cobrems->setCoherentEdge(coherent_peak);
	cobrems->setBeamErms(beam_energy_rms);
	cobrems->setBeamEmittance(beam_emittance);
	cobrems->setCollimatorDistance(collimator_distance);
	cobrems->setCollimatorDiameter(collimator_diameter);
	cobrems->setTargetThickness(radiator_thickness);
	cobrems->setCollimatedFlag(true);

	cobrems->printBeamlineInfo();

	double Emin = endpoint_energy_low;
	double Emax = endpoint_energy_high;

	double x0 = Emin / beam_energy;
	double x1 = Emax / beam_energy;

	double xvals[nbins];
	double yvals[nbins];

	if (coherent_peak > 0.0) 
	{
		cout<<"Polarized BGRate"<<endl;
		for(int i=0; i < nbins; i++)
		{
			xvals[i] = x0 + (i +0.5) * (x1 - x0) / nbins;
			yvals[i] = cobrems->Rate_dNtdx(xvals[i]) * beam_on_current / 1.6e-13;
		}
	}
	else
	{
		cout<<"Amorphous BGRate"<<endl;
		for(int i=0; i < nbins; i++)
                {
                        xvals[i] = x0 + (i +0.5) * (x1 - x0) / nbins;
                        yvals[i] = cobrems->Rate_dNidx(xvals[i]) * beam_on_current / 1.6e-13;
                }
	}

	cobrems->applyBeamCrystalConvolution(nbins, xvals, yvals);

        TH1D* dRtdkH1 = new TH1D("dRtdkH1", "", nbins, Emin, Emax);
        dRtdkH1->GetXaxis()->SetRangeUser(Emin + (Emax - Emin)/10., Emax);

	for(int i=0; i < nbins; i++)
	{
		dRtdkH1->Fill(xvals[i]*beam_energy, yvals[i]/beam_energy);
	}

	double persec = (Emax - Emin) * 1./nbins;
	double erate = dRtdkH1->Integral(dRtdkH1->FindBin(endpoint_energy_low), dRtdkH1->FindBin(endpoint_energy_high) - 1) * persec;

	if (coherent_peak == 0.0)
	{
		erate = erate * getTargetRadiationLength_Schiff(13, 4, 404.95e-12) / cobrems->getTargetRadiationLength_Schiff();
	}

	if (write == true)
        {
		char charBuff[50];
	        sprintf(charBuff, "\\mbox{photon beam spectrum vs }E_\\gamma \\mbox{ (/GeV/s)} runNo: %i", runNo);
		dRtdkH1->SetTitle(charBuff);
	        sprintf(charBuff, "BGRate_%i.root", runNo);
		cout<<"Saving file"<<endl;
		TFile* f = new TFile(charBuff,"recreate");
		dRtdkH1->Write();
		f->Write();
		f->Close();
        }


	cout<<"BGRate GHz = "<<erate*pow(10,-9)<<endl;

	return 0;
}
