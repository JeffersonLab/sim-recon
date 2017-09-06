// Author: David Lawrence  June 25, 2004
//
//
// hd_ana.cc
//

#include "DANA/DApplication.h"
using namespace std;

#include <DVector3.h>
#include <HDGEOMETRY/DRootGeom.h>
#include <HDGEOMETRY/DGeometry.h>

#include <TROOT.h>
#include <TFile.h>
#include <TH2D.h>

int32_t RUN_NUMBER = -1;

void ParseCommandLineArguments(int narg, char *argv[]);
void Usage(string mess="");


//-----------
// main
//-----------
int main(int narg, char *argv[])
{

	ParseCommandLineArguments(narg, argv);


	// open ROOT file
	TFile *f = new TFile("material.root","RECREATE","Produced by material2root");
	cout<<"Opened ROOT file \""<<"material.root"<<"\" ..."<<endl;

	DApplication *dapp = new DApplication(narg, argv);
	DRootGeom *rg = new DRootGeom(dapp);
	DGeometry *geom = dapp->GetDGeometry(RUN_NUMBER);

	TH2D *radlen_LL = new TH2D("radlen_LL", "Radiation Length;z (cm); r(cm); rad. length (cm)", 2500, -100.0, 1170.0, 1000, 0.0, 125.0);
	TH2D *radlen_table = (TH2D*)radlen_LL->Clone("radlen_table");
	TH2D *radlen_LL_xy = new TH2D("radlen_LL_xy", "Radiation Length;x (cm); y(cm); rad. length (cm)", 900, -10.0, 10.0, 900, -10.0, 10.0);
	TH2D *radlen_table_xy = (TH2D*)radlen_LL_xy->Clone("radlen_table_xy");
	TH2D *A_LL = (TH2D*)radlen_LL->Clone("A_LL");
	TH2D *A_table = (TH2D*)radlen_LL->Clone("A_table");
	TH2D *Z_LL = (TH2D*)radlen_LL->Clone("Z_LL");
	TH2D *Z_table = (TH2D*)radlen_LL->Clone("Z_table");
	TH2D *density_LL = (TH2D*)radlen_LL->Clone("density_LL");
	density_LL->SetTitle("Density");
	density_LL->SetXTitle("z (cm)");
	density_LL->SetYTitle("r (cm)");
	density_LL->SetZTitle("density (g/cm^3)");
	TH2D *density_table = (TH2D*)radlen_LL->Clone("density_table");

	for(int ir = 1; ir<=radlen_LL->GetNbinsY(); ir++){
		double r = radlen_LL->GetYaxis()->GetBinCenter(ir);
		for(int iz = 1; iz<=radlen_LL->GetNbinsX(); iz++){
			double z = radlen_LL->GetXaxis()->GetBinCenter(iz);
			
			DVector3 pos(r, 0.0, z);
			double density, A, Z, RadLen;

			density = A = Z = RadLen = 0.0;

			rg->FindMatLL(pos, density, A, Z, RadLen);
			radlen_LL->Fill(z, r, RadLen);
			A_LL->Fill(z, r, A);
			Z_LL->Fill(z, r, Z);
			density_LL->Fill(z, r, density);

			density = A = Z = RadLen = 0.0;
			
			//rg->FindMatTable(pos, density, A, Z, RadLen);
			geom->FindMat(pos, density, A, Z, RadLen);
			radlen_table->Fill(z, r, RadLen);
			A_table->Fill(z, r, A);
			Z_table->Fill(z, r, Z);
			density_table->Fill(z, r, density);
		}
	}

	for(int ix = 1; ix<=radlen_LL_xy->GetNbinsX(); ix++){
		double x = radlen_LL_xy->GetXaxis()->GetBinCenter(ix);
		for(int iy = 1; iy<=radlen_LL_xy->GetNbinsX(); iy++){
			double y = radlen_LL_xy->GetYaxis()->GetBinCenter(iy);
			
			DVector3 pos(x, y, 65.0);
			double density, A, Z, RadLen;

			rg->FindMatLL(pos, density, A, Z, RadLen);
			radlen_LL_xy->Fill(x, y, RadLen);

			//rg->FindMatTable(pos, density, A, Z, RadLen);
			geom->FindMat(pos, density, A, Z, RadLen);
			radlen_table_xy->Fill(x, y, RadLen);
		}
	}
	
	f->Write();

	return 0;
}

//-----------
// ParseCommandLineArguments
//-----------
void ParseCommandLineArguments(int narg, char *argv[])
{
	if( narg < 2 ) Usage();
	
	for(int i=1; i<narg; i++){
		string arg = argv[i];
		string next = (i+1) < narg ? argv[i+1]:"";
		bool missing_arg = false;
		
		if( arg == "-h" || arg =="--help" ) Usage();
		if( arg == "-r" ){
			if( next.find("-") != 0){
				RUN_NUMBER = atoi(next.c_str());
			}else{ missing_arg = true; }
		}
		
		if(missing_arg) Usage( "argument " + arg + " requires and argument!" );
	}
	
	if(RUN_NUMBER<0) {
		Usage("Run number MUST be specified with -r option!");
	}
}


//-----------
// Usage
//-----------
void Usage(string mess)
{
	cout << endl;
	cout << "Usage:" << endl;
	cout << "    material2root -r RUN" << endl;
	cout << endl;
	cout << mess << endl;
	cout << endl;
	
	exit(0);
}

