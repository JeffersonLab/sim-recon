// Author: David Lawrence  June 25, 2004
//
//
// hd_ana.cc
//

#include "DANA/DApplication.h"
using namespace std;

#include <DVector3.h>
#include <HDGEOMETRY/DRootGeom.h>

#include <TROOT.h>
#include <TFile.h>
#include <TH2D.h>

void Usage(JApplication &app);


//-----------
// main
//-----------
int main(int narg, char *argv[])
{
	// open ROOT file
	TFile *f = new TFile("material.root","RECREATE","Produced by hd_ana");
	cout<<"Opened ROOT file \""<<"material.root"<<"\" ..."<<endl;

	DApplication *dapp = new DApplication(narg, argv);
	DRootGeom *rg = new DRootGeom(dapp);

	TH2D *radlen_LL = new TH2D("radlen_LL", "Radiation Length", 500, -100.0, 650.0, 500, 0.0, 125.0);
	TH2D *radlen_table = (TH2D*)radlen_LL->Clone("radlen_table");
	TH2D *A_LL = (TH2D*)radlen_LL->Clone("A_LL");
	TH2D *A_table = (TH2D*)radlen_LL->Clone("A_table");
	TH2D *Z_LL = (TH2D*)radlen_LL->Clone("Z_LL");
	TH2D *Z_table = (TH2D*)radlen_LL->Clone("Z_table");
	TH2D *density_LL = (TH2D*)radlen_LL->Clone("density_LL");
	TH2D *density_table = (TH2D*)radlen_LL->Clone("density_table");
	
	for(int ir = 1; ir<=radlen_LL->GetNbinsX(); ir++){
		double r = radlen_LL->GetYaxis()->GetBinCenter(ir);
		for(int iz = 1; iz<=radlen_LL->GetNbinsY(); iz++){
			double z = radlen_LL->GetXaxis()->GetBinCenter(iz);
			
			DVector3 pos(r, 0.0, z);
			double density, A, Z, RadLen;

			rg->FindMatLL(pos, density, A, Z, RadLen);
			radlen_LL->Fill(z, r, RadLen);
			A_LL->Fill(z, r, A);
			Z_LL->Fill(z, r, Z);
			density_LL->Fill(z, r, density);

			rg->FindMatTable(pos, density, A, Z, RadLen);
			radlen_table->Fill(z, r, RadLen);
			A_table->Fill(z, r, A);
			Z_table->Fill(z, r, Z);
			density_table->Fill(z, r, density);
		}
	}
	
	f->Write();
	
	return 0;
}

//-----------
// Usage
//-----------
void Usage(JApplication &app)
{
	cout<<endl;
	cout<<"Usage:"<<endl;
	cout<<"    material2root"<<endl;
	cout<<endl;
	app.Usage();
	cout<<endl;
	
	exit(0);
}

