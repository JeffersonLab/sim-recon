//
//
// This needs to be compiled rather than run via cint.
//
// root [0] .x dBtot_vs_r_vs_z.C+
//

//#include "StandardLabels.C"
#include "GlueX_boundaries.C"


void dBtot_vs_r_vs_z(void)
{
	gROOT->Reset();
	//gStyle->SetPadRightMargin(0.15);
	
	TCanvas *c1 = new TCanvas("c1");
	c1->SetTicks();
	c1->SetGrid();
	c1->SetRightMargin(0.125);
	
	new TFile("bfield.root");
	TH2D *dBtot_vs_r_vs_z = (TH2D*)gROOT->FindObject("dBtot_vs_r_vs_z");

	// Convert from Tesla/cm to Gauss/cm
	dBtot_vs_r_vs_z->Scale(10000.0);
	
	TH2D *axes = new TH2D("axes", "B-field Gradient magnitude", 100, -100.0, 450.0, 100, 0.0, 200.0);
	axes->SetXTitle("z (cm)");
	axes->SetYTitle("r (cm)");
	axes->SetStats(0);
	axes->GetZaxis()->SetRangeUser(0.0, 1000.0);
	axes->Draw();
	
	dBtot_vs_r_vs_z->GetZaxis()->SetRangeUser(0.0, 1000.0);
	dBtot_vs_r_vs_z->GetZaxis()->SetTitle("Gauss/cm");
	dBtot_vs_r_vs_z->Draw("same zcol");

	DrawGlueXBoundaries();
	//StandardLabels1D(grad,"solenoid_1500_poisson_20090814_01");
	
	// Save
	c1->SaveAs("dBtot_vs_r_vs_z.pdf");
	c1->SaveAs("dBtot_vs_r_vs_z.gif");
}



