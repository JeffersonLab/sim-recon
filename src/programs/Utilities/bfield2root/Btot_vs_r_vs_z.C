//
//
// This needs to be compiled rather than run via cint.
//
// root [0] .x Btot_vs_r_vs_z.C+
//

//#include "StandardLabels.C"
#include "GlueX_boundaries.C"


void Btot_vs_r_vs_z(void)
{
	gROOT->Reset();
	//gStyle->SetPadRightMargin(0.15);
	
	TCanvas *c1 = new TCanvas("c1");
	c1->SetTicks();
	c1->SetGrid();
	c1->SetRightMargin(0.125);
	
	new TFile("bfield.root");
	TH2D *Btot_vs_r_vs_z = (TH2D*)gROOT->FindObject("Btot_vs_r_vs_z");

	TH2D *axes = new TH2D("axes", "B-field magnitude", 100, -100.0, 450.0, 100, 0.0, 200.0);
	axes->SetXTitle("z (cm)");
	axes->SetYTitle("r (cm)");
	axes->SetStats(0);
	axes->GetZaxis()->SetRangeUser(-0.1, 3.0);
	axes->Draw();
	
	Btot_vs_r_vs_z->GetZaxis()->SetRangeUser(-0.1, 3.0);
	Btot_vs_r_vs_z->GetZaxis()->SetTitle("Tesla");
	Btot_vs_r_vs_z->Draw("same zcol");

	DrawGlueXBoundaries();
	//StandardLabels1D(grad,"solenoid_1500_poisson_20090814_01");
	
	// Save
	c1->SaveAs("Btot_vs_r_vs_z.pdf");
	c1->SaveAs("Btot_vs_r_vs_z.gif");
	
}



