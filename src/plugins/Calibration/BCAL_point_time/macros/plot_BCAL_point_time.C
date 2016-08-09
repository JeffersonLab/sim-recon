
// Plot the various histograms produced by the BCAL_point_time plugin
#include <stdio.h>


void plot_BCAL_point_time(string filename, string prefix = "05GeVgamma")
{
	char plot_filename[255];

	gStyle->SetOptStat(0);

	TDirectory *main = gDirectory;

	TFile *_file0 = new TFile(filename.c_str());
	TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("bcalpointtime");
	if(dir) dir->cd();

	TCanvas *canvas_point_NormVsZ = new TCanvas("canvas_point_NormVsZ","canvas_point_NormVsZ");
	gPad->SetGridx();
	gPad->SetGridy();
	TH1I *point_NormVsZ_layer1 = (TH1I*)gDirectory->FindObjectAny("point_NormVsZ_layer1");
	TH1I *point_NormVsZ_layer2 = (TH1I*)gDirectory->FindObjectAny("point_NormVsZ_layer2");
	TH1I *point_NormVsZ_layer3 = (TH1I*)gDirectory->FindObjectAny("point_NormVsZ_layer3");
	TH1I *point_NormVsZ_layer4 = (TH1I*)gDirectory->FindObjectAny("point_NormVsZ_layer4");
	//	TH1I *thrown_NVsZ = (TH1I*)gDirectory->FindObjectAny("thrown_NVsZ");
	point_NormVsZ_layer2->Draw();
	point_NormVsZ_layer1->Draw("same");
	point_NormVsZ_layer3->Draw("same");
	point_NormVsZ_layer4->Draw("same");
	//thrown_NVsZ->Draw("same");
	TLegend *leg = new TLegend(0.85,0.8,0.99,0.99);
	leg->AddEntry(point_NormVsZ_layer1,"Layer 1");
	leg->AddEntry(point_NormVsZ_layer2,"Layer 2");
	leg->AddEntry(point_NormVsZ_layer3,"Layer 3");
	leg->AddEntry(point_NormVsZ_layer4,"Layer 4");
	leg->Draw();
	sprintf(plot_filename,"plots/point_NormVsZ_%s.png",prefix.c_str());
	canvas_point_NormVsZ->Print(plot_filename);
	sprintf(plot_filename,"plots/point_NormVsZ_%s.pdf",prefix.c_str());
	canvas_point_NormVsZ->Print(plot_filename);

	TCanvas *canvas_point_NormVsTheta = new TCanvas("canvas_point_NormVsTheta","canvas_point_NormVsTheta");
	gPad->SetGridx();
	gPad->SetGridy();
	TH1I *point_NormVsTheta_layer1 = (TH1I*)gDirectory->FindObjectAny("point_NormVsTheta_layer1");
	TH1I *point_NormVsTheta_layer2 = (TH1I*)gDirectory->FindObjectAny("point_NormVsTheta_layer2");
	TH1I *point_NormVsTheta_layer3 = (TH1I*)gDirectory->FindObjectAny("point_NormVsTheta_layer3");
	TH1I *point_NormVsTheta_layer4 = (TH1I*)gDirectory->FindObjectAny("point_NormVsTheta_layer4");
	//	TH1I *thrown_NVsTheta = (TH1I*)gDirectory->FindObjectAny("thrown_NVsTheta");
	point_NormVsTheta_layer2->Draw();
	point_NormVsTheta_layer1->Draw("same");
	point_NormVsTheta_layer3->Draw("same");
	point_NormVsTheta_layer4->Draw("same");
	//thrown_NVsTheta->Draw("same");
	leg = new TLegend(0.85,0.8,0.99,0.99);
	leg->AddEntry(point_NormVsTheta_layer1,"Layer 1");
	leg->AddEntry(point_NormVsTheta_layer2,"Layer 2");
	leg->AddEntry(point_NormVsTheta_layer3,"Layer 3");
	leg->AddEntry(point_NormVsTheta_layer4,"Layer 4");
	leg->Draw();
	sprintf(plot_filename,"plots/point_NormVsTheta_%s.png",prefix.c_str());
	canvas_point_NormVsTheta->Print(plot_filename);
	sprintf(plot_filename,"plots/point_NormVsTheta_%s.pdf",prefix.c_str());
	canvas_point_NormVsTheta->Print(plot_filename);

	TCanvas *canvas_point_TimeVsZ = new TCanvas("canvas_point_TimeVsZ","canvas_point_TimeVsZ");
	gPad->SetGridx();
	gPad->SetGridy();
	TH1I *point_TimeVsZ_layer1_prof = (TH1I*)gDirectory->FindObjectAny("point_TimeVsZ_layer1_prof");
	TH1I *point_TimeVsZ_layer2_prof = (TH1I*)gDirectory->FindObjectAny("point_TimeVsZ_layer2_prof");
	TH1I *point_TimeVsZ_layer3_prof = (TH1I*)gDirectory->FindObjectAny("point_TimeVsZ_layer3_prof");
	TH1I *point_TimeVsZ_layer4_prof = (TH1I*)gDirectory->FindObjectAny("point_TimeVsZ_layer4_prof");
	point_TimeVsZ_layer1_prof->Draw();
	point_TimeVsZ_layer2_prof->Draw("same");
	point_TimeVsZ_layer3_prof->Draw("same");
	point_TimeVsZ_layer4_prof->Draw("same");
	leg = new TLegend(0.85,0.8,0.99,0.99);
	leg->AddEntry(point_TimeVsZ_layer1_prof,"Layer 1");
	leg->AddEntry(point_TimeVsZ_layer2_prof,"Layer 2");
	leg->AddEntry(point_TimeVsZ_layer3_prof,"Layer 3");
	leg->AddEntry(point_TimeVsZ_layer4_prof,"Layer 4");
	leg->Draw();
	sprintf(plot_filename,"plots/point_TimeVsZ_prof_%s.png",prefix.c_str());
	canvas_point_TimeVsZ->Print(plot_filename);
	sprintf(plot_filename,"plots/point_TimeVsZ_prof_%s.pdf",prefix.c_str());
	canvas_point_TimeVsZ->Print(plot_filename);

	TCanvas *canvas_hit_TimeVsZ = new TCanvas("canvas_hit_TimeVsZ","canvas_hit_TimeVsZ");
	gPad->SetGridx();
	gPad->SetGridy();
	TH1I *hitus_TimeVsZ_layer1_prof = (TH1I*)gDirectory->FindObjectAny("hitus_TimeVsZ_layer1_prof");
	TH1I *hitus_TimeVsZ_layer2_prof = (TH1I*)gDirectory->FindObjectAny("hitus_TimeVsZ_layer2_prof");
	TH1I *hitus_TimeVsZ_layer3_prof = (TH1I*)gDirectory->FindObjectAny("hitus_TimeVsZ_layer3_prof");
	TH1I *hitus_TimeVsZ_layer4_prof = (TH1I*)gDirectory->FindObjectAny("hitus_TimeVsZ_layer4_prof");
	hitus_TimeVsZ_layer1_prof->Draw();
	hitus_TimeVsZ_layer2_prof->Draw("same");
	hitus_TimeVsZ_layer3_prof->Draw("same");
	hitus_TimeVsZ_layer4_prof->Draw("same");
	TH1I *hitds_TimeVsZ_layer1_prof = (TH1I*)gDirectory->FindObjectAny("hitds_TimeVsZ_layer1_prof");
	TH1I *hitds_TimeVsZ_layer2_prof = (TH1I*)gDirectory->FindObjectAny("hitds_TimeVsZ_layer2_prof");
	TH1I *hitds_TimeVsZ_layer3_prof = (TH1I*)gDirectory->FindObjectAny("hitds_TimeVsZ_layer3_prof");
	TH1I *hitds_TimeVsZ_layer4_prof = (TH1I*)gDirectory->FindObjectAny("hitds_TimeVsZ_layer4_prof");
	hitds_TimeVsZ_layer1_prof->Draw("same");
	hitds_TimeVsZ_layer2_prof->Draw("same");
	hitds_TimeVsZ_layer3_prof->Draw("same");
	hitds_TimeVsZ_layer4_prof->Draw("same");
	leg->Draw();
	sprintf(plot_filename,"plots/hit_TimeVsZ_prof_%s.png",prefix.c_str());
	canvas_hit_TimeVsZ->Print(plot_filename);
	sprintf(plot_filename,"plots/hit_TimeVsZ_prof_%s.pdf",prefix.c_str());
	canvas_hit_TimeVsZ->Print(plot_filename);

	TCanvas *canvas_hit_TimediffVsZ = new TCanvas("canvas_TimediffVsZ","canvas_TimediffVsZ");
	gPad->SetGridx();
	gPad->SetGridy();
	TH1I *hit_TimediffVsZ_layer1_prof = (TH1I*)gDirectory->FindObjectAny("hit_TimediffVsZ_layer1_prof");
	TH1I *hit_TimediffVsZ_layer2_prof = (TH1I*)gDirectory->FindObjectAny("hit_TimediffVsZ_layer2_prof");
	TH1I *hit_TimediffVsZ_layer3_prof = (TH1I*)gDirectory->FindObjectAny("hit_TimediffVsZ_layer3_prof");
	TH1I *hit_TimediffVsZ_layer4_prof = (TH1I*)gDirectory->FindObjectAny("hit_TimediffVsZ_layer4_prof");
	hit_TimediffVsZ_layer1_prof->Draw();
	hit_TimediffVsZ_layer2_prof->Draw("same");
	hit_TimediffVsZ_layer3_prof->Draw("same");
	hit_TimediffVsZ_layer4_prof->Draw("same");
	leg->Draw();
	sprintf(plot_filename,"plots/hit_TimediffVsZ_prof_%s.png",prefix.c_str());
	canvas_hit_TimediffVsZ->Print(plot_filename);
	sprintf(plot_filename,"plots/hit_TimediffVsZ_prof_%s.pdf",prefix.c_str());
	canvas_hit_TimediffVsZ->Print(plot_filename);

	TCanvas *canvas_hit_TimesumVsZ = new TCanvas("canvas_hit_TimesumVsZ","canvas_hit_TimesumVsZ");
	gPad->SetGridx();
	gPad->SetGridy();
	TH1I *hit_TimesumVsZ_layer1_prof = (TH1I*)gDirectory->FindObjectAny("hit_TimesumVsZ_layer1_prof");
	TH1I *hit_TimesumVsZ_layer2_prof = (TH1I*)gDirectory->FindObjectAny("hit_TimesumVsZ_layer2_prof");
	TH1I *hit_TimesumVsZ_layer3_prof = (TH1I*)gDirectory->FindObjectAny("hit_TimesumVsZ_layer3_prof");
	TH1I *hit_TimesumVsZ_layer4_prof = (TH1I*)gDirectory->FindObjectAny("hit_TimesumVsZ_layer4_prof");
	hit_TimesumVsZ_layer4_prof->Draw();
	hit_TimesumVsZ_layer4_prof->SetMinimum(hit_TimesumVsZ_layer4_prof->GetMinimum(0.01)-1);
	hit_TimesumVsZ_layer1_prof->Draw("same");
	hit_TimesumVsZ_layer2_prof->Draw("same");
	hit_TimesumVsZ_layer3_prof->Draw("same");
	leg->Draw();
	sprintf(plot_filename,"plots/hit_TimesumVsZ_prof_%s.png",prefix.c_str());
	canvas_hit_TimesumVsZ->Print(plot_filename);
	sprintf(plot_filename,"plots/hit_TimesumVsZ_prof_%s.pdf",prefix.c_str());
	canvas_hit_TimesumVsZ->Print(plot_filename);


	TCanvas *canvashit_TimeVsZ_scat = new TCanvas("canvashit_TimeVsZ_scat","canvashit_TimeVsZ_scat",800,1000);
	TH1I *hitus_TimeVsZ_layer1 = (TH1I*)gDirectory->FindObjectAny("hitus_TimeVsZ_layer1");
	TH1I *hitus_TimeVsZ_layer2 = (TH1I*)gDirectory->FindObjectAny("hitus_TimeVsZ_layer2");
	TH1I *hitus_TimeVsZ_layer3 = (TH1I*)gDirectory->FindObjectAny("hitus_TimeVsZ_layer3");
	TH1I *hitus_TimeVsZ_layer4 = (TH1I*)gDirectory->FindObjectAny("hitus_TimeVsZ_layer4");
	TH1I *hitds_TimeVsZ_layer1 = (TH1I*)gDirectory->FindObjectAny("hitds_TimeVsZ_layer1");
	TH1I *hitds_TimeVsZ_layer2 = (TH1I*)gDirectory->FindObjectAny("hitds_TimeVsZ_layer2");
	TH1I *hitds_TimeVsZ_layer3 = (TH1I*)gDirectory->FindObjectAny("hitds_TimeVsZ_layer3");
	TH1I *hitds_TimeVsZ_layer4 = (TH1I*)gDirectory->FindObjectAny("hitds_TimeVsZ_layer4");
	canvashit_TimeVsZ_scat->Divide(2,4,0.001,0.001);
	canvashit_TimeVsZ_scat->cd(1);
	hitus_TimeVsZ_layer1->Draw("colz");
	canvashit_TimeVsZ_scat->cd(2);
	hitus_TimeVsZ_layer2->Draw("colz");
	canvashit_TimeVsZ_scat->cd(3);
	hitus_TimeVsZ_layer3->Draw("colz");
	canvashit_TimeVsZ_scat->cd(4);
	hitus_TimeVsZ_layer4->Draw("colz");
	canvashit_TimeVsZ_scat->cd(5);
	hitds_TimeVsZ_layer1->Draw("colz");
	canvashit_TimeVsZ_scat->cd(6);
	hitds_TimeVsZ_layer2->Draw("colz");
	canvashit_TimeVsZ_scat->cd(7);
	hitds_TimeVsZ_layer3->Draw("colz");
	canvashit_TimeVsZ_scat->cd(8);
	hitds_TimeVsZ_layer4->Draw("colz");
	sprintf(plot_filename,"plots/hit_TimeVsZ_%s.png",prefix.c_str());
	canvashit_TimeVsZ_scat->Print(plot_filename);
	sprintf(plot_filename,"plots/hit_TimeVsZ_%s.pdf",prefix.c_str());
	canvashit_TimeVsZ_scat->Print(plot_filename);


	TCanvas *canvas_hit_TimediffVsZ_scat = new TCanvas("canvas_hit_TimediffVsZ_scat","canvas_hit_TimediffVsZ_scat",800,600);
	TH1I *hit_TimediffVsZ_layer1 = (TH1I*)gDirectory->FindObjectAny("hit_TimediffVsZ_layer1");
	TH1I *hit_TimediffVsZ_layer2 = (TH1I*)gDirectory->FindObjectAny("hit_TimediffVsZ_layer2");
	TH1I *hit_TimediffVsZ_layer3 = (TH1I*)gDirectory->FindObjectAny("hit_TimediffVsZ_layer3");
	TH1I *hit_TimediffVsZ_layer4 = (TH1I*)gDirectory->FindObjectAny("hit_TimediffVsZ_layer4");
	canvas_hit_TimediffVsZ_scat->Divide(2,2,0.001,0.001);
	canvas_hit_TimediffVsZ_scat->cd(1);
	hit_TimediffVsZ_layer1->Draw("colz");
	canvas_hit_TimediffVsZ_scat->cd(2);
	hit_TimediffVsZ_layer2->Draw("colz");
	canvas_hit_TimediffVsZ_scat->cd(3);
	hit_TimediffVsZ_layer3->Draw("colz");
	canvas_hit_TimediffVsZ_scat->cd(4);
	hit_TimediffVsZ_layer4->Draw("colz");
	sprintf(plot_filename,"plots/hit_TimediffVsZ_%s.png",prefix.c_str());
	canvas_hit_TimediffVsZ_scat->Print(plot_filename);
	sprintf(plot_filename,"plots/hit_TimediffVsZ_%s.pdf",prefix.c_str());
	canvas_hit_TimediffVsZ_scat->Print(plot_filename);

	TCanvas *canvas_hit_TimesumVsZ_scat = new TCanvas("canvas_hit_TimesumVsZ_scat","canvas_hit_TimesumVsZ_scat",800,600);
	TH1I *hit_TimesumVsZ_layer1 = (TH1I*)gDirectory->FindObjectAny("hit_TimesumVsZ_layer1");
	TH1I *hit_TimesumVsZ_layer2 = (TH1I*)gDirectory->FindObjectAny("hit_TimesumVsZ_layer2");
	TH1I *hit_TimesumVsZ_layer3 = (TH1I*)gDirectory->FindObjectAny("hit_TimesumVsZ_layer3");
	TH1I *hit_TimesumVsZ_layer4 = (TH1I*)gDirectory->FindObjectAny("hit_TimesumVsZ_layer4");
	canvas_hit_TimesumVsZ_scat->Divide(2,2,0.001,0.001);
	canvas_hit_TimesumVsZ_scat->cd(1);
	hit_TimesumVsZ_layer1->Draw("colz");
	canvas_hit_TimesumVsZ_scat->cd(2);
	hit_TimesumVsZ_layer2->Draw("colz");
	canvas_hit_TimesumVsZ_scat->cd(3);
	hit_TimesumVsZ_layer3->Draw("colz");
	canvas_hit_TimesumVsZ_scat->cd(4);
	hit_TimesumVsZ_layer4->Draw("colz");
	sprintf(plot_filename,"plots/hit_TimesumVsZ_%s.png",prefix.c_str());
	canvas_hit_TimesumVsZ_scat->Print(plot_filename);
	sprintf(plot_filename,"plots/hit_TimesumVsZ_%s.pdf",prefix.c_str());
	canvas_hit_TimesumVsZ_scat->Print(plot_filename);


	// back to main dir
	main->cd();

}
