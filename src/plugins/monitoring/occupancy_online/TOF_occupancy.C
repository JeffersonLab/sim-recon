
// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//
// hnamepath: /occupancy/tdcOccS
// hnamepath: /occupancy/tdcOccN
// hnamepath: /occupancy/tdcOccU
// hnamepath: /occupancy/tdcOccD
// hnamepath: /occupancy/adcOccS
// hnamepath: /occupancy/adcOccN
// hnamepath: /occupancy/adcOccU
// hnamepath: /occupancy/adcOccD

{
	TDirectory *savedir = gDirectory;

	TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("occupancy");
	if(dir) dir->cd();

	TH1I *ha_s = (TH1I*)gDirectory->FindObjectAny("adcOccS");
	TH1I *ha_n = (TH1I*)gDirectory->FindObjectAny("adcOccN");
	TH1I *ha_u = (TH1I*)gDirectory->FindObjectAny("adcOccU");
	TH1I *ha_d = (TH1I*)gDirectory->FindObjectAny("adcOccD");

	TH1I *h_s = (TH1I*)gDirectory->FindObjectAny("tdcOccS");
	TH1I *h_n = (TH1I*)gDirectory->FindObjectAny("tdcOccN");
	TH1I *h_u = (TH1I*)gDirectory->FindObjectAny("tdcOccU");
	TH1I *h_d = (TH1I*)gDirectory->FindObjectAny("tdcOccD");

	TLegend *legend_sa = new TLegend(0.1,0.85,0.3,0.9);
	TLegend *legend_na = new TLegend(0.3,0.85,0.5,0.9);
	TLegend *legend_ua = new TLegend(0.5,0.85,0.7,0.9);
	TLegend *legend_da = new TLegend(0.7,0.85,0.9,0.9);

	TLegend *legend_s = new TLegend(0.1,0.85,0.3,0.9);
	TLegend *legend_n = new TLegend(0.3,0.85,0.5,0.9);
	TLegend *legend_u = new TLegend(0.5,0.85,0.7,0.9);
	TLegend *legend_d = new TLegend(0.7,0.85,0.9,0.9);

	ha_s->SetBarWidth(0.5);
	ha_s->SetBarOffset(0);
	ha_s->SetFillColor(2);
	ha_s->SetStats(0);
	ha_s->SetXTitle("Module number");
	ha_s->SetYTitle("fADC occupancy");
	ha_s->SetTitleSize(0.05,"X");
	ha_s->GetXaxis()->CenterTitle();
	ha_s->SetTitleSize(0.05,"Y");
	ha_s->GetYaxis()->CenterTitle();

	ha_n->SetBarWidth(0.5);
	ha_n->SetBarOffset(0.5);
	ha_n->SetFillColor(3);
	ha_n->SetStats(0);

	ha_u->SetBarWidth(0.5);
	ha_u->SetBarOffset(1.0);
	ha_u->SetFillColor(4);
	ha_u->SetStats(0);

	ha_d->SetBarWidth(0.5);
	ha_d->SetBarOffset(1.5);
	ha_d->SetFillColor(6);
	ha_d->SetStats(0);

	h_s->SetBarWidth(0.5);
	h_s->SetBarOffset(0);
	h_s->SetFillColor(2);
	h_s->SetStats(0);
	h_s->SetXTitle("Module number");
	h_s->SetYTitle("TDC occupancy");
	h_s->SetTitleSize(0.05,"X");
	h_s->GetXaxis()->CenterTitle();
	h_s->SetTitleSize(0.05,"Y");
	h_s->GetYaxis()->CenterTitle();

	h_n->SetBarWidth(0.5);
	h_n->SetBarOffset(0.5);
	h_n->SetFillColor(3);
	h_n->SetStats(0);

	h_u->SetBarWidth(0.5);
	h_u->SetBarOffset(1.0);
	h_u->SetFillColor(4);
	h_u->SetStats(0);

	h_d->SetBarWidth(0.5);
	h_d->SetBarOffset(1.5);
	h_d->SetFillColor(6);
	h_d->SetStats(0);

	legend_sa->AddEntry(ha_s,"South","f");
	legend_na->AddEntry(ha_n,"North","f");
	legend_ua->AddEntry(ha_u,"Up","f");
	legend_da->AddEntry(ha_d,"Down","f");

	legend_s->AddEntry(h_s,"South","f");
	legend_n->AddEntry(h_n,"North","f");
	legend_u->AddEntry(h_u,"Up","f");
	legend_d->AddEntry(h_d,"Down","f");


	// Just for testing
	if(gPad == NULL){
		TCanvas *c1 = new TCanvas("c1");
		c1->cd(0);
		c1->Draw();
		c1->Update();
	}
	if(!gPad) {savedir->cd(); return;}

	TCanvas *c1 = gPad->GetCanvas();
	c1->cd(0);
	c1->Clear();

	c1->Divide(1,2);
	TVirtualPad *pad1 = c1->cd(1);
	pad1->SetTicks();
	pad1->SetGridy();
	if(ha_s) ha_s->Draw("BAR");
	if(ha_n) ha_n->Draw("BAR sames");
	if(ha_u) ha_u->Draw("BAR sames");
	if(ha_d) ha_d->Draw("BAR sames");

	legend_sa->Draw();
	legend_na->Draw();
	legend_ua->Draw();
	legend_da->Draw();

	TVirtualPad *pad2 = c1->cd(2);
	pad2->SetTicks();
	pad2->SetGridy();
	if(h_s) h_s->Draw("BAR");
	if(h_n) h_n->Draw("BAR sames");
	if(h_u) h_u->Draw("BAR sames");
	if(h_d) h_d->Draw("BAR sames");

	legend_s->Draw();
	legend_n->Draw();
	legend_u->Draw();
	legend_d->Draw();

	savedir->cd();
}
