// hnamepath: /highlevel/EventVertexZ
// hnamepath: /highlevel/EventVertexYVsX

//EventVertexZ
//EventVertexYVsX

{
	TDirectory *locTopDirectory = gDirectory;

	//Goto Beam Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("highlevel");
	if(!locDirectory)
		return;
	locDirectory->cd();

	TH1* locHist_EventVertexZ = (TH1*)gDirectory->Get("EventVertexZ");
	TH2* locHist_EventVertexYVsX = (TH2*)gDirectory->Get("EventVertexYVsX");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("Vertex", "Vertex", 1200, 600); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 1);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_EventVertexZ != NULL)
	{
		locHist_EventVertexZ->GetXaxis()->SetTitleSize(0.05);
		locHist_EventVertexZ->GetYaxis()->SetTitleSize(0.05);
		locHist_EventVertexZ->GetXaxis()->SetLabelSize(0.05);
		locHist_EventVertexZ->GetYaxis()->SetLabelSize(0.035);
		locHist_EventVertexZ->Draw();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_EventVertexYVsX != NULL)
	{
		locHist_EventVertexYVsX->GetXaxis()->SetTitleSize(0.05);
		locHist_EventVertexYVsX->GetYaxis()->SetTitleSize(0.045);
		locHist_EventVertexYVsX->GetXaxis()->SetLabelSize(0.05);
		locHist_EventVertexYVsX->GetYaxis()->SetLabelSize(0.05);
		locHist_EventVertexYVsX->Draw("colz");
		gPad->SetLogz();
		gPad->Update();
	}
}
