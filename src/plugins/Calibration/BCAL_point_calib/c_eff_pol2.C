#include <TRandom.h>

void GetCCDBConstants(TString path, Int_t run, TString variation, vector<double>& container, Int_t column = 1){
  char command[256];

  sprintf(command, "ccdb dump %s:%i:%s", path.Data(), run, variation.Data());
  FILE* inputPipe = gSystem->OpenPipe(command, "r");
  if(inputPipe == NULL)
    return;
  //get the first (comment) line
  char buff[1024];
  if(fgets(buff, sizeof(buff), inputPipe) == NULL)
    return;
  //get the remaining lines
  double entry;
  int counter = 0;
  while(fgets(buff, sizeof(buff), inputPipe) != NULL){
    istringstream locConstantsStream(buff);
    while (locConstantsStream >> entry){
      counter++;
      if (counter % column == 0) container.push_back(entry);
    }
  }
  //Close the pipe
  gSystem->ClosePipe(inputPipe);
}


// macro that gets TGraphs from the ROOT file generated from the BCAL_point_calibs plugin and returns the effective velocities for each channel, using a quadratic fit
void c_eff_pol2(TString fileName = "hd_root.root", int runNumber = 22016, TString variation = "default")
{
        ofstream outfile("c_eff_pol2.out");

        //gROOT->Reset();

	// change this to the location and name of your own ROOT file
        TFile *in = TFile::Open( fileName , "READ");
	if (in == 0) {
	  cout << "Unable to open file " << fileName.Data() << "...Exiting" << endl;
	  return;
	}

	// histogram to hold the effective velocities
        TH1D* h1_c_eff2 = NULL;

	// histogram to hold the p2 fit parameter
        TH1D* h1_p2_fit_parameter = NULL;

	// function to hold the fit results
        TF1* func2 = NULL;

	// graph for storing the graphs returned from the ROOT file
        TGraph* h2_tgraph = NULL;

        h1_c_eff2 = new TH1D("h1_c_eff2","Effective Velocity per Channel",800,0,800);
        h1_c_eff2->SetXTitle("Channel #");
        h1_c_eff2->SetYTitle("Effective Velocity (cm/ns)");
        h1_c_eff2->GetYaxis()->SetRangeUser(15,19);
        h1_c_eff2->SetMarkerStyle(20);
        h1_c_eff2->SetOption("E1");


        h1_p2_fit_parameter = new TH1D("h1_p2_fit_parameter", "p2 fit parameter", 800, 0, 800);
        h1_p2_fit_parameter->SetXTitle("Channel #");
        h1_p2_fit_parameter->SetYTitle("p2 fit parameter (cm^{-1)");
        h1_p2_fit_parameter->SetMarkerStyle(20);
        h1_p2_fit_parameter->SetOption("E1");

        int channels = 768;
	
	vector<double> effective_velocity;
	GetCCDBConstants("/BCAL/effective_velocities",runNumber, variation, effective_velocity);

	outfile.open(prefix + "effective_velocities.txt");

        // variables for saving the parameters of the fit
        float eff_velocities_graphs = 0;
        float eff_velocities_graphs_errors = 0;
        float p2_graphs = 0;
        float p1_graphs = 0;
        float p2_graphs_error = 0;
        float p1_graphs_error = 0;
        char string[20];

	// loop over the channels, apply a quadratic fit for each one and calculate the effective velocities
        for(int m = 0; m < channels; ++m){

                sprintf(string,"bcal_point_calibs/h2_tgraph[%d]", m+1);

                h2_tgraph = (TGraph*)in->Get(string);

                if(h2_tgraph->GetN() > 2){
			// quadratic fit
                        h2_tgraph->Fit("pol2");
                        func2 = h2_tgraph->GetFunction("pol2");
                        p2_graphs = func2->GetParameter(2);
                        p2_graphs_error = func2->GetParError(2);
                        p1_graphs = func2->GetParameter(1);
                        p1_graphs_error = func2->GetParError(1);
                        eff_velocities_graphs = effective_velocity[m] / (p1_graphs);
                        eff_velocities_graphs_errors = effective_velocity[m] * p1_graphs_error / ( p1_graphs * p1_graphs );

			outfile << eff_velocities_graphs << endl;
                }

		// fill the histograms
                h1_c_eff2->Fill(m+1,eff_velocities_graphs);
                h1_c_eff2->SetBinError(m+2,eff_velocities_graphs_errors);
                h1_p2_fit_parameter->Fill(m+1, p2_graphs);
                h1_p2_fit_parameter->SetBinError(m+2, p2_graphs_error);
	
		h2_tgraph->Delete();

        } // end of loop over channels


        TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",323,71,953,603);
        Canvas_1->Range(-100,12.25,900,19.75);
        Canvas_1->SetFillColor(0);
        Canvas_1->SetBorderMode(0);
        Canvas_1->SetBorderSize(2);
        Canvas_1->SetFrameBorderMode(0);
        Canvas_1->SetFrameBorderMode(0);


        h1_c_eff2->Draw("");

	// save results
        Canvas_1->SaveAs("h1_c_eff2.png");
        Canvas_1->SaveAs("h1_c_eff2.C");

        TCanvas *Canvas_4 = new TCanvas("Canvas_4", "Canvas_4",323,71,953,603);
        Canvas_4->Range(-100,12.25,900,19.75);
        Canvas_4->SetFillColor(0);
        Canvas_4->SetBorderMode(0);
        Canvas_4->SetBorderSize(2);
        Canvas_4->SetFrameBorderMode(0);
        Canvas_4->SetFrameBorderMode(0);

        h1_p2_fit_parameter->Draw("");

	// save results
        Canvas_4->SaveAs("h1_p2_fit_parameter.png");
        Canvas_4->SaveAs("h1_p2_fit_parameter.C");

	in->Close();
	outfile.close();
}
