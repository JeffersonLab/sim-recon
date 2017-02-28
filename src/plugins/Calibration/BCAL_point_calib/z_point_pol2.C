
// macro that gets TGraphs from the ROOT file generated from the BCAL_point_calibs plugin and returns the effective velocities for each channel, using a quadratic fit
void z_point_pol2(TString fileName = "hd_root.root", int runNumber = 22016, TString variation = "default")
{
        ofstream outfile("z_point_pol2.out");

        //gROOT->Reset();

	// change this to the location and name of your own ROOT file
        TFile *in = TFile::Open( fileName , "READ");
	if (in == 0) {
	  cout << "Unable to open file " << fileName.Data() << "...Exiting" << endl;
	  return;
	}

	// function to hold the fit results
        TF1* func2 = NULL;

	// graph for storing the graphs returned from the ROOT file
        TGraph* h2_tgraph = NULL;


        int channels = 768;
	
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

                sprintf(string,"bcal_point_calibs/h2_tgraph_deltat[%d]", m+1);

                h2_tgraph = (TGraph*)in->Get(string);

                if(h2_tgraph->GetN() > 2){
			// quadratic fit
                        h2_tgraph->Fit("pol2");
                        func2 = h2_tgraph->GetFunction("pol2");
                        //p2_graphs = func2->GetParameter(2);
                        //p2_graphs_error = func2->GetParError(2);
                        //p1_graphs = func2->GetParameter(1);
                        //p1_graphs_error = func2->GetParError(1);

			outfile << func2->GetParameter(0) << " "  << func2->GetParameter(1) << " "  << func2->GetParameter(2) << endl;
                }

		//h2_tgraph->Delete();

        } // end of loop over channels


	in->Close();
	outfile.close();
}
