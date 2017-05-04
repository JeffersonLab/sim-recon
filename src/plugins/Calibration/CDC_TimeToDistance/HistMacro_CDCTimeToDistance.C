double long_drift_func[3][3];
double short_drift_func[3][3];
double magnet_correction[2][2];

vector<double> cdc_drift_table;

// Set values for the region cut
float deltaMin = -0.175, deltaMax = 0.175, tMin = 300, tMax = 1200;

unsigned int Locate(vector<double>&xx, double x){
   int n=xx.size();
   if (x==xx[0]) return 0;
   else if (x==xx[n-1]) return n-2;

   int jl=-1;
   int ju=n;
   int ascnd=(xx[n-1]>=xx[0]);
   while(ju-jl>1){
      int jm=(ju+jl)>>1;
      if ( (x>=xx[jm])==ascnd)
         jl=jm;
      else
         ju=jm;
   }
   return jl;
}

Double_t TimeToDistance( Double_t *x, Double_t *par){
   Double_t d=0.0;
   double delta = x[1]; // yAxis
   double t = x[0]; // xAxis

   // Cut out region in  fit.
   if (delta > deltaMax || delta < deltaMin) return 0.0;
   if (delta < (((deltaMax - deltaMin) / (tMax - tMin))*(t - tMin) + deltaMin)) return 0.0;
   // Variables to store values for time-to-distance functions for delta=0
   // and delta!=0
   double f_0=0.;
   double f_delta=0.;
   if (t > 0){
      if (delta>0){
         double a1=par[0];
         double a2=par[1];
         double a3=par[2];
         double b1=par[3];
         double b2=par[4];
         double b3=par[5];
         double c1=par[6];
         double c2=par[7];
         double c3=par[8];
         // use "long side" functional form
         double my_t=0.001*t;
         double sqrt_t=sqrt(my_t);
         double t3=my_t*my_t*my_t;
         double delta_mag=fabs(delta);
         double a=a1+a2*delta_mag;
         double b=b1+b2*delta_mag;
         double c=c1+c2*delta_mag+c3*delta*delta;
         f_delta=a*sqrt_t+b*my_t+c*t3;
         f_0=a1*sqrt_t+b1*my_t+c1*t3;

      }
      else{
         double my_t=0.001*t;
         double sqrt_t=sqrt(my_t);
         double delta_mag=fabs(delta);

         // use "short side" functional form
         double a1=par[9];
         double a2=par[10];
         double a3=par[11];
         double b1=par[12];
         double b2=par[13];
         double b3=par[14];
         double c1=par[15];
         double c2=par[16];
         double c3=par[17];
         double delta_sq=delta*delta;
         double a=a1+a2*delta_mag+a3*delta_sq;
         double b=b1+b2*delta_mag+b3*delta_sq;
         f_delta=a*sqrt_t+b*my_t;
         f_0=a1*sqrt_t+b1*my_t;

      }

      unsigned int max_index=cdc_drift_table.size()-1;
      if (t>cdc_drift_table[max_index]){
         d=f_delta;
         return d;
      }
      // Drift time is within range of table -- interpolate...
      unsigned int index=0;
      index=Locate(cdc_drift_table,t);
      double dt=cdc_drift_table[index+1]-cdc_drift_table[index];
      double frac=(t-cdc_drift_table[index])/dt;
      double d_0=0.01*(double(index)+frac);

      double P=0.;
      double tcut=250.0; // ns
      if (t<tcut) {
         P=(tcut-t)/tcut;
      }
      d=f_delta*(d_0/f_0*P+1.-P);
   }
   return d;
}


void HistMacro_CDCTimeToDistance(){
   TF2 *f = new TF2("f",TimeToDistance, 10, 1500, -0.3, 0.3, 18);

   TProfile *constants = (TProfile *) gDirectory->Get("/CDC_TimeToDistance/CDC_TD_Constants");
   long_drift_func[0][0] = constants->GetBinContent(101);
   long_drift_func[0][1] = constants->GetBinContent(102);
   long_drift_func[0][2] = constants->GetBinContent(103);
   long_drift_func[1][0] = constants->GetBinContent(104);
   long_drift_func[1][1] = constants->GetBinContent(105);
   long_drift_func[1][2] = constants->GetBinContent(106);
   long_drift_func[2][0] = constants->GetBinContent(107);
   long_drift_func[2][1] = constants->GetBinContent(108);
   long_drift_func[2][2] = constants->GetBinContent(109);
   magnet_correction[0][0] = constants->GetBinContent(110);
   magnet_correction[0][1] = constants->GetBinContent(111);
   short_drift_func[0][0] = constants->GetBinContent(112);
   short_drift_func[0][1] = constants->GetBinContent(113);
   short_drift_func[0][2] = constants->GetBinContent(114);
   short_drift_func[1][0] = constants->GetBinContent(115);
   short_drift_func[1][1] = constants->GetBinContent(116);
   short_drift_func[1][2] = constants->GetBinContent(117);
   short_drift_func[2][0] = constants->GetBinContent(118);
   short_drift_func[2][1] = constants->GetBinContent(119);
   short_drift_func[2][2] = constants->GetBinContent(120);
   magnet_correction[1][0] = constants->GetBinContent(121);
   magnet_correction[1][1] = constants->GetBinContent(122);

   for (unsigned int i=1; i<=78; i++){
      cdc_drift_table.push_back(constants->GetBinContent(i));
   }

   Double_t parameters[18] =
   {long_drift_func[0][0], long_drift_func[0][1], long_drift_func[0][2],
      long_drift_func[1][0], long_drift_func[1][1], long_drift_func[1][2],
      long_drift_func[2][0], long_drift_func[2][1], long_drift_func[2][2],
      short_drift_func[0][0], short_drift_func[0][1], short_drift_func[0][2],
      short_drift_func[1][0], short_drift_func[1][1], short_drift_func[1][2],
      short_drift_func[2][0], short_drift_func[2][1], short_drift_func[2][2]};

   f->SetParameters(parameters);

   TCanvas *c = new TCanvas("c_TD","CDC TD", 1200, 550);
   c->Divide(2,1);

   c->cd(1);

   TProfile2D *profile = (TProfile2D *) gDirectory->Get("/CDC_TimeToDistance/Predicted Drift Distance Vs Delta Vs t_drift");

   gStyle->SetOptStat(0);

   Double_t contours[21] =
   { 0.00, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
      0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};

   profile->SetContour(21, contours);
   profile->Draw("colz");
   profile->GetXaxis()->SetRangeUser(0,1200);
   profile->GetYaxis()->SetRangeUser(-0.22,0.22);
   profile->GetYaxis()->SetTitleOffset(1.15);
   f->SetContour(21, contours);
   profile->Draw("cont2 list same");
   f->Draw("cont2 list same");

   c->cd(2);
   TH2I *hist = (TH2I *) gDirectory->Get("/CDC_TimeToDistance/Residual Vs. Drift Time"); 
   hist->Draw("colz");
   hist->GetXaxis()->SetRangeUser(0,1200);
   hist->GetYaxis()->SetTitleOffset(1.40);
   c->cd(2)->SetGridx();
   c->cd(2)->SetGridy();

   c->cd();
   c->Update();
}
