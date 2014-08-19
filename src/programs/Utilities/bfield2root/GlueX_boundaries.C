//-----------------------------------------------------------------
// This macro is to be included and called from an external file 
// to draw the various GlueX detectors on a TH2D. Call this
// macro with:
//
// DrawGlueXBoundaries([full_detector=false], [fill_in=false],
//                     [clip_TOF=1000], [color=-1 (multi-colored]);
// full_detector: boolean which toggles drawing the detectors
//                above or both above and below the beamline.
//       fill_in: boolean which toggles between drawing outlines
//                or filled shapes
//      clip_TOF: double which will reduce the size of the TOF
//                detector in the r-direction
//         color: Int_t which will draw all detectors with the
//                specified color
//
// A minimal example of an macro calling GlueX_boundaries:
//-----------------------------------------------------------------
//#include "GlueX_boundaries.C"
#include <TH2D.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TFile.h>

#if 0
void boundaries_plot()
{
	TColor::CreateColorWheel();

	TCanvas *c1 = new TCanvas("boundaries canvas");
	c1->SetTicks();

	Int_t zmin = -100;
	Int_t zmax=700, rmax=200;
	TH2D *boundaries = new TH2D("boundaries", "The GlueX Detectors",
                              zmax-zmin, zmin, zmax, rmax*2, -rmax, rmax);
	boundaries->SetStats(0);
	boundaries->SetXTitle("Z (cm)");
	boundaries->SetYTitle("R (cm)");
	
	boundaries->Draw();
	DrawGlueXBoundaries(true, true);
}
#endif
//-----------------------------------------------------------------

#include <TBox.h>
#include <TLatex.h>
#include <TPolyLine.h>

Int_t full=false;   // used to toggle between drawing full detector or only above the beamline
Int_t fill_detectors=false;   // used to toggle outlined or full color detectors
Int_t target_color = kBlack;
Int_t start_counter_color = kGreen+3;
Int_t FDC_color = kBlue;
Int_t CDC_color = kPink+6;
Int_t TOF_color = kGreen+2;
Int_t BCAL_color = kOrange+1;
Int_t FCAL_color = kOrange+2;
Int_t Magnet_color = kBlack;

void DrawPoly(int N, double *z, double *r, int color);

//--------------------------
// DrawTargetBoundaries
//--------------------------
void DrawTargetBoundaries(int color=target_color)
{
  double Zlo = 50.0;
  double Zhi = Zlo + 30.0;
  double Rlo = 0.0;
  double Rhi = Rlo + 1.5;
  if (full) { Rlo=-1.5; Rhi=Rlo+3; }

  TBox *box = new TBox(Zlo, Rlo, Zhi, Rhi);
  box->SetLineWidth(2.0);
  if (fill_detectors) { box->SetFillColor(color); box->SetFillStyle(1001); }
  else box->SetFillStyle(0);
  box->SetLineColor(color);
  box->Draw();

  TLatex *lab = new TLatex(Zlo-5, 0, "target");
  if (full) lab->SetTextAlign(32); else lab->SetTextAlign(31);
  lab->SetTextSize(0.02);
  lab->SetTextColor(color);
  lab->Draw();
}

//--------------------------
// DrawStartCounterBoundaries
//--------------------------
void DrawStartCounterBoundaries(int color=start_counter_color)
{
  // Values from Richard's old spreadsheet "start_geom.xls"
  const int Npoints = 15;
  double r_hi[] = {8.675, 8.675, 8.675, 7.759, 7.354, 4.353, 3.134,
          2.290, 1.915, 5.593, 5.893, 6.907, 6.951, 6.951, 8.675};
  double r_lo[Npoints];
  for (Int_t i=0; i<Npoints; i++) r_lo[i] = -r_hi[i];
  
  double z[] = {   0.0, 51.383, 51.726, 53.938, 54.914, 57.915, 59.134,
          58.442, 58.290, 54.360, 54.238, 51.531, 51.488,    0.0, 0.0};
  
  // shift z to proper location in lab system
  for(int i=0; i<Npoints; i++)z[i] += 38.75;
  
  TPolyLine *pol = new TPolyLine(Npoints, z, r_hi);
  pol->SetLineWidth(2.0);
  if (fill_detectors) { pol->SetFillColor(color); pol->SetFillStyle(1001); }
  else pol->SetFillStyle(0);
  pol->SetLineColor(color);
  pol->Draw();

  TLatex *lab = new TLatex(z[6]+(z[6]-z[0])*0.1, 0, "start counter");
  if (full) lab->SetTextAlign(12); else lab->SetTextAlign(11);
  lab->SetTextSize(0.02);
  lab->SetTextColor(color);
  lab->Draw();

  if (full) {
    pol = new TPolyLine(Npoints, z, r_lo);
    pol->SetLineWidth(2.0);
    if (fill_detectors) { pol->SetFillColor(color); pol->SetFillStyle(1001); }
    else pol->SetFillStyle(0);
    pol->SetLineColor(color);
    pol->Draw();
  }
}

//--------------------------
// DrawFDCBoundaries
//--------------------------
void DrawFDCBoundaries(int color=FDC_color)
{
  double Rlo = 2.0, Rhi=48.5;

  double Zlo[7], Zhi[7];

  Zlo[0] = 176.1586;
  Zhi[0] = 187.1614;
  
  Zlo[1] = 233.7186;
  Zhi[1] = 244.7214;

  Zlo[2] = 291.2786;
  Zhi[2] = 302.2814;

  Zlo[3] = 348.8386;
  Zhi[3] = 359.8414;

  for(int sl=0; sl<4; sl++){
  
    TBox *box_hi = new TBox(Zlo[sl], Rlo, Zhi[sl], Rhi);
    box_hi->SetLineWidth(2.0);
    if (fill_detectors) { box_hi->SetFillColor(color); box_hi->SetFillStyle(1001); }
    else box_hi->SetFillStyle(0);
    box_hi->SetLineColor(color);
    box_hi->Draw();

    TLatex *lab_hi = new TLatex((Zlo[sl]+Zhi[sl])/2.0+2.0, Rhi+1, "FDC");
    lab_hi->SetTextAlign(21);
    lab_hi->SetTextSize(0.025);
    lab_hi->SetTextColor(color);
    lab_hi->Draw();
    
    if (full) {
      TBox *box_lo = new TBox(Zlo[sl], -Rlo, Zhi[sl], -Rhi);
      box_lo->SetLineWidth(2.0);
      if (fill_detectors) { box_lo->SetFillColor(color); box_lo->SetFillStyle(1001); }
      else box_lo->SetFillStyle(0);
      box_lo->SetLineColor(color);
      box_lo->Draw();
    }
  }       
}

//--------------------------
// DrawCDCBoundaries
//--------------------------
void DrawCDCBoundaries(int color=CDC_color)
{
  double Zlo=17.0, Zhi=Zlo+150.0;
  double Rlo[7], Rhi[7];

  Rlo[0] = 10.7219;
  Rhi[0] = 15.1621;
  
  Rlo[1] = 16.9321;
  Rhi[1] = sqrt(pow(21.8912,2.0) + pow(0.860106,2.0));

  Rlo[2] = 23.8544;
  Rhi[2] = sqrt(pow(28.5658,2.0) + pow(0.846871,2.0));

  Rlo[3] = 31.3799;
  Rhi[3] = 35.8301;

  Rlo[4] = 37.4446;
  Rhi[4] = sqrt(pow(41.9225,2.0) + pow(0.833676,2.0));

  Rlo[5] = 43.6152;
  Rhi[5] = sqrt(pow(48.0733,2.0) + pow(0.829899,2.0));

  Rlo[6] = 50.3747;
  Rhi[6] = 54.7617;
  
  TBox *box_hi = new TBox(Zlo, Rlo[0], Zhi, Rhi[6]);
  box_hi->SetLineWidth(2.0);
  if (fill_detectors) { box_hi->SetFillColor(color); box_hi->SetFillStyle(1001); }
  else box_hi->SetFillStyle(0);
  box_hi->SetLineColor(color);
  box_hi->Draw();
  
  TLatex *lab_hi = new TLatex(Zlo+5, Rhi[6]-1, "CDC");
  lab_hi->SetTextAlign(13);
  lab_hi->SetTextSize(0.03);
  if (fill_detectors) lab_hi->SetTextColor(kWhite); else lab_hi->SetTextColor(color);
  lab_hi->Draw();

  if (full) {
    TBox *box_lo = new TBox(Zlo, -Rlo[0], Zhi, -Rhi[6]);
    box_lo->SetLineWidth(2.0);
    if (fill_detectors) { box_lo->SetFillColor(color); box_lo->SetFillStyle(1001); }
    else box_lo->SetFillStyle(0);
    box_lo->SetLineColor(color);
    box_lo->Draw();
  
    TLatex *lab_lo = new TLatex(Zlo+5, -Rlo[0]-1, "CDC");
    lab_lo->SetTextAlign(13);
    lab_lo->SetTextSize(0.03);
    if (fill_detectors) lab_lo->SetTextColor(kWhite); else lab_lo->SetTextColor(color);
    lab_lo->Draw();
  }
}

//--------------------------
// DrawTOFBoundaries
//--------------------------
void DrawTOFBoundaries(int color=TOF_color, double clip_tof=1000.0)
{
  double Zlo = 616.25;
  double Zhi = Zlo + 2.0*2.54;
  double Rlo = 3.0;
  double Rhi = 252.0/2.0;
  
  // Clip TOF at top of histogram area
  if(Rhi>clip_tof)Rhi = clip_tof;
  
  TBox *box_hi = new TBox(Zlo, Rlo, Zhi, Rhi);
  box_hi->SetLineWidth(2.0);
  if (fill_detectors) { box_hi->SetFillColor(color); box_hi->SetFillStyle(1001); }
  else box_hi->SetFillStyle(0);
  box_hi->SetLineColor(color);
  box_hi->Draw();

  TLatex *lab_hi = new TLatex(Zlo-5, Rhi-1, "TOF");
  lab_hi->SetTextAlign(31);
  lab_hi->SetTextSize(0.03);
  lab_hi->SetTextAngle(90.0);
  lab_hi->SetTextColor(color);
  lab_hi->Draw();

  if (full) {
    TBox *box_lo = new TBox(Zlo, -Rlo, Zhi, -Rhi);
    box_lo->SetLineWidth(2.0);
    box_lo->SetLineColor(color);
    if (fill_detectors) { box_lo->SetFillColor(color); box_lo->SetFillStyle(1001);}
    else box_lo->SetFillStyle(0);
    box_lo->Draw();
  }
}

//--------------------------
// DrawBCALBoundaries
//--------------------------
void DrawBCALBoundaries(int color=BCAL_color)
{
  double Zlo = 17.0;
  double Zhi = Zlo + 390.0;
  double Rlo = 64.2; // includes aluminum plate
  double Rhi = 90.0;
  
  TBox *box_hi = new TBox(Zlo, Rlo, Zhi, Rhi);
  box_hi->SetLineWidth(2.0);
  if (fill_detectors) { box_hi->SetFillColor(color); box_hi->SetFillStyle(1001); }
  else box_hi->SetFillStyle(0);
  box_hi->SetLineColor(color);
  box_hi->Draw();
  
  TLatex *lab_hi = new TLatex(Zlo+5, Rhi-1, "BCAL");
  lab_hi->SetTextAlign(13);
  lab_hi->SetTextSize(0.03);
  if (fill_detectors) lab_hi->SetTextColor(kWhite); else lab_hi->SetTextColor(color);
  lab_hi->Draw();

  if (full) {
    TBox *box_lo = new TBox(Zlo, -Rlo, Zhi, -Rhi);
    box_lo->SetLineWidth(2.0);
    if (fill_detectors) { box_lo->SetFillColor(color); box_lo->SetFillStyle(1001); }
    else box_lo->SetFillStyle(0);
    box_lo->SetLineColor(color);
    box_lo->Draw();
  
    TLatex *lab_lo = new TLatex(Zlo+5, -Rlo-1, "BCAL");
    lab_lo->SetTextAlign(13);
    lab_lo->SetTextSize(0.03);
    if (fill_detectors) lab_lo->SetTextColor(kWhite); else lab_lo->SetTextColor(color);
    lab_lo->Draw();
  }
}

//--------------------------
// DrawBCALBoundaries
//--------------------------
void DrawFCALBoundaries(int color=FCAL_color)
{
  double Zlo = 622.8;
  double Zhi = Zlo + 45.;
  double Rlo = 6.0;
  double Rhi = 106.;
  
  TBox *box_hi = new TBox(Zlo, Rlo, Zhi, Rhi);
  box_hi->SetLineWidth(2.0);
  if (fill_detectors) { box_hi->SetFillColor(color); box_hi->SetFillStyle(1001); }
  else box_hi->SetFillStyle(0);
  box_hi->SetLineColor(color);
  box_hi->Draw();

  TLatex *lab_hi = new TLatex(Zlo, Rhi-1, "FCAL");
  lab_hi->SetTextAlign(13);
  lab_hi->SetTextSize(0.025);
  if (fill_detectors) lab_hi->SetTextColor(kWhite); else lab_hi->SetTextColor(color);
  lab_hi->Draw();

  if (full) {
    TBox *box_lo = new TBox(Zlo, -Rlo, Zhi, -Rhi);
    box_lo->SetLineWidth(2.0);
    box_lo->SetLineColor(color);
    if (fill_detectors) { box_lo->SetFillColor(color); box_lo->SetFillStyle(1001); }
    else box_lo->SetFillStyle(0);
    box_lo->Draw();
  }
}

//--------------------------
// DrawMagnet
//--------------------------
void DrawMagnet(int yoke_color=Magnet_color, int coil_color=Magnet_color)
{
	// Yoke
	int N = 0;
	double z[50], r[50];
	z[N] =   0.000; r[N++] =   89.540;
	z[N] =   0.000; r[N++] =  147.320;
	z[N] =  53.820; r[N++] =  147.320;
	z[N] =  53.820; r[N++] =  150.500;
	z[N] =  56.200; r[N++] =  150.500;
	z[N] =  56.200; r[N++] =  154.940;
	z[N] =  60.960; r[N++] =  154.940;
	z[N] =  60.960; r[N++] =  147.320;
	z[N] = 187.800; r[N++] =  147.320;
	z[N] = 187.800; r[N++] =  150.500;
	z[N] = 190.180; r[N++] =  150.500;
	z[N] = 190.180; r[N++] =  154.940;
	z[N] = 194.950; r[N++] =  154.940;
	z[N] = 194.950; r[N++] =  147.320;
	z[N] = 264.000; r[N++] =  147.320;
	z[N] = 264.000; r[N++] =  150.500;
	z[N] = 266.380; r[N++] =  150.500;
	z[N] = 266.380; r[N++] =  154.940;
	z[N] = 271.150; r[N++] =  154.940;
	z[N] = 271.150; r[N++] =  147.320;
	z[N] = 340.200; r[N++] =  147.320;
	z[N] = 340.200; r[N++] =  150.500;
	z[N] = 342.580; r[N++] =  150.500;
	z[N] = 342.580; r[N++] =  154.940;
	z[N] = 347.350; r[N++] =  154.940;
	z[N] = 347.350; r[N++] =  147.320;
	z[N] = 354.970; r[N++] =  147.320;
	z[N] = 354.970; r[N++] =   92.710;
	z[N] = 428.630; r[N++] =   92.710;
	z[N] = 428.630; r[N++] =  187.960;
	z[N] = -50.800; r[N++] =  187.960;
	z[N] = -50.800; r[N++] =   89.540;
	z[N] =   0.000; r[N++] =   89.540;
	DrawPoly(N, z, r, yoke_color);
	 
	// Baffles
	N = 0;
	z[N] =   76.200; r[N++] =  147.320;
	z[N] =   76.200; r[N++] =   92.710;
	z[N] =   84.460; r[N++] =   92.710;
	z[N] =   84.460; r[N++] =  147.320;
	z[N] =   76.200; r[N++] =  147.320;
	DrawPoly(N, z, r, yoke_color);

	N = 0;
	z[N] =  202.570; r[N++] =  147.320;
	z[N] =  202.570; r[N++] =   92.710;
	z[N] =  210.190; r[N++] =   92.710;
	z[N] =  210.190; r[N++] =  147.320;
	z[N] =  202.570; r[N++] =  147.320;
	DrawPoly(N, z, r, yoke_color);

	// Coils
	//&reg mat=1, cur= 432000.0; ; Coil 1A 
	N = 0;
	z[N] =   92.853; r[N++] =  101.881;
	z[N] =  104.548; r[N++] =  101.881;
	z[N] =  104.548; r[N++] =  116.952;
	z[N] =   92.853; r[N++] =  116.952;
	z[N] =   92.853; r[N++] =  101.881;
	DrawPoly(N, z, r, coil_color);

	//&reg mat=1, cur= 414000.0; ; Coil 1B 
	N = 0;
	z[N] =  104.548; r[N++] =  101.881;
	z[N] =  116.511; r[N++] =  101.881;
	z[N] =  116.511; r[N++] =  116.320;
	z[N] =  104.548; r[N++] =  116.320;
	z[N] =  104.548; r[N++] =  101.881;
	DrawPoly(N, z, r, coil_color);

	//&reg mat=1, cur= 180000.0; ; Coil 1C 
	N = 0;
	z[N] =  124.961; r[N++] =  101.881;
	z[N] =  132.679; r[N++] =  101.881;
	z[N] =  132.679; r[N++] =  111.263;
	z[N] =  124.961; r[N++] =  111.263;
	z[N] =  124.961; r[N++] =  101.881;
	DrawPoly(N, z, r, coil_color);

	//&reg mat=1, cur= 120000.0; ; Coil 1D 
	N = 0;
	z[N] =  140.845; r[N++] =  101.881;
	z[N] =  148.563; r[N++] =  101.881;
	z[N] =  148.563; r[N++] =  108.607;
	z[N] =  140.845; r[N++] =  108.607;
	z[N] =  140.845; r[N++] =  101.881;
	DrawPoly(N, z, r, coil_color);

	//&reg mat=1, cur= 324000.0; ; Coil 1E 
	N = 0;
	z[N] =  156.729; r[N++] =  101.881;
	z[N] =  168.423; r[N++] =  101.881;
	z[N] =  168.423; r[N++] =  113.159;
	z[N] =  156.729; r[N++] =  113.159;
	z[N] =  156.729; r[N++] =  101.881;
	DrawPoly(N, z, r, coil_color);

	//&reg mat=1, cur= 324000.0; ; Coil 1F 
	N = 0;
	z[N] =  172.416; r[N++] =  101.881;
	z[N] =  180.134; r[N++] =  101.881;
	z[N] =  180.134; r[N++] =  118.849;
	z[N] =  172.416; r[N++] =  118.849;
	z[N] =  172.416; r[N++] =  101.881;
	DrawPoly(N, z, r, coil_color);

	//&reg mat=1, cur= 348000.0; ; Coil 1G 
	N = 0;
	z[N] =  180.134; r[N++] =  101.881;
	z[N] =  188.087; r[N++] =  101.881;
	z[N] =  188.087; r[N++] =  120.113;
	z[N] =  180.134; r[N++] =  120.113;
	z[N] =  180.134; r[N++] =  101.881;
	DrawPoly(N, z, r, coil_color);

	//&reg mat=1, cur= 468000.0; ; Coil 2A 
	N = 0;
	z[N] =    7.256; r[N++] =  101.881;
	z[N] =   18.950; r[N++] =  101.881;
	z[N] =   18.950; r[N++] =  118.217;
	z[N] =    7.256; r[N++] =  118.217;
	z[N] =    7.256; r[N++] =  101.881;
	DrawPoly(N, z, r, coil_color);

	//&reg mat=1, cur= 204000.0; ; Coil 2B 
	N = 0;
	z[N] =   22.124; r[N++] =  101.881;
	z[N] =   29.842; r[N++] =  101.881;
	z[N] =   29.842; r[N++] =  112.527;
	z[N] =   22.124; r[N++] =  112.527;
	z[N] =   22.124; r[N++] =  101.881;
	DrawPoly(N, z, r, coil_color);

	//&reg mat=1, cur= 324000.0; ; Coil 2C 
	N = 0;
	z[N] =   29.842; r[N++] =  101.881;
	z[N] =   41.772; r[N++] =  101.881;
	z[N] =   41.772; r[N++] =  112.527;
	z[N] =   29.842; r[N++] =  112.527;
	z[N] =   29.842; r[N++] =  101.881;
	DrawPoly(N, z, r, coil_color);

	//&reg mat=1, cur= 396000.0; ; Coil 2D 
	N = 0;
	z[N] =   42.007; r[N++] =  101.881;
	z[N] =   53.702; r[N++] =  101.881;
	z[N] =   53.702; r[N++] =  115.688;
	z[N] =   42.007; r[N++] =  115.688;
	z[N] =   42.007; r[N++] =  101.881;
	DrawPoly(N, z, r, coil_color);

	//&reg mat=1, cur= 468000.0; ; Coil 3A 
	N = 0;
	z[N] =  217.441; r[N++] =  101.881;
	z[N] =  229.135; r[N++] =  101.881;
	z[N] =  229.135; r[N++] =  118.217;
	z[N] =  217.441; r[N++] =  118.217;
	z[N] =  217.441; r[N++] =  101.881;
	DrawPoly(N, z, r, coil_color);

	//&reg mat=1, cur= 306000.0; ; Coil 3B 
	N = 0;
	z[N] =  232.309; r[N++] =  101.881;
	z[N] =  244.003; r[N++] =  101.881;
	z[N] =  244.003; r[N++] =  112.527;
	z[N] =  232.309; r[N++] =  112.527;
	z[N] =  232.309; r[N++] =  101.881;
	DrawPoly(N, z, r, coil_color);

	//&reg mat=1, cur= 192000.0; ; Coil 3C 
	N = 0;
	z[N] =  244.003; r[N++] =  101.881;
	z[N] =  251.957; r[N++] =  101.881;
	z[N] =  251.957; r[N++] =  112.527;
	z[N] =  244.003; r[N++] =  112.527;
	z[N] =  244.003; r[N++] =  101.881;
	DrawPoly(N, z, r, coil_color);

	//&reg mat=1, cur= 198000.0; ; Coil 3D 
	N = 0;
	z[N] =  251.957; r[N++] =  101.881;
	z[N] =  263.887; r[N++] =  101.881;
	z[N] =  263.887; r[N++] =  108.734;
	z[N] =  251.957; r[N++] =  108.734;
	z[N] =  251.957; r[N++] =  101.881;
	DrawPoly(N, z, r, coil_color);

	//&reg mat=1, cur= 324000.0; ; Coil 4A/B 
	N = 0;
	z[N] =  293.929; r[N++] =  101.881;
	z[N] =  305.623; r[N++] =  101.881;
	z[N] =  305.623; r[N++] =  113.159;
	z[N] =  293.929; r[N++] =  113.159;
	z[N] =  293.929; r[N++] =  101.881;
	DrawPoly(N, z, r, coil_color);

	//&reg mat=1, cur=1890000.0; ; Coil 4C/D 
	N = 0;
	z[N] =  310.709; r[N++] =  101.881;
	z[N] =  340.299; r[N++] =  101.881;
	z[N] =  340.299; r[N++] =  128.331;
	z[N] =  310.709; r[N++] =  128.331;
	z[N] =  310.709; r[N++] =  101.881;
	DrawPoly(N, z, r, coil_color);

	// Label
	TLatex *lab_hi = new TLatex(-48.0, 180.0, "solenoid");
	lab_hi->SetTextAlign(12);
	lab_hi->SetTextSize(0.03);
	lab_hi->SetTextColor(coil_color);
	lab_hi->Draw();

}

//--------------------------
// DrawPoly
//--------------------------
void DrawPoly(int N, double *z, double *r, int color)
{
	TPolyLine *pl = new TPolyLine(N, z, r);
	pl->SetLineColor(color);
	if (fill_detectors) { pl->SetFillColor(color); pl->SetFillStyle(1001); }
	pl->SetLineWidth(2);
	pl->Draw();
	
	 if (full) {
		for(int i=0; i<N; i++)r[i] = -r[i];
		pl = new TPolyLine(N, z, r);
		pl->SetLineColor(color);
		if (fill_detectors) { pl->SetFillColor(color); pl->SetFillStyle(1001); }
		pl->SetLineWidth(2);
		pl->Draw();
	 }
}

//--------------------------
// DrawGlueXBoundaries
//--------------------------
void DrawGlueXBoundaries(Int_t full_detector=false, Int_t fill_in=false, Double_t clip_tof=1000.0, Int_t color=-1)
{
  full = full_detector;
  fill_detectors = fill_in;
  if (color==-1) DrawBCALBoundaries(); else DrawBCALBoundaries(color);
  if (color==-1) DrawCDCBoundaries(); else DrawCDCBoundaries(color);
  if (color==-1) DrawFDCBoundaries(); else DrawFDCBoundaries(color);
  if (color==-1) DrawStartCounterBoundaries(); else DrawStartCounterBoundaries(color);
  if (color==-1) DrawTOFBoundaries(); else DrawTOFBoundaries(color, clip_tof);
  if (color==-1) DrawTargetBoundaries(); else DrawTargetBoundaries(color);
  if (color==-1) DrawFCALBoundaries(); else DrawFCALBoundaries(color);
  if (color==-1) DrawMagnet(); else DrawMagnet(color, color);
}

