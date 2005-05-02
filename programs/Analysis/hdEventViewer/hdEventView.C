{
  gROOT->Reset();

    gROOT->LoadMacro("hdgeant.C");

  TCanvas *c1 = new TCanvas("c1","Viewer",10,10,800,800);
  c1->SetFillColor(0);
  c1->Divide(1,1);

  int colorPad = 41;
  int colorBCAL = 20;
  int colorCDC = 4;
  p1 = new TPad("p1","p1",0.05,0.02,0.95,0.82,colorPad,3,1);
  p1->Draw();
  p1->cd();


  float xmin[3] = {-100, -100, -100};
  float xmax[3] = {100, 100, 100};   
  TView *view = new TView(xmin, xmax, 1);
  //  TView *view = new TView(2);
  //view->SetParallel();
    hdgeant();

    gGeoManager->SetVisLevel(3);
    gGeoManager->SetVisLevel(3);
    //gGeoManager->GetVolume("BCAL")->SetLineColor(colorBCAL);
    gGeoManager->GetVolume("CDC")->SetLineColor(colorCDC);
    gGeoManager->GetMasterVolume()->Draw();
    view->Side();

    TGeoVolume *CDC = gGeoManager->GetVolume("CDC");  
    gGeoManager->SetTopVolume(CDC);  
    gGeoManager->SetVisLevel(1);
    gGeoManager->GetTopVolume()->Draw();


    




}
