#include<TGClient.h>
#include<TCanvas.h>
#include<TF1.h>
#include<TRandom.h>
#include<TGButton.h>
#include<TGFrame.h>
#include<TRootEmbeddedCanvas.h>
#include<RQ_OBJECT.h>

class MyMainFrame{

  RQ_OBJECT("MyMainFrame")

  private: 
    TGMainFrame     *fMain;
    TRootEmbeddedCanvas *fEcanvas;
    TPad *p1;
    TView *fView;
    Int_t fThetaView = -60;
    Int_t fPhiView = 0;
    Int_t fPsiView = 90;
    Double_t fDview = 100;

    int colorPad = 41;
    int colorBCAL = 20;
    int colorCDC = 4;

  public:
    MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h);
    TView *View() {return fView;};
    virtual ~MyMainFrame();
    void DrawSITE();
    void DrawCDC();
    void DrawBCAL();
    void DrawFTOF();
    void DrawFCAL();
    void DrawCERE();
    void DrawFDC();
    void DrawLASS();
    void DrawBOTH();
    void ThetaPlus();
    void PhiPlus();
    void ThetaMinus();
    void PhiMinus();
    void SideView();
    void SetSidewaysView();
    void UpdateView();
};

MyMainFrame::MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h)
{
  //Create a main frame
  fMain = new TGMainFrame(p, w, h);

  TGHorizontalFrame *mFrame = new TGHorizontalFrame(fMain,600,600);

  TGVerticalFrame *vframe1 = new TGVerticalFrame(mFrame, 40, 40);
  TGVerticalFrame *vframe2 = new TGVerticalFrame(mFrame, 40, 40);


  TGTextButton *exit2 = new TGTextButton(vframe1, "&Exit2", "gApplication->Terminate(0)");
  vframe1->AddFrame(exit2, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  //Create canvas widget
  fEcanvas = new TRootEmbeddedCanvas("Ecanvas", mFrame, 600, 600);
  p1 = new TPad("p1","p1",0.05,0.02,0.95,0.82,colorPad,3,1);
  p1->Draw();
  p1->cd();


  //Create a horizontal frame widget with buttons

  TGTextButton *drawSITE = new TGTextButton(vframe2, "&SITE");
  drawSITE->Connect("Clicked()", "MyMainFrame", this, "DrawSITE()");
  vframe2->AddFrame(drawSITE, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *drawCDC = new TGTextButton(vframe2, "&CDC");
  drawCDC->Connect("Clicked()", "MyMainFrame", this, "DrawCDC()");
  vframe2->AddFrame(drawCDC, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *drawBCAL = new TGTextButton(vframe2, "&BCAL");
  drawBCAL->Connect("Clicked()", "MyMainFrame", this, "DrawBCAL()");
  vframe2->AddFrame(drawBCAL, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *drawFTOF = new TGTextButton(vframe2, "&FTOF");
  drawFTOF->Connect("Clicked()", "MyMainFrame", this, "DrawFTOF()");
  vframe2->AddFrame(drawFTOF, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *drawFCAL = new TGTextButton(vframe2, "&FCAL");
  drawFCAL->Connect("Clicked()", "MyMainFrame", this, "DrawFCAL()");
  vframe2->AddFrame(drawFCAL, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *drawCERE = new TGTextButton(vframe2, "&CERE");
  drawCERE->Connect("Clicked()", "MyMainFrame", this, "DrawCERE()");
  vframe2->AddFrame(drawCERE, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *drawFDC = new TGTextButton(vframe2, "&FDC");
  drawFDC->Connect("Clicked()", "MyMainFrame", this, "DrawFDC()");
  vframe2->AddFrame(drawFDC, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *drawLASS = new TGTextButton(vframe2, "&LASS");
  drawLASS->Connect("Clicked()", "MyMainFrame", this, "DrawLASS()");
  vframe2->AddFrame(drawLASS, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *drawBOTH = new TGTextButton(vframe2, "&BOTH");
  drawBOTH->Connect("Clicked()", "MyMainFrame", this, "DrawBOTH()");
  vframe2->AddFrame(drawBOTH, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *phiM = new TGTextButton(vframe2, "&phi-");
  phiM->Connect("Clicked()", "MyMainFrame", this, "PhiMinus()");
  vframe2->AddFrame(phiM, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *phiP = new TGTextButton(vframe2, "&phi+");
  phiP->Connect("Clicked()", "MyMainFrame", this, "PhiPlus()");
  vframe2->AddFrame(phiP, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *thetaM = new TGTextButton(vframe2, "&theta-");
  thetaM->Connect("Clicked()", "MyMainFrame", this, "ThetaMinus()");
  vframe2->AddFrame(thetaM, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *thetaP = new TGTextButton(vframe2, "&theta+");
  thetaP->Connect("Clicked()", "MyMainFrame", this, "ThetaPlus()");
  vframe2->AddFrame(thetaP, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *sideV = new TGTextButton(vframe2, "&Side View");
  sideV->Connect("Clicked()", "MyMainFrame", this, "SideView()");
  vframe2->AddFrame(sideV, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *setSidewaysView = new TGTextButton(vframe2, "&Set Side View");
  setSidewaysView->Connect("Clicked()", "MyMainFrame", this, "SetSidewaysView()");
  vframe2->AddFrame(setSidewaysView, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *exit = new TGTextButton(vframe2, "&Exit", "gApplication->Terminate(0)");
  vframe2->AddFrame(exit, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  mFrame->AddFrame(vframe1, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));
  mFrame->AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 1));
  mFrame->AddFrame(vframe2, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));

  fMain->AddFrame(mFrame, new TGLayoutHints(kLHintsCenterX|kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 1));

  cerr << "Here.....0" << endl;

  // Set a name to the main frame
  fMain->SetWindowName("Simple Example");

  // Map all subwindows of main frame
  fMain->MapSubwindows();


  // Initialize the layout algorithm
  fMain->Resize(fMain->GetDefaultSize());

  // Map main frame
  fMain->MapWindow();

  //gROOT->LoadMacro("hdgeant.C");
  gROOT->LoadMacro("hddsroot.C");
  float xmin[3] = {-100, -100, -100};
  float xmax[3] = {100, 100, 100};
  fView = new TView(xmin, xmax, 1);

  //hdgeant();
  hddsroot();

  this->DrawSITE();


}

void MyMainFrame::DrawSITE()
{
  // Draws function graphics in randomly choosen interval
  gGeoManager->SetVisLevel(3);
  gGeoManager->GetMasterVolume()->Draw();

  UpdateView();

}

void MyMainFrame::DrawBCAL()
{
  // Draws function graphics in randomly choosen interval
  TGeoVolume *BCAL = gGeoManager->GetVolume("BCAL");
  gGeoManager->SetTopVolume(BCAL);
  gGeoManager->SetVisLevel(1);
  gGeoManager->GetTopVolume()->Draw();

  UpdateView();
}

void MyMainFrame::DrawFTOF()
{
  // Draws function graphics in randomly choosen interval
  TGeoVolume *FTOF = gGeoManager->GetVolume("FTOF");
  gGeoManager->SetTopVolume(FTOF);
  gGeoManager->SetVisLevel(3);
  gGeoManager->GetTopVolume()->Draw();

  UpdateView();
}

void MyMainFrame::DrawFCAL()
{
  // Draws function graphics in randomly choosen interval
  TGeoVolume *FCAL = gGeoManager->GetVolume("FCAL");
  gGeoManager->SetTopVolume(FCAL);
  gGeoManager->SetVisLevel(2);
  gGeoManager->GetTopVolume()->Draw();


  UpdateView();
}

void MyMainFrame::DrawCERE()
{
  // Draws function graphics in randomly choosen interval
  TGeoVolume *CERE = gGeoManager->GetVolume("CERE");
  gGeoManager->SetTopVolume(CERE);
  gGeoManager->SetVisLevel(2);
  gGeoManager->GetTopVolume()->Draw();


  UpdateView();
}

void MyMainFrame::DrawFDC()
{
  // Draws function graphics in randomly choosen interval
  TGeoVolume *FDC = gGeoManager->GetVolume("FDC");
  gGeoManager->SetTopVolume(FDC);
  gGeoManager->SetVisLevel(2);
  gGeoManager->GetTopVolume()->Draw();


  UpdateView();
}

void MyMainFrame::DrawLASS()
{
  // Draws function graphics in randomly choosen interval
  TGeoVolume *LASS = gGeoManager->GetVolume("LASS");
  gGeoManager->SetTopVolume(LASS);
  gGeoManager->SetVisLevel(2);
  gGeoManager->GetTopVolume()->Draw();


  UpdateView();
}

void MyMainFrame::DrawCDC()
{
  // Draws function graphics in randomly choosen interval
  TGeoVolume *CDC = gGeoManager->GetVolume("CDC");
  gGeoManager->SetTopVolume(CDC);
  gGeoManager->SetVisLevel(1);
  gGeoManager->GetTopVolume()->Draw();


  UpdateView();
}

void MyMainFrame::DrawBOTH()
{
  // Draws function graphics in randomly choosen interval
  TGeoMedium *medium = 0;
  TGeoVolume *vol = gGeoManager->MakeBox("TOP",medium,100,250,250);
  gGeoManager->SetTopVolume(vol);

  TGeoVolume *BCAL = gGeoManager->GetVolume("BCAL");
  TGeoVolume *CDC = gGeoManager->GetVolume("CDC");
  TGeoVolume *FTOF = gGeoManager->GetVolume("FTOF");
  TGeoVolume *CERE = gGeoManager->GetVolume("CERE");

  vol->AddNode(CDC, 1, gGeoIdentity);
  vol->AddNode(BCAL, 1, gGeoIdentity);
  vol->AddNode(FTOF, 1, gGeoIdentity);
  vol->AddNode(CERE, 1, gGeoIdentity);

  gGeoManager->SetTopVolume(vol);
  gGeoManager->SetVisLevel(1);
  //  gGeoManager->CloseGeometry();
  gGeoManager->GetTopVolume()->Draw();

   UpdateView();
}

void MyMainFrame::PhiPlus()
{
  fPhiView+=10;
  UpdateView();
}

void MyMainFrame::ThetaPlus()
{
  fThetaView+=10;
  UpdateView();
}

void MyMainFrame::PhiMinus()
{
  fPhiView-=10;
  UpdateView();
}

void MyMainFrame::ThetaMinus()
{
  fThetaView-=10;
  UpdateView();
}

void MyMainFrame::UpdateView()
{
  const Float_t kPI = Float_t (TMath::Pi());

  Int_t irep;

  Float_t longitude_deg = fPhiView;// * 180.0/kPI - 90;
  Float_t  latitude_deg = fThetaView;//  * 180.0/kPI + 90;
  Float_t       psi_deg = fPsiView ;//      * 180.0/kPI;

  p1->GetView()->SetView(longitude_deg, latitude_deg, 90, irep);

//  p1->SetPhi(-90 - longitude_deg);
//  p1->SetTheta(90 - latitude_deg);

  p1->Modified(kTRUE);
  p1->Update();

}

void MyMainFrame::SetSidewaysView()
{
  const Float_t kPI = Float_t (TMath::Pi());

  Int_t irep;

  Float_t longitude_deg = fPhiView;// * 180.0/kPI - 90;
  Float_t  latitude_deg = fThetaView;//  * 180.0/kPI + 90;
  Float_t       psi_deg = fPsiView ;//      * 180.0/kPI;

  p1->GetView()->SetView(longitude_deg, latitude_deg, 90, irep);

  p1->SetPhi(-90 - longitude_deg);
  p1->SetTheta(90 - latitude_deg);

  p1->Modified(kTRUE);
  p1->Update();

  cerr << "Presing SidewaysView" << endl;
}

void MyMainFrame::SideView()
{
  fView->Side();
  this->fPhiView = 0;
  this->fThetaView = 90;
  this->fPsiView = 0;

}

MyMainFrame::~MyMainFrame()
{
  // Clean up used widgets: frames, buttons, layouthints
  fMain->Cleanup();
  delete fMain;
}

void hdEventView()
{
  // Popup the GUI.....
  MyMainFrame *mFrame = new MyMainFrame(gClient->GetRoot(), 800, 800);
}
