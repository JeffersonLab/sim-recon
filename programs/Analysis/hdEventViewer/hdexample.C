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
    TView *fView;
    Int_t fThetaView = -60;
    Int_t fPhiView = 0;
    Int_t fPsiView = 0;

    int colorPad = 41;
    int colorBCAL = 20;
    int colorCDC = 4;

  public:
    MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h);
    virtual ~MyMainFrame();
    void DrawSITE();
    void DrawCDC();
    void DrawBCAL();
    void ThetaPlus();
    void PhiPlus();
    void ThetaMinus();
    void PhiMinus();
    void SideView();
};

MyMainFrame::MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h)
{
  //Create a main frame
  fMain = new TGMainFrame(p, w, h);

  //Create canvas widget
  fEcanvas = new TRootEmbeddedCanvas("Ecanvas", fMain, 600, 600);
  fMain->AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 1));

  //Create a horizontal frame widget with buttons
  TGHorizontalFrame *hframe = new TGHorizontalFrame(fMain, 200, 40);

  TGTextButton *drawSITE = new TGTextButton(hframe, "&SITE");
  drawSITE->Connect("Clicked()", "MyMainFrame", this, "DrawSITE()");
  hframe->AddFrame(drawSITE, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *drawCDC = new TGTextButton(hframe, "&CDC");
  drawCDC->Connect("Clicked()", "MyMainFrame", this, "DrawCDC()");
  hframe->AddFrame(drawCDC, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *drawBCAL = new TGTextButton(hframe, "&BCAL");
  drawBCAL->Connect("Clicked()", "MyMainFrame", this, "DrawBCAL()");
  hframe->AddFrame(drawBCAL, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *phiM = new TGTextButton(hframe, "&phi-");
  phiM->Connect("Clicked()", "MyMainFrame", this, "PhiMinus()");
  hframe->AddFrame(phiM, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *phiP = new TGTextButton(hframe, "&phi+");
  phiP->Connect("Clicked()", "MyMainFrame", this, "PhiPlus()");
  hframe->AddFrame(phiP, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *thetaM = new TGTextButton(hframe, "&theta-");
  thetaM->Connect("Clicked()", "MyMainFrame", this, "ThetaMinus()");
  hframe->AddFrame(thetaM, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *thetaP = new TGTextButton(hframe, "&theta+");
  thetaP->Connect("Clicked()", "MyMainFrame", this, "ThetaPlus()");
  hframe->AddFrame(thetaP, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *sideV = new TGTextButton(hframe, "&Side View");
  sideV->Connect("Clicked()", "MyMainFrame", this, "SideView()");
  hframe->AddFrame(sideV, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton *exit = new TGTextButton(hframe, "&Exit", "gApplication->Terminate(0)");
  hframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  fMain->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));

  cerr << "Here.....0" << endl;

  // Set a name to the main frame
  fMain->SetWindowName("Simple Example");

  // Map all subwindows of main frame
  fMain->MapSubwindows();


  // Initialize the layout algorithm
  fMain->Resize(fMain->GetDefaultSize());

  // Map main frame
  fMain->MapWindow();

  gROOT->LoadMacro("hdgeant.C");
  float xmin[3] = {-100, -100, -100};
  float xmax[3] = {100, 100, 100};
  fView = new TView(xmin, xmax, 1);

  hdgeant();


}

void MyMainFrame::DrawSITE()
{
  // Draws function graphics in randomly choosen interval
  gGeoManager->SetVisLevel(3);
  gGeoManager->GetMasterVolume()->Draw();

  fView->RotateView(fPhiView, fThetaView);
  fView->SetPsi(this->fPsiView);      
  fView->Zoom();    

  TCanvas *fCanvas = fEcanvas->GetCanvas();
  fCanvas->cd();
  fCanvas->Update();

}

void MyMainFrame::DrawBCAL()
{
  // Draws function graphics in randomly choosen interval
  TGeoVolume *BCAL = gGeoManager->GetVolume("BCAL");
  gGeoManager->SetTopVolume(BCAL);
  gGeoManager->SetVisLevel(1);
  gGeoManager->GetTopVolume()->Draw();

  fView->RotateView(fPhiView, fThetaView);
  fView->SetPsi(this->fPsiView);      
  fView->Zoom();    

  TCanvas *fCanvas = fEcanvas->GetCanvas();
  fCanvas->cd();
  fCanvas->Update();

}

void MyMainFrame::DrawCDC()
{
  // Draws function graphics in randomly choosen interval
  TGeoVolume *CDC = gGeoManager->GetVolume("CDC");
  gGeoManager->SetTopVolume(CDC);
  gGeoManager->SetVisLevel(1);
  gGeoManager->GetTopVolume()->Draw();

  fView->RotateView(fPhiView, fThetaView);
  fView->SetPsi(this->fPsiView);      
  fView->Zoom();    

  TCanvas *fCanvas = fEcanvas->GetCanvas();
  fCanvas->cd();
  fCanvas->Update();

}

void MyMainFrame::PhiPlus()
{
  fPhiView+=10;
  fView->RotateView(fPhiView, fThetaView);
  fView->SetPsi(this->fPsiView);      
  fView->Zoom();    
}

void MyMainFrame::ThetaPlus()
{
  fThetaView+=10;
  fView->RotateView(fPhiView, fThetaView);
  fView->SetPsi(this->fPsiView);      
  fView->Zoom();    
}

void MyMainFrame::PhiMinus()
{
  fPhiView-=10;
  fView->RotateView(fPhiView, fThetaView);
  fView->SetPsi(this->fPsiView);      
  fView->Zoom();    
}

void MyMainFrame::ThetaMinus()
{
  fThetaView-=10;
  fView->RotateView(fPhiView, fThetaView);
  fView->SetPsi(this->fPsiView);      
  fView->Zoom();    
}

void MyMainFrame::SideView()
{
  //  fView->Side();
  //  this->fPhiView = 0;
  //  this->fThetaView = 90;

  this->fPhiView = 0;
  this->fThetaView = 90;
  this->fPsiView = 90;
  fView->RotateView(this->fPhiView,this->fThetaView);
  fView->SetPsi(this->fPsiView);      
  fView->Zoom();    
}

MyMainFrame::~MyMainFrame()
{
  // Clean up used widgets: frames, buttons, layouthints
  fMain->Cleanup();
  delete fMain;
}

void hdexample()
{
  // Popup the GUI.....
  new MyMainFrame(gClient->GetRoot(), 600, 600);
}
