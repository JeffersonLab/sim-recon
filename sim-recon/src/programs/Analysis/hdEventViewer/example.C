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

  public:
    MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h);
    virtual ~MyMainFrame();
    void DoDraw();
};

MyMainFrame::MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h)
{
  //Create a main frame
  fMain = new TGMainFrame(p, w, h);

  //Create canvas widget
  fEcanvas = new TRootEmbeddedCanvas("Ecanvas", fMain, 200, 200);
  fMain->AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 1));

  //Create a horizontal frame widget with buttons
  TGHorizontalFrame *hframe = new TGHorizontalFrame(fMain, 200, 40);

  TGTextButton *draw = new TGTextButton(hframe, "&Draw");
  draw->Connect("Clicked()", "MyMainFrame", this, "DoDraw()");
  hframe->AddFrame(draw, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

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

}

void MyMainFrame::DoDraw()
{
  // Draws function graphics in randomly choosen interval
  TF1 *f1 = new TF1("f1","sin(x)/x", 0, gRandom->Rndm()*10);
  f1->SetFillColor(19);
  f1->SetFillStyle(1);
  f1->SetLineWidth(3);
  f1->Draw();
  TCanvas *fCanvas = fEcanvas->GetCanvas();
  fCanvas->cd();
  fCanvas->Update();
}

MyMainFrame::~MyMainFrame()
{
  // Clean up used widgets: frames, buttons, layouthints
  fMain->Cleanup();
  delete fMain;
}

void example()
{
  // Popup the GUI.....
  new MyMainFrame(gClient->GetRoot(), 200, 200);
}
