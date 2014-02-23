

#include <iostream>
#include <iomanip>
using namespace std;

#include <TFile.h>
#include <TTree.h>

#include "cern_c.h"
#define MEMH 8000000
#define LREC 1024		/* record length of hbook direct access file in WORDS */
#define LUN 3			/* logical unit number of hbook file */
extern "C" {
	float pawc_[MEMH];
	int quest_[100];
};


#define MAX_PARTS 50

class event_t{
	public:
		int n;
		float E[MAX_PARTS];
		float px[MAX_PARTS];
		float py[MAX_PARTS];
		float pz[MAX_PARTS];
		float neuETot;
		int nGen;
		float genE[MAX_PARTS];
		float genPx[MAX_PARTS];
		float genPy[MAX_PARTS];
		float genPz[MAX_PARTS];
		int type[MAX_PARTS];
		float genNeuE;
};

int main(int narg, char *argv[])
{
	const char *fname ="primex_event.root";
	if(narg>1)fname = argv[1];

	// Open ROOT file for input
	TFile *f = new TFile(fname);
	if(!f || !f->IsOpen()){
		cerr<<"Unable to open file \""<<fname<<"\"!"<<endl;
		return -1;
	}	
	cout<<"Opened input file \""<<fname<<"\""<<endl;
	
	// Find TTree
	TTree *tree = (TTree*)f->Get("event");
	if(!tree){
		cerr<<"Unable to find TTree \"event\" in ROOT file!"<<endl;
		return -2;
	}
	
	// Set up branches
	event_t event;
	tree->SetBranchAddress("n", &event.n);
	tree->SetBranchAddress("E", &event.E);
	tree->SetBranchAddress("px", &event.px);
	tree->SetBranchAddress("py", &event.py);
	tree->SetBranchAddress("pz", &event.pz);
	tree->SetBranchAddress("neuETot", &event.neuETot);
	tree->SetBranchAddress("nGen", &event.nGen);
	tree->SetBranchAddress("genE", &event.genE);
	tree->SetBranchAddress("genPx", &event.genPx);
	tree->SetBranchAddress("genPy", &event.genPy);
	tree->SetBranchAddress("genPz", &event.genPz);
	tree->SetBranchAddress("type", &event.type);
	tree->SetBranchAddress("genNeuE", &event.genNeuE);

	// Initialize cernlib (monkey shines!)
	quest_[9] = 65000;
	int memh = MEMH;
	hlimit(memh);
	
	// Open HBOOK file for writing
	hropen(LUN, "lun", "primex_event.hbook" , "N", LREC, 0);
	
	// Create Ntuple
	hbnt(10,"event","");
	char ntp_str[256];
	sprintf(ntp_str, "n[0,%d]:I,E(n):R,px(n),py(n),pz(n),neuETot,nGen[0,%d]:I,genE(nGen):R,genPx(nGen),genPy(nGen),genPz(nGen),type(nGen):I,genNeuE:R", MAX_PARTS, MAX_PARTS);
	cout<<ntp_str<<endl;
	hbname(10,"EVNT", &event.n, ntp_str);
	
	// Loop over entries in ROOT file
	int Nevents = tree->GetEntries();
	for(int i=1; i<=Nevents; i++){
		// Read in event
		tree->GetEntry(i);

		// Write event to Ntuple
		hfnt(10);
		
		// Update ticker
		if(i%100 == 0)cout<<"  "<<i<<" events processed     \r";
	}
	
	cout<<endl;
	cout<<Nevents<<" events processed total."<<endl;
	
	// Close hbook file
	int icycle=0;
	hrout(0,icycle,"T");
	hrend("lun");
	
	// Close ROOT file
	f->Close();

	return 0;
}


#if 0
******************************************************************************
*Tree    :event     : Event Reconstruction                                   *
*Entries :    10000 : Total =         1882594 bytes  File  Size =    1363534 *
*        :          : Tree compression factor =   1.37                       *
******************************************************************************
*Br    0 :n         : n/I                                                    *
*Entries :    10000 : Total  Size=      40677 bytes  File Size  =       4586 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   6.98     *
*............................................................................*
*Br    1 :E         : E[n]/F                                                 *
*Entries :    10000 : Total  Size=     228491 bytes  File Size  =     183329 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   1.19     *
*............................................................................*
*Br    2 :px        : px[n]/F                                                *
*Entries :    10000 : Total  Size=     228505 bytes  File Size  =     187776 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   1.16     *
*............................................................................*
*Br    3 :py        : py[n]/F                                                *
*Entries :    10000 : Total  Size=     228505 bytes  File Size  =     187842 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   1.16     *
*............................................................................*
*Br    4 :pz        : pz[n]/F                                                *
*Entries :    10000 : Total  Size=     228505 bytes  File Size  =     184061 *
*Baskets :        8 : Basket Size=      32000 bytes  Compression=   1.18     *
*............................................................................*
*Br    5 :neuETot   : neuETot/F                                              *
*Entries :    10000 : Total  Size=      40719 bytes  File Size  =      26283 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.22     *
*............................................................................*
*Br    6 :nGen      : nGen/I                                                 *
*Entries :    10000 : Total  Size=      40698 bytes  File Size  =        263 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression= 121.66     *
*............................................................................*
*Br    7 :genE      : genE[nGen]/F                                           *
*Entries :    10000 : Total  Size=     161217 bytes  File Size  =     118567 *
*Baskets :        6 : Basket Size=      32000 bytes  Compression=   1.30     *
*............................................................................*
*Br    8 :genPx     : genPx[nGen]/F                                          *
*Entries :    10000 : Total  Size=     161229 bytes  File Size  =     122566 *
*Baskets :        6 : Basket Size=      32000 bytes  Compression=   1.25     *
*............................................................................*
*Br    9 :genPy     : genPy[nGen]/F                                          *
*Entries :    10000 : Total  Size=     161229 bytes  File Size  =     122616 *
*Baskets :        6 : Basket Size=      32000 bytes  Compression=   1.25     *
*............................................................................*
*Br   10 :genPz     : genPz[nGen]/F                                          *
*Entries :    10000 : Total  Size=     161229 bytes  File Size  =     118466 *
*Baskets :        6 : Basket Size=      32000 bytes  Compression=   1.30     *
*............................................................................*
*Br   11 :type      : type[nGen]/I                                           *
*Entries :    10000 : Total  Size=     161210 bytes  File Size  =      14124 *
*Baskets :        6 : Basket Size=      32000 bytes  Compression=  10.88     *
*............................................................................*
*Br   12 :genNeuE   : genNeuE/F                                              *
*Entries :    10000 : Total  Size=      40719 bytes  File Size  =      20990 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.52     *
*............................................................................*
root [3] 
#endif


