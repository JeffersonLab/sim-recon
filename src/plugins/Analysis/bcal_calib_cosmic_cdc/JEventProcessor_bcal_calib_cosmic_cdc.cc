// $Id$
//
//    File: JEventProcessor_bcal_calib_cosmic_cdc.cc
// Created: Tue Jul  1 13:11:51 EDT 2014
// Creator: dalton (on Linux gluon05.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

#include "JEventProcessor_bcal_calib_cosmic_cdc.h"
using namespace jana;

#include <TRACKING/DTrackCandidate.h>
#include <DAQ/Df250PulseIntegral.h>
#include <DAQ/Df250PulseIntegral.h>
#include <BCAL/DBCALHit.h>
#include "BCAL/DBCALGeometry.h"

float getx(float r, float phi) {
	return r*cos(phi);
}
float gety(float r, float phi) {
	return r*sin(phi);
}
float getr(float x, float y) {
	return sqrt(x*x + y*y);
}
// float getphi(float x, float y) {
// 	return atan(y/x);
// }
float getposphi(float x, float y) {
	float phi=0;
	if (x!=0) phi = atan2(y,x);
	if (phi<0) phi+=2*PI;
	return phi;
}
float getdistance(float x1, float y1, float x2, float y2) {
	float deltax2 = (x1-x2)*(x1-x2);
	float deltay2 = (y1-y2)*(y1-y2);
	return sqrt(deltax2+deltay2);
}


// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
	void InitPlugin(JApplication *app){
		InitJANAPlugin(app);
		app->AddProcessor(new JEventProcessor_bcal_calib_cosmic_cdc());
	}
} // "C"


//------------------
// JEventProcessor_bcal_calib_cosmic_cdc (Constructor)
//------------------
JEventProcessor_bcal_calib_cosmic_cdc::JEventProcessor_bcal_calib_cosmic_cdc()
{

}

//------------------
// ~JEventProcessor_bcal_calib_cosmic_cdc (Destructor)
//------------------
JEventProcessor_bcal_calib_cosmic_cdc::~JEventProcessor_bcal_calib_cosmic_cdc()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_bcal_calib_cosmic_cdc::init(void)
{
	// This is called once at program startup. If you are creating
	// and filling historgrams in this plugin, you should lock the
	// ROOT mutex like this:
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();
	//

	/// Setup the parameters
	/**
	  VERBOSE=1  Output every event
	  VERBOSE=2  Output every intersection
	  VERBOSE=3  Output every sector
	  VERBOSE=4  Output every hit
	*/
	VERBOSE=0; 
	gPARMS->SetDefaultParameter("BCCC:VERBOSE",VERBOSE);
	
	/// Create the root tree
	bcal_calib_cosmic_cdc_tree = new TTree("bcal_calib_cosmic_cdc",
										   "tree of DBCALHit energies and the length of the track through each cell.");
	bcal_calib_cosmic_cdc_tree->Branch("eventnum",&eventnum,"eventnum/i");
	bcal_calib_cosmic_cdc_tree->Branch("cell",&cell,"cell/i");
//	bcal_calib_cosmic_cdc_tree->Branch("",&,"/i");
	bcal_calib_cosmic_cdc_tree->Branch("layer",&tlayer,"layer/i");
	bcal_calib_cosmic_cdc_tree->Branch("module",&tmodule,"module/i");
	bcal_calib_cosmic_cdc_tree->Branch("sector",&tsector,"sector/i");
	bcal_calib_cosmic_cdc_tree->Branch("globalsect",&tglobalsect,"globalsect/i");
	bcal_calib_cosmic_cdc_tree->Branch("numcells",&numcells,"numcells/i");
	bcal_calib_cosmic_cdc_tree->Branch("dist",&tdist,"dist/f");
	bcal_calib_cosmic_cdc_tree->Branch("use",&use,"use/f");
	bcal_calib_cosmic_cdc_tree->Branch("dse",&dse,"dse/f");
	bcal_calib_cosmic_cdc_tree->Branch("track_m",&track_m,"track_m/f");
	bcal_calib_cosmic_cdc_tree->Branch("track_c",&track_c,"track_c/f");
	bcal_calib_cosmic_cdc_tree->Branch("chisq",&chisq,"chisq/f");
	bcal_calib_cosmic_cdc_tree->Branch("Ndof",&Ndof,"Ndof/i");

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_bcal_calib_cosmic_cdc::brun(JEventLoop *eventLoop, int runnumber)
{
	// This is called whenever the run number changes
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_bcal_calib_cosmic_cdc::evnt(JEventLoop *loop, int eventnumber)
{
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// loop->Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.
	// Here's an example:
	//
	// vector<const MyDataClass*> mydataclasses;
	// loop->Get(mydataclasses);
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();


	/// DTrackCandidate:CDCCOSMIC
	/// Get a vector of DTrackCandidate:CDCCOSMIC objects for this event 
	vector<const DTrackCandidate*> CDCCOSMIC_vec;
	loop->Get(CDCCOSMIC_vec,"CDCCOSMIC");
	unsigned int num_CDCCOSMIC = CDCCOSMIC_vec.size();
	/// Create vectors to store the ID of the hit cells and the traversed distance
	vector<int> CellId_vec;
	vector<float> distance_vec;
	/// Only do anything if there is a track found
	if (num_CDCCOSMIC>0) {
		/// Get the parameters of the fitted track
		const DTrackCandidate* CDCCOSMIC = CDCCOSMIC_vec[0];
		const DVector3 cdcmom = CDCCOSMIC->momentum();
		const DVector3 cdcpos = CDCCOSMIC->position();
		chisq = CDCCOSMIC->chisq;
		Ndof = CDCCOSMIC->Ndof;

		/// Extract the slope and intercept
		track_m = cdcmom.y()/cdcmom.x();
		track_c = cdcpos.y() - track_m*cdcpos.x();
		if (VERBOSE>=1)  
			printf("BCCC >> pos (%6.2f,%6.2f,%6.2f)   mom (%6.2f,%6.2f,%6.2f)   m=%6.2f,  c=%6.2f, chisq=%6.2f, Ndof=%i\n",
				   cdcpos.x(), cdcpos.y(), cdcpos.z(), cdcmom.x(), cdcmom.y(), cdcmom.z(),track_m,track_c,chisq,Ndof);
		/// Store the parameters for both track intersections with the 5 layers
		float r[5], phi[2][5], x[2][5], y[2][5];
		/// For each layer boundary, calculate the intersection of the DTrackCandidate
		/// with the circle and store the r and phi value for both intersections.
		for (int laybound=0; laybound<=4; laybound++) {

			r[laybound] =  DBCALGeometry::fADC_radius[laybound];
			
			float A = 1 + track_m*track_m;
			float B = 2 * track_m * track_c;
			float C = track_c*track_c - r[laybound]*r[laybound];
			float disc = (B*B - 4*A*C);
			
			x[0][laybound] = 0, y[0][laybound] = 0, x[1][laybound] = 0, y[1][laybound] = 0;
			phi[0][laybound]=0 , phi[1][laybound]=0;
			if (disc > 0) {
				x[0][laybound] = (-B + sqrt(disc))/(2*A);
				x[1][laybound] = (-B - sqrt(disc))/(2*A);
				y[0][laybound] = track_m*x[0][laybound] + track_c;
				y[1][laybound] = track_m*x[1][laybound] + track_c;
				phi[0][laybound] = getposphi(x[0][laybound],y[0][laybound]);
				phi[1][laybound] = getposphi(x[1][laybound],y[1][laybound]);
			}
			if (VERBOSE>=2)  {
				for (int track=0; track<=1; track++) {
					int fADC_cellId = DBCALGeometry::fADCcellId_rphi( r[laybound], phi[track][laybound]);
					int module = DBCALGeometry::module( fADC_cellId );
					int layer = DBCALGeometry::layer( fADC_cellId );
					int sector = DBCALGeometry::sector( fADC_cellId );
					int glosector = DBCALGeometry::getglobalsector(module,sector);
					printf("BCCC >>  intersection: boundary=%i (x,y)=(%6.2f,%6.2f)  (r,phi)=(%6.2f,%6.3f)  cellId 0x%4x  (mod,lay,sec,glosec) = (%2i,%2i,%2i,%3i)\n",
						   laybound, x[track][laybound], y[track][laybound], r[laybound], phi[track][laybound], 
						   fADC_cellId, module, layer, sector, glosector);
				}
			}
			// fADC_cellId = DBCALGeometry::fADCcellId_rphi( r[laybound], phi[1][laybound]);
			// module = DBCALGeometry::module( fADC_cellId );
			// layer = DBCALGeometry::layer( fADC_cellId );
			// sector = DBCALGeometry::sector( fADC_cellId );
			// glosector = DBCALGeometry::getglobalsector(module,sector);
			// if (VERBOSE>=2)  {				
			// 	printf("BCCC >>  intersection: boundary=%i (x,y)=(%6.2f,%6.2f)  (r,phi)=(%6.2f,%6.3f)  cellId 0x%4x  (mod,lay,sec,glosec) = (%2i,%2i,%2i,%3i)\n",
			// 		   laybound, x[1][laybound], y[1][laybound], r[laybound], phi[1][laybound], fADC_cellId, module, layer, sector, glosector);
		}

		/// For each intersection set, step through each layer
		for (int track=0; track<=1; track++) {
			float total_distance = 0;
			/// For each layer, step through each sector
			for (int layer=0; layer<=3; layer++) {
				/// IN refers to inner radial boundary, OUT refers to outer radial boundary
				/// Find the global sector number for the entrance and exit sectors in this layer
				int fADC_cellIdIN = DBCALGeometry::fADCcellId_rphi( r[layer], phi[track][layer]);
				int moduleIN = DBCALGeometry::module( fADC_cellIdIN );
				int sectorIN = DBCALGeometry::sector( fADC_cellIdIN );
				int globalsectorIN = DBCALGeometry::getglobalsector(moduleIN,sectorIN);
				int fADC_cellIdOUT = DBCALGeometry::fADCcellId_rphi( r[layer+1], phi[track][layer+1]);
				int moduleOUT = DBCALGeometry::module( fADC_cellIdOUT );
				int sectorOUT = DBCALGeometry::sector( fADC_cellIdOUT );
				int globalsectorOUT = DBCALGeometry::getglobalsector(moduleOUT,sectorOUT);
				int globalsectormin, globalsectormax; // These are the ordered edge sectors, smallest and largets
				/// Check that the cells are valid
				if (fADC_cellIdIN<=0 || fADC_cellIdOUT<=0) {
					printf("BCCC >>Cells are not valid, event %i\n",eventnumber);
				} else {
					if (globalsectorIN == globalsectorOUT) {
						/// If the track only goes through 1 sector in this layer, then you are done, 
						/// get the distance from the (x,y) values of the layer inner and outer boundary intersections
						float dist = getdistance(x[track][layer], y[track][layer], x[track][layer+1], y[track][layer+1]);
						if (VERBOSE>=3)  
							printf("BCCC >>   single sector layer:   point (%6.3f,%6.3f)  point (%6.3f,%6.3f)\n",
								   x[track][layer], y[track][layer], x[track][layer+1], y[track][layer+1]);
						total_distance+=dist;
						if (VERBOSE>=2)  
							printf("BCCC >>  (mod,lay,sec,glosec) = (%2i,%2i,%2i,%3i)  fADC_cellId 0x%4x  philess %6.3f, phimore %6.3f   dist = %6.2f\n",
								   moduleIN,layer+1,sectorIN,globalsectorIN,fADC_cellIdIN,0.0,0.0,dist);
						CellId_vec.push_back(fADC_cellIdIN);
						distance_vec.push_back(dist);
					} else {
						/// If the track goes through multiple sectors then find the intersection of the DTrackCandidate
						/// with the sector boundary.
						if (globalsectorIN > globalsectorOUT) {
							globalsectormin = globalsectorOUT;
							globalsectormax = globalsectorIN;
						} else  {
							globalsectormax = globalsectorOUT;
							globalsectormin = globalsectorIN;
						}
						if (VERBOSE>=2) {
							printf("BCCC >>  layer=%i  phi=%5.2f  fADC_cellId 0x%4x  (mod,sec) = (%2i,%2i)  %3i\n", 
								   layer, phi[track][layer], fADC_cellIdIN, moduleIN, sectorIN, globalsectorIN);
							printf("BCCC >>  layer=%i  phi=%5.2f  fADC_cellId 0x%4x  (mod,sec) = (%2i,%2i)  %3i\n", 
								   layer, phi[track][layer+1], fADC_cellIdOUT, moduleOUT, sectorOUT, globalsectorOUT);
						}

						for (int glosect=globalsectormin; glosect<=globalsectormax; glosect++) {
							float dist = 0;
							int sector = DBCALGeometry::getsector(glosect);
							int module = DBCALGeometry::getmodule(glosect);
							int fADCId = DBCALGeometry::cellId(module,layer+1,sector); // layers are 1 to 4
							float phi = DBCALGeometry::phi(fADCId);
							float phihalfSize = DBCALGeometry::phiSize(fADCId)/2.;
							float philess=phi-phihalfSize, phimore=phi+phihalfSize;
							float mless, xless, yless, mmore, xmore, ymore;
							/// For each sector, get line boundary on each side and find the (x,y) point
							/// at which the track crosses that line
							mless = tan(philess);
							mmore = tan(phimore);
							xless = track_c/(mless-track_m);
							yless = (mless*track_c)/(mless-track_m);
							xmore = track_c/(mmore-track_m);
							ymore = (mmore*track_c)/(mmore-track_m);
							// printf("BCCC >>mless %6.2f, xless %6.2f, yless %6.2f, mmore %6.2f, xmore %6.2f, ymore %6.2f\n",
							// 	   mless, xless, yless, mmore, xmore, ymore);
							/// Edge sectors in each layer have 
							if (glosect==globalsectorIN) {
								if (glosect==globalsectormin) {
									dist = getdistance(x[track][layer],y[track][layer],xmore,ymore);
								} else {
									dist = getdistance(x[track][layer],y[track][layer],xless,yless);
								}
							} else {
								if (glosect==globalsectorOUT) {
									if (glosect==globalsectormin) {
										dist = getdistance(x[track][layer+1],y[track][layer+1],xmore,ymore);
									} else {
										dist = getdistance(x[track][layer+1],y[track][layer+1],xless,yless);
									}
								} else {
									/// For a central sector, get line boundary on each side and find the distance
									/// the phi less side and the phi more side							
									// printf("BCCC >>philess %6.2f, mless %6.2f, xless %6.2f, yless %6.2f, phimore %6.2f, mmore %6.2f, xmore %6.2f, ymore %6.2f\n",
									// 	   philess, mless, xless, yless, phimore, mmore, xmore, ymore);
									dist = getdistance(xless, yless, xmore, ymore);
								}
							}
							CellId_vec.push_back(fADCId);
							distance_vec.push_back(dist);
							total_distance+=dist;
							if (VERBOSE>=2) 
								printf("BCCC >>  (mod,lay,sec,glosec) = (%2i,%2i,%2i,%3i)  fADC_cellId 0x%4x  philess %6.3f, phimore %6.3f   dist = %6.2f\n",
									   module,layer+1,sector,glosect,fADCId,philess,phimore,dist);
						} 			
					} // end multiple sector logic
				} // end check that cells are valid
			} // end loop over layers
			if (VERBOSE>=1) {
				float dist = getdistance(x[track][0],y[track][0],x[track][4],y[track][4]);
				printf("BCCC >> total distance in BCAL: %f cm, expected %f\n",total_distance,dist);
			}
		}
	}


	/// DBCALHit
	/// Get a vector of DBCALHit objects for this event (1 object for each crate/slot/channel above threshold)
	vector<const DBCALHit*> BCALHit_vec;
	loop->Get(BCALHit_vec);
	unsigned int num_BCALHit = BCALHit_vec.size();
	/// Creat maps from the BCAL cell ID to the DBCALHit object pointer
	
	map<int, const DBCALHit*>::iterator myiter;
	map<int, const DBCALHit*> upstream;
	map<int, const DBCALHit*> downstream;
	/// Loop over all DBCALHit objects in this event
	for(unsigned int c_chan=0; c_chan<num_BCALHit; c_chan++){
		const DBCALHit *BCALHit = BCALHit_vec[c_chan];
		//int fADCId = DBCALGeometry::fADCId_fADC( BCALHit->module, BCALHit->layer, BCALHit->sector );
		int fADCId = DBCALGeometry::cellId( BCALHit->module, BCALHit->layer, BCALHit->sector );
		if (VERBOSE>=4) 
			printf("BCCC >>    Hit: (module, layer, sector) = (%2i,%2i,%2i)  fADCId 0x%x\n",\
				   BCALHit->module, BCALHit->layer, BCALHit->sector, fADCId);
		/// Fill the maps with the hits for this event.
		if (BCALHit->end == DBCALGeometry::kUpstream ) {
			upstream[fADCId] = BCALHit;
		} else  {
			if (BCALHit->end == DBCALGeometry::kDownstream ) {
				downstream[fADCId] = BCALHit;
			} else {
				printf("BCCC >>Bad BCAL enum\n");
			}
		}
	}
	if (VERBOSE>=1) 
		printf("BCCC >> found %ld upstream and %ld downstream hits for map.\n",upstream.size(),downstream.size());
	//cout << "found " << upstream.size() << " upstream and " << downstream.size() << " downstream hits for map.\n";

	/// Trees are filled with data
	japp->RootWriteLock();
	
	eventnum = eventnumber;
	/// Loop over the intersected cells
	unsigned int num_CellId = CellId_vec.size();
	for (unsigned int cellnum=0; cellnum<num_CellId; cellnum++) {
		cell = CellId_vec[cellnum];
		tdist = distance_vec[cellnum];
		tmodule = DBCALGeometry::module( cell );
		tlayer = DBCALGeometry::layer( cell );
		tsector = DBCALGeometry::sector( cell );
		tglobalsect = DBCALGeometry::getglobalsector(tmodule,tsector);
		numcells = num_CellId;
		use=0; 
		dse=0;
		myiter = upstream.find(cell);
		if (myiter!=upstream.end()) {
			const DBCALHit* US = myiter->second;
			use = US->E;
		} 
		myiter = downstream.find(cell);
		if (myiter!=downstream.end()) {
			const DBCALHit* DS = myiter->second;
			dse = DS->E;
		} 
		if (VERBOSE>=2) 
			printf("BCCC >>  eventnum %4i  cellnum %2i  0x%4x  (mod,lay,sec) = (%2i,%2i,%2i)  %3i   %6.2f %6.2f %6.2f \n",
				   eventnum, cellnum, cell,  tmodule, tlayer, tsector, tglobalsect, tdist, use, dse);


		bcal_calib_cosmic_cdc_tree->Fill();
	}


	japp->RootUnLock();

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_bcal_calib_cosmic_cdc::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_bcal_calib_cosmic_cdc::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}

 
/* emacs
 * Local Variables:
 * mode:C++
 * mode:font-lock
 * c-file-style: "stroustrup"
 * tab-width: 4
 * End:
 */
