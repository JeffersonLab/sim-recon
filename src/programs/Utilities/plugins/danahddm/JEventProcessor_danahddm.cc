// JEventProcessor_danahddm.cc
//
//
// JANA event processor plugin writes out hddm event to file
//
//
//  David Lawrence, 7-May-2010

#include <JANA/JApplication.h>
#include <HDDM/DEventSourceHDDM.h>
#include <TRACKING/DTrackTimeBased.h>

#include "JEventProcessor_danahddm.h"


// hddm output file name, use hddm:FILENAME configuration parameter to override
static string hddmFileName = "dana_events.hddm";


// mutex for serializing writing to file
static pthread_mutex_t hddmMutex = PTHREAD_MUTEX_INITIALIZER;

// Make us a plugin
// for initializing plugins
extern "C" {
  void InitPlugin(JApplication *app) {
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_danahddm(), true);
  }
} // "extern C"

//-------------------------------
// Constructor
//-------------------------------
JEventProcessor_danahddm::JEventProcessor_danahddm()
{

  jout << endl << "  Default JEventProcessor_danahddm invoked" << endl << endl;

  // Check for hddm:FILENAME output file name parameter
  gPARMS->SetDefaultParameter("hddm:FILENAME",hddmFileName);
  jout << endl << "  hddm output file name is " << hddmFileName << endl << endl;
  
  file = NULL;
  Nevents_written = 0;
}  

//-------------------------------
// Destructor
//-------------------------------
JEventProcessor_danahddm::~JEventProcessor_danahddm()
{
	
}

//-------------------------------
// init
//-------------------------------
jerror_t JEventProcessor_danahddm::init(void)
{
	return NOERROR;
}

//-------------------------------
// brun
//-------------------------------
jerror_t JEventProcessor_danahddm::brun(JEventLoop *loop, int runnumber)
{
	// If file is already open, don't reopen it. Just keep adding to it.
	if(file)return NOERROR;

	// We wait until here to open the output so that we can check if the 
	// input is hddm. If it's not, tell the user and exit immediately
	JEvent& event = loop->GetJEvent();
	JEventSource *source = event.GetJEventSource();
	DEventSourceHDDM *hddm_source = dynamic_cast<DEventSourceHDDM*>(source);
	if(!hddm_source){
		cerr<<" This program MUST be used with an HDDM file as input!"<<endl;
		exit(-1);
	}

	// If we got here, it must be an HDDM source. Open a new file.
	file = init_s_HDDM((char*)hddmFileName.c_str());
	Nevents_written = 0;

	return NOERROR;
}

//-------------------------------
// evnt
//-------------------------------
jerror_t JEventProcessor_danahddm::evnt(JEventLoop *loop, int eventnumber)
{
	// This is a little complicated. We need to get a hold of the s_HDDM_t
	// structure pointer for this event so we can pass it to flush_s_HDDM()
	// along with our ouput stream pointer. The flush routine frees up the
	// memory in the s_HDDM_t structure. When the framework tries "flush"ing
	// a second time, we get a seg fault. To prevent the framework from
	// flushing, we have to clear the free_on_flush flag (by default set
	// to true). This means we need to get the DEventSource pointer and
	// downcast to a DEventSourceHDDM structure. It's a little strange setting
	// this for every event, but we have no way of knowing when the event
	// source changes and this at least guarantees it for all event sources.
	JEvent& event = loop->GetJEvent();
	JEventSource *source = event.GetJEventSource();
	DEventSourceHDDM *hddm_source = dynamic_cast<DEventSourceHDDM*>(source);
	if(!hddm_source){
		cerr<<" This program MUST be used only with HDDM files as inputs!"<<endl;
		exit(-1);
	}
	s_HDDM_t *hddm = (s_HDDM_t*)event.GetRef();
	if(!hddm)return NOERROR;
	
	// Delete any data in the reconView branch of the event. This may be a little
	// confusing. We want to delete the reconView branch and any memory allocated
	// to its sub-braches. The easiest way to do this is to create another
	// event that contains a pointer ONLY to the reconView branch and then
	// delete that event. The original event will need to have it's reconView
	// point set to HDDM_NULL to indicate it's empty. Note that one would
	// expect for most cases the input file not to have any reconView data.
	s_ReconView_t *recon = (s_ReconView_t*)HDDM_NULL;
	s_PhysicsEvents_t* PE = hddm->physicsEvents;
	if(PE){
		for(unsigned int i=0; i<PE->mult; i++){
			s_ReconView_t *my_recon = PE->in[i].reconView;
			if(my_recon != HDDM_NULL){
				// Create a new, temporary event
				s_HDDM_t *tmp_hddm = make_s_HDDM();
				s_PhysicsEvents_t *tmp_PE = make_s_PhysicsEvents(1);
				
				// Move recon branch over to temporary tree
				tmp_PE->mult=1;
				tmp_PE->in[0].reconView = my_recon;
				tmp_PE->in[0].reactions = (s_Reactions_t*)HDDM_NULL;
				tmp_PE->in[0].hitView = (s_HitView_t*)HDDM_NULL;
				PE->in[i].reconView = (s_ReconView_t*)HDDM_NULL;
				
				// Delete temporary tree and any memory in the recon branch along with it
				flush_s_HDDM(tmp_hddm, NULL);
			}
			
			// Create a new reconView branch to hang our data from (only for
			// first physics event).
			if(recon==HDDM_NULL)recon = PE->in[i].reconView = make_s_ReconView();
		}
	}
	
	// In order to do anything worthwhile here, we need to have a valid recon
	// pointer.
	if(recon==NULL || recon==HDDM_NULL)return NOERROR;

	// Fill in reconstructed banks, replacing any that are already there
	Add_DTrackTimeBased(loop, recon);

	// get write lock
	pthread_mutex_lock(&hddmMutex);

	// Write event to file and update counter
	flush_s_HDDM(hddm, file);
	Nevents_written++;

	// unlock
	pthread_mutex_unlock(&hddmMutex);

	// Tell the JEventSourceHDDM object not to free this event a second time
	hddm_source->flush_on_free = false;

	return NOERROR;
}

//-------------------------------
// erun
//-------------------------------
jerror_t JEventProcessor_danahddm::erun(void)
{
	return NOERROR;
}

//-------------------------------
// fini
//-------------------------------
jerror_t JEventProcessor_danahddm::fini(void)
{
	if(file){
		close_s_HDDM(file);
		cout<<endl<<"Closed HDDM file"<<endl;
	}
	cout<<" "<<Nevents_written<<" event written to "<<hddmFileName<<endl;

	return NOERROR;
}

//-------------------------------
// Add_DTrackTimeBased
//-------------------------------
void JEventProcessor_danahddm::Add_DTrackTimeBased(JEventLoop *loop, s_ReconView_t *recon)
{
	// Get objects to write out
	vector<const DTrackTimeBased*> tracktimebaseds;
	loop->Get(tracktimebaseds);
	if(tracktimebaseds.size()==0)return;

	// Allocate memory for all time based tracks
	s_Tracktimebaseds_t *tbt = recon->tracktimebaseds = make_s_Tracktimebaseds(tracktimebaseds.size());
	tbt->mult = 0;

	for(unsigned int i=0; i<tracktimebaseds.size(); i++, tbt->mult++){
		const DTrackTimeBased *tbt_dana = tracktimebaseds[i];
		s_Tracktimebased_t *tbt_hddm = &(tbt->in[tbt->mult]);
		
		DVector3 pos = tbt_dana->position();
		DVector3 mom = tbt_dana->momentum();
		
		tbt_hddm->FOM = tbt_dana->FOM;
		tbt_hddm->candidateid = tbt_dana->candidateid;
		tbt_hddm->trackid = tbt_dana->trackid;
		tbt_hddm->id = tbt_dana->id;
		tbt_hddm->chisq = tbt_dana->chisq;
		tbt_hddm->Ndof = tbt_dana->Ndof;

		tbt_hddm->momentum = make_s_Momentum();
		tbt_hddm->properties = make_s_Properties();
		tbt_hddm->origin = make_s_Origin();
		tbt_hddm->errorMatrix = make_s_ErrorMatrix();
		tbt_hddm->TrackingErrorMatrix = make_s_TrackingErrorMatrix();

		tbt_hddm->momentum->E = tbt_dana->energy();
		tbt_hddm->momentum->px = mom.x();
		tbt_hddm->momentum->py = mom.y();
		tbt_hddm->momentum->pz = mom.z();

		tbt_hddm->properties->charge = tbt_dana->charge();
		tbt_hddm->properties->mass = tbt_dana->mass();
		
		tbt_hddm->origin->t = 0.0;
		tbt_hddm->origin->vx = pos.x();
		tbt_hddm->origin->vy = pos.y();
		tbt_hddm->origin->vz = pos.z();
		
		string vals = DMatrixDSymToString(tbt_dana->errorMatrix());
		tbt_hddm->errorMatrix->Ncols = 7;
		tbt_hddm->errorMatrix->Nrows = 7;
		tbt_hddm->errorMatrix->type = strdup("DMatrixDSym"); // HDDM always frees strings automatically
		tbt_hddm->errorMatrix->vals = strdup(vals.c_str()); // HDDM always frees strings automatically

		tbt_hddm->TrackingErrorMatrix->Ncols = 5;
		tbt_hddm->TrackingErrorMatrix->Nrows = 5;
		tbt_hddm->TrackingErrorMatrix->type = strdup("DMatrixDSym"); // HDDM always frees strings automatically
		tbt_hddm->TrackingErrorMatrix->vals = strdup(DMatrixDSymToString(tbt_dana->TrackingErrorMatrix()).c_str());

	}
}

//-------------------------------
// DMatrixDSymToString
//-------------------------------
string JEventProcessor_danahddm::DMatrixDSymToString(const DMatrixDSym &mat)
{
	// Convert the given symmetric matrix into a single string that
	// can be used in an HDDM file.

	stringstream ss;
	for(int irow=0; irow<mat.GetNrows(); irow++) {
		for(int icol=irow; icol<mat.GetNcols(); icol++) {
			ss << mat[irow][icol] << " ";
		}
	}
	
	return ss.str();
}
