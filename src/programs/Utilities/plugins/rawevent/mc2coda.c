/* mc2coda - function to accept a list of data and
 convert it into a CODA3 Raw Event format
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdint.h>

#include "mc2coda.h"


/* Global mc2coda Library Variables */
int mc2coda_inited = MC2CINIT_NULL;
CODA_EXP_INFO mc2coda_expid;
int mc2coda_maxevsize = 0;


/* Local variables */
static int mc2coda_ncrates_defined = 0;
static unsigned int *dabufp, *StartOfRocBank;
static unsigned int RUN_NUMBER = 1;


/* Include module specific definitions */
#include "mc2coda_modules.h"


/* Allow user to set the run number (added 8/21/2013 DL)*/
void
mc2codaSetRunNumber(unsigned int run_number)
{
	RUN_NUMBER = run_number;
}

/* Initialize the Crate/Module maps for crates in the system
 Must be called prior to calling the mc2coda event packer.
 
 - Returns an expID (experiment setup structure pointer) */

CODA_EXP_INFO *
mc2codaInitExp(int nCrates, const char *name)
{
	
	int ii;
	
	if(nCrates>MAX_CRATES) {
		printf("mc2codaInitExp: ERROR: Too many crates defined for the system\n");
		return (NULL);
	}
	
	/* Initialize a new experiment config structure */
	if(mc2coda_inited) {
		printf("mc2codaInitExp: ERROR: Experiment already Initialized - Call mc2codaFree(expid) first\n");
		return(NULL);
	}
	
	printf("Exp name = %s\n",name);
	
	bzero((char *) &mc2coda_expid,sizeof(CODA_EXP_INFO));
	
	strncpy(mc2coda_expid.expname,name,64);
	
	mc2coda_expid.ncrates = nCrates;
	
	for(ii=0;ii<nCrates;ii++) {
		mc2coda_expid.rocid[ii] = (ii+1);
		mc2coda_expid.crate[ii] = (CODA_CRATE_MAP *) malloc(sizeof(CODA_CRATE_MAP));
		bzero(mc2coda_expid.crate[ii],sizeof(CODA_CRATE_MAP));
		mc2coda_expid.crate[ii]->crate_id = ii;
	}
	
	mc2coda_expid.inited = MC2CINIT_DEFINED;       /* Only partially initialized */
	mc2coda_expid.openevents = 0;   /* No open events using this experiment config */
	mc2coda_ncrates_defined = 0;
	printf("Initialized: Experiment %s (id = 0x%lx) consisting of %d crates/ROCs\n",
		   mc2coda_expid.expname, (unsigned long)&mc2coda_expid,mc2coda_expid.ncrates);
	
	return(&mc2coda_expid);
	
}


/* Set Crate specific info for a given expID -
 All crates must be defined - with module types before an experiment config
 can generate RAW events. This function must be called at least once for each
 crate defined by mc2codaInitExp
 
 pass number of modules, and module and detector id arrays
 
 returns Number of crates initialized  or Error (-1)    */

int
mc2codaSetCrate(CODA_EXP_INFO *expID, int crateid, int nmod, int *modules, int *detid)
{
	
	int ii;
	
	/* check that parameters are valid */
	
	if ((crateid<0)||(crateid > expID->ncrates)) {
		printf("mc2codaSetCrate: ERROR: crate ID out of range %d > %d (ncrates) \n",
			   crateid, expID->ncrates);
		return(-1);
	}
	if(nmod >= MAX_SLOTS) {
		printf("mc2codaSetCrate: ERROR: Number of modules out of range %d\n",nmod);
		return(-1);
	}
	
	if(mc2coda_ncrates_defined >= expID->ncrates) {
		mc2coda_ncrates_defined = expID->ncrates - 1;
		printf("mc2codaSetCrate: WARN: overwriting Crate %d module map\n",crateid);
	}
	
	/* Fill the Crate structure */
	printf("mc2codSetCrate: INFO: Crate %d has %d modules for readout\n",crateid, nmod);
	expID->crate[(crateid-1)]->nModules = nmod;
	expID->crate[(crateid-1)]->moduleMask = 0;
	for (ii=0;ii<MAX_SLOTS;ii++) {
		expID->crate[(crateid-1)]->module_map[ii] = modules[ii];
		expID->crate[(crateid-1)]->det_map[ii] = detid[ii];
		if(modules[ii] != 0) expID->crate[(crateid-1)]->moduleMask |= (1<<ii);
	}
	
	mc2coda_ncrates_defined++;
	
	/* Initialization is Complete ?? */
	if(mc2coda_ncrates_defined == expID->ncrates) {
		printf("mc2codaSetCrate: All %d crates defined: Initialization complete.\n",mc2coda_ncrates_defined);
		mc2coda_inited = MC2CINIT_COMPLETE;
	}
	
	return(mc2coda_ncrates_defined);
	
}



/* Open a RAW event for population
 
 - Returns an eventID */
CODA_EVENT_INFO *
mc2codaOpenEvent(CODA_EXP_INFO *expID, uint64_t eventNum, uint64_t trigTime, unsigned short eventType, int maxSize)
{
	
	int ii, jj, has_tt = 0;
	unsigned int *evbuf;
	CODA_EVENT_INFO *evinfo;
	
	if(expID == 0) {
		printf("m2codaOpenEvent: ERROR: expID not defined\n");
		return(NULL);
	}else{  /*check if the event number is Valid */
		if(eventNum == 0) {
			printf("mc2codaOpenEvent: ERROR: event Number invalid (0)\n");
			return(NULL);
		}
		if(trigTime == 0)
			has_tt=0;
		else
			has_tt=1;
	}
	
	/* Allocate an Event Info structure */
	evinfo = (CODA_EVENT_INFO *) malloc(sizeof(CODA_EVENT_INFO));
	evinfo->nhits   = 0;
	evinfo->eventid = eventNum;
	evinfo->trigtime = trigTime;
	evinfo->evtype = eventType&0x0000ffff;
	evinfo->expid   = expID;
	
	/* Allocate Hit Arrays for each valid crate/slot */
	for (ii=0; ii<(expID->ncrates); ii++) {
		for(jj=1; jj<(MAX_SLOTS); jj++) {  /* Skip CPU slot */
			
			evinfo->hcount[ii][jj] = 0;
			/* Check if there is a valid module for this crate/slot pair  and
			 allocate the MAX number of HIT structures for that slot */
			if((expID->crate[ii]->moduleMask)&(1<<jj)) {
				/* printf("DEBUG: Allocating hit array for crate,slot = %d,%d\n",ii,jj); */
				evinfo->hits[ii][jj] = (CODA_HIT_INFO *) malloc(MAX_HITS_PER_SLOT * sizeof(CODA_HIT_INFO));
			}else{
				evinfo->hits[ii][jj] = NULL;
			}
			
		}
	}
	
	/* Allocate an Event buffer */
	if(maxSize) {
		evbuf = (unsigned int *) malloc (maxSize);
		mc2coda_maxevsize = maxSize;
	}else{
		evbuf = (unsigned int *) malloc (MAX_EVENT_SIZE);
		mc2coda_maxevsize = MAX_EVENT_SIZE;
	}
	
	bzero((char *)evbuf, mc2coda_maxevsize);
	evinfo->evbuf = (unsigned int *)evbuf;
	evinfo->maxBytes = mc2coda_maxevsize;
	
	
	/* Populate the Event header and trigger Bank */
	if(has_tt) {
		/* Event header if we have a trigger time */
		evbuf[0]  =         12;
		evbuf[1]  = 0xff511001;
		evbuf[2]  =         10;
		evbuf[3]  = 0xff232000 | ((expID->ncrates)&0xff); /* changed from ff21 to include run number 8/21/2013 DL */
		evbuf[4]  = 0x010a0006;  /* segment of 64 bit uints */
		memcpy((char *)&evbuf[5],(char *)&eventNum,8);
		memcpy((char *)&evbuf[7],(char *)&trigTime,8);
		evbuf[ 9] = 0x01;       /* run type */
		evbuf[10] = RUN_NUMBER; /* This goes into high 32 bits which seems backwards ?? */
		evbuf[11] = 0x01050001; /* segment of shorts header with 1 value */
		evbuf[12] = (eventType);
	}else{
		/* Event header if we don't have a trigger time */
		evbuf[0]  =         10;
		evbuf[1]  = 0xff501001;
		evbuf[2]  =          8;
		evbuf[3]  = 0xff222000 | ((expID->ncrates)&0xff); /* changed from ff20 to include run number 8/21/2013 DL */
		evbuf[4]  = 0x010a0004; /* segment of 64 bit uints */
		memcpy((char *)&evbuf[5],(char *)&eventNum,8);
		evbuf[7]  = 0x01;       /* run type */
		evbuf[8]  = RUN_NUMBER; /* This goes into high 32 bits which seems backwards ?? */
		evbuf[9]  = 0x01050001; /* segment of shorts header with 1 value */
		evbuf[10] = (eventType);
	}
	
	
	/* Bump the openevents counter */
	expID->openevents++;
	
	
	return(evinfo);
}

/* Write Monte Carlo hit(s) info into the event.
 
 This routine allocates local hit structures and copies hit infomation
 into them and then updates a list. Once the function returns the original
 hit structures can be freed or cleared.
 
 This routine must be reentrant so that different threads can write hits
 asynchonously. No attempt is made to reorder hits. Only copy them and
 updates list totals.
 
 Returns: # of hits written to the Event
 
 */

int
mc2codaWrite(CODA_EVENT_INFO *event, int nHits, struct coda_hit_info *codaHits)
{
	
	int ii, cnt, lcnt, crate, slot;
	CODA_EXP_INFO *exp;
	CODA_HIT_INFO *tmpH;
	static int bad_crate_slot_warning_issued = 0;
	
	/* check for valid hits */
	if((nHits <= 0) || (codaHits == NULL)) {
		return(0);
	}
	
	/* Check for a valid event pointer */
	if(event) {
		exp = event->expid;
		if(exp->openevents <= 0) {
			printf("mc2codaWrite: ERROR: Invalid or Corrupted event ID\n");
			return(-1);
		}
	}else{
		printf("mc2codaWrite: ERROR: Null event ID\n");
		return(-1);
	}
	
	
	/* Now loop through all new hits and sort/save them*/
	lcnt = 0;
	for (ii=0; ii<nHits; ii++) {
		crate = codaHits[ii].crate_id - 1;
		slot  = codaHits[ii].slot_id - 1;
		/* printf("DEBUG: Writing hit %d for crate,slot = %d,%d\n",cnt,crate,slot); */
		if(crate<0 || slot<0){
			if(!bad_crate_slot_warning_issued){
				printf("mc2codaWrite: ERROR: Invalid crate(%d) or slot(%d)! This could be due to\n", crate, slot);
				printf("                     data for a channel not defined in the crate map.\n");
				printf("                     This is only reported once.\n");
				bad_crate_slot_warning_issued = 1;
			}
		}else{
		
			cnt   = event->hcount[crate][slot];
			if(cnt > MAX_HITS_PER_SLOT) {
				printf("mc2codaWrite: ERROR: No available space to store hit %d for crate/slot = %d/%d\n",
					   codaHits->hit_id,crate,slot);
			} else {
				tmpH = (CODA_HIT_INFO *)&event->hits[crate][slot][cnt];
				if(tmpH == NULL){
					printf("%s:%d ERROR!! no CODA_HIT_INFO structure allocated for crate=%d, slot=%d\n", __FILE__, __LINE__, crate+1, slot+1);
					continue;
				}

				memcpy((char *)&(tmpH->hit_id),(char *)&(codaHits[ii].hit_id),sizeof(CODA_HIT_INFO)) ;
				/* printf("DEBUG: malloc data array\n"); */
				tmpH->hdata = (uint32_t *) malloc((codaHits[ii].nwords)<<2);
				memcpy((char *)(tmpH->hdata), (char *)(codaHits[ii].hdata),(codaHits[ii].nwords)<<2);
				
				event->hcount[crate][slot] += 1;
			}
		} // crate<0 || slot<0
		
		event->nhits++;
		lcnt++;
		
	}
	
	return(lcnt);
}



/* Close out the event packer and return pointer to Event Buffer for writing
 function also returns number of total words in the event (eventID[0] */
unsigned int
mc2codaCloseEvent(CODA_EVENT_INFO *event)
{
	int ii, jj;
	unsigned int roc, det, mod;
	unsigned int *StartofEvent;
	CODA_EXP_INFO *expID = event->expid;
	CODA_CRATE_MAP *crate;
	
	if(event == NULL) {
		printf("mc2codaCloseEvent: ERROR: Null Event Pointer!!\n");
		return(0);
	}
	
	/* Setup pointers */
	StartofEvent = &(event->evbuf[0]);
	dabufp = StartofEvent + (event->evbuf[0] + 1);  /* Set data pointer to end of the event */
	
	
	/* Loop through all ROCs  */
	for(ii=0; ii<expID->ncrates; ii++) {
		
		roc = expID->rocid[ii];
		crate = expID->crate[ii];

		// This needs to be fixed:
		// Each ROC may put out one or more DATA BLOCK banks(DBB). If the ROC
		// contains only JLab modules, it may concatentate the data from
		// all of them into a single DBB. It may also split things up such
		// that multiple DMA transfers are used, each going to a different
		// DBB, but from the same ROC.
		//
		// Another possibility is that there are modules that are not JLab
		// modules so the data itself cannot be interpreted to determine
		// the exact module type. In this case, the detid bits (lower 8bits
		// of "tag") must be used to specify the module type. This means that
		// these modules must go in a separate DBB.
		//
		// To handle this, the DBB should be opened using the module type
		// (0 for non-digitizing, 1 for JLab, 20 for CAEN1290, ....) of
		// the first module in the crate. If a module with a different
		// det_id is encountered, the DBB should be closed and another
		// opened using that new det_id. This means we need to find the
		// first non-empty slot, not counting the first slot which is
		// the CPU.
		det = 0;
		for(jj=0; jj<MAX_SLOTS; jj++){
			if(((1<<jj)&(crate->moduleMask)) == 0) continue;
			det = crate->det_map[jj];
			if(det != 0) break;
		}

		// If no digitization modules exist in this crate then skip it
		if(det == 0) continue;

		/* Open a ROC data Bank */
		ROC_BANK_OPEN(0,roc,1);
				
		/* Open a Data Bank */
		// Note that we can't use the DATA_BANK_OPEN and DATA_BANK_CLOSE macros
		// since they open and close a curly bracket{} and we want to close
		// an old bank and open a new one inside the for loop below.
		//DATA_BANK_OPEN(0,det,1);
		uint32_t *StartOfBank = dabufp;
		uint32_t status = 0;
		uint32_t nevents = 1;
		*(++dabufp) = (((status) << 28) | (det) << 16) | 0x0100 | (nevents);
		(dabufp)++;
		
		/* Loop through all modules */
		for(jj=0; jj<MAX_SLOTS; jj++) {
			if(((1<<jj)&(crate->moduleMask)) == 0) continue;
			
			// If module type changes, then we need to open a new
			// data bank
			if(det != crate->det_map[jj]){
				// Close existing bank
				*StartOfBank = (uint32_t) (dabufp - StartOfBank - 1);
				
				// Open new bank with the new detid (but not for detid==0)
				det = crate->det_map[jj];
				StartOfBank = dabufp;
				if(det != 0){
					*(++dabufp) = (((status) << 28) | (det) << 16) | 0x0100 | (nevents);
					(dabufp)++;
				}
			}
			
			mod = crate->module_map[jj];
			/*      printf(" slot %d  module %d \n",(jj+1),mod); */
			
			
			switch (mod) {
				case FADC250:
					fadc250_write_data (event, roc,(jj+1) , FADC250_MODE_IP);
					break;
				case FADC125:
					fadc125_write_data (event, roc,(jj+1) , FADC125_MODE_IP);
					break;
				case F1TDC32:
					f1tdc32_write_data (event, roc,(jj+1) ,0);
					break;
				case F1TDC48:
					f1tdc48_write_data (event, roc,(jj+1), 0);
					break;
				case CAEN1290:
					caen1290_write_data (event, roc,(jj+1), 0);
					break;
				default:
					break;
			}
			
		}
		
		// Close data bank (see note above regarding DATA_BANK_CLOSE)
		//DATA_BANK_CLOSE;
		*StartOfBank = (uint32_t) (dabufp - StartOfBank - 1);
	}
	
	ROC_BANK_CLOSE;
	*StartofEvent = (uint32_t) (dabufp - StartofEvent - 1);
	
	return((*StartofEvent+1));
}


/* Reset the Existing event for use again. This is equivilent to OpenEvent
 so that new event Number, trigger time and event type must be specified
 but it avoids mallocing the necessary memory and structures all over again.
 
 returns 0 - OK  or -1 Error
 
 */
int
mc2codaResetEvent(CODA_EVENT_INFO *eventID, uint64_t eventNum, uint64_t trigTime, unsigned short eventType)
{
	
	CODA_EXP_INFO  *exp;
	CODA_HIT_INFO  *tmpH;
	int ii, jj, kk, ccnt, hcnt, lcnt;
	
	
	if(eventID != NULL) {
		exp = eventID->expid;
		if(exp->openevents <= 0) {
			printf("mc2codaResetEvent: ERROR: Invalid or Corrupted event ID\n");
			return(-1);
		}
		ccnt = exp->ncrates;
		hcnt = eventID->nhits;
	} else {
		printf("mc2codaResetEvent: ERROR: NULL pointer to Event Buffer\n");
		return(-1);
	}
	
	/* First free the data arrays that have been allocated for individual hits */
	for (ii=0; ii<ccnt; ii++) {
		for(jj=1; jj<(MAX_SLOTS); jj++) {  /* Skip CPU slot */
			if(eventID->hits[ii][jj] != NULL) {
				/* crate,slot array exists, now check if it has hits */
				lcnt = eventID->hcount[ii][jj];
				if(lcnt) {
					/* free all the hit data */
					for(kk=0;kk<lcnt;kk++) {
						tmpH = (CODA_HIT_INFO *) &eventID->hits[ii][jj][kk];
						free(tmpH->hdata);
					}
					eventID->hcount[ii][jj] = 0;  /* Set Hit count to 0 */
				}
			}
		}
	}
	
	/* Now clear the existing event structure and buffer and reset with the new info */
	eventID->nhits   = 0;
	eventID->eventid = eventNum;
	eventID->trigtime = trigTime;
	eventID->evtype = eventType&0x0000ffff;
	
	if(eventID->maxBytes > 0) {
		bzero((char *)eventID->evbuf, eventID->maxBytes);
	}else{
		printf("mc2codaResetEvent: ERROR: Event buffer size is invalid (%d)\n",eventID->maxBytes);
		return(-1);
	}
	
	
	/* Populate the Event header and trigger Bank */
	if(trigTime) {
		eventID->evbuf[0]  =         10;
		eventID->evbuf[1]  = 0xff511001;
		eventID->evbuf[2]  =          8;
		eventID->evbuf[3]  = 0xff212000 | ((exp->ncrates)&0xff);
		eventID->evbuf[4]  = 0x010a0004;
		memcpy((char *)&eventID->evbuf[5],(char *)&eventNum,8);
		memcpy((char *)&eventID->evbuf[7],(char *)&trigTime,8);
		eventID->evbuf[9]  = 0x01050001;
		eventID->evbuf[10] = (eventType);
	}else{
		eventID->evbuf[0]  =          8;
		eventID->evbuf[1]  = 0xff501001;
		eventID->evbuf[2]  =          6;
		eventID->evbuf[3]  = 0xff202000 | ((exp->ncrates)&0xff);
		eventID->evbuf[4]  = 0x010a0004;
		memcpy((char *)&eventID->evbuf[5],(char *)&eventNum,8);
		eventID->evbuf[7]  = 0x01050001;
		eventID->evbuf[8] = (eventType);
	}
	
	
	return(0);
}



/* Free the event buffer for the specified event */
int
mc2codaFreeEvent(CODA_EVENT_INFO *eventID)
{
	
	CODA_EXP_INFO  *exp;
	CODA_HIT_INFO  *tmpH;
	int ii, jj, kk, ccnt, hcnt, lcnt;
	
	/* Get the crate and Hit counts */
	exp = eventID->expid;
	ccnt = exp->ncrates;
	hcnt = eventID->nhits;
	
	if(eventID != NULL) {
		/* First free the Event Buffer */
		if(eventID->evbuf != NULL) {
			free(eventID->evbuf);
			eventID->evbuf = NULL;
		} else {
			printf("mc2codaFreeEvent: ERROR invalid pointer to Event Buffer\n");
		}
		
		/* Free all the allocated Arrays of Hit structures and their associated data array*/
		for (ii=0; ii<ccnt; ii++) {
			for(jj=1; jj<(MAX_SLOTS); jj++) {  /* Skip CPU slot */
				
				if(eventID->hits[ii][jj] != NULL) {
					/* crate,slot array exists, now check if it has hits */
					lcnt = eventID->hcount[ii][jj];
					if(lcnt) {
						/* free all the hit data */
						for(kk=0;kk<lcnt;kk++) {
							/* printf("DEBUG: freeing data\n"); */
							tmpH = (CODA_HIT_INFO *) &eventID->hits[ii][jj][kk];
							free(tmpH->hdata);
						}
					}
					/* Now free the hit structure array */
					/* printf("DEBUG: freeing Hit array for crate,slot = %d,%d\n",ii,jj); */
					free(eventID->hits[ii][jj]);
					eventID->hits[ii][jj] = NULL;
				}
			}
		}
		
		/* Free the Event Structure */
		free(eventID);
		eventID = NULL;
		
		/* Decrement the openevents counter */
		exp->openevents--;
		
	}else{
		printf("mc2codaFreeEvent: ERROR invalid eventID \n");
		return(ERROR);
	}
	
	return(0);
}


/* Free the Experiment Config (including any open events for it) */
void
mc2codaFree(CODA_EXP_INFO *expID)
{
	int ii;
	
	if(expID == NULL) {
		printf("mc2codaFree: ERROR invalid expID\n");
		return;
	}
	
	/* Check if there are open events */
	if(expID->openevents) {
		printf("mc2codaFree: ERROR: Cannot Free experiment. There are currently %d open events\n",
			   expID->openevents);
		return;
	}
	
	if(expID->inited > MC2CINIT_NULL) {
		printf("mc2codaFree: WARN: Deleting Experiment %s\n",expID->expname);
		for(ii=0;ii<(expID->ncrates);ii++)
			free(expID->crate[ii]);
		
		bzero(expID,sizeof(CODA_EXP_INFO));
	}
}




/* Print out various statistics about the event */
void
mc2codaStats(CODA_EVENT_INFO *eventID, int sflag)
{
	
	CODA_EXP_INFO  *exp;
	//CODA_HIT_INFO  *tmpH;
	int ii, jj, ccnt, hcnt;
	int chits [MAX_CRATES];
	
	if(eventID == NULL) {
		printf("mc2codaStats: ERROR: Null Event ID\n");
		return;
	}
	
	/* Get the crate and Hit counts */
	exp = eventID->expid;
	ccnt = exp->ncrates;
	hcnt = eventID->nhits;
	
	
	/* Print out Experiment Info */
	printf("Experiment: %s   # of Open Events: %d\n",exp->expname, exp->openevents);
	printf("  Total Crates: %d    Total Hits: %d\n",ccnt, hcnt);
	
	/* Count up Hits per Crate and print out Crate map with hit info*/
	printf("Crate Map:  moduleID( hit_count)  for each slot in each crate\n");
	for(ii=0;ii<ccnt;ii++) {
		chits[ii]=0;
		printf("  Crate %3d: ",(ii+1));
		for(jj=0;jj<MAX_SLOTS;jj++) {
			chits[ii] += eventID->hcount[ii][jj];
			printf("%1d(%4d)  ",exp->crate[ii]->module_map[jj],eventID->hcount[ii][jj]);
		}
		printf("\n");
	}
	
	
	return;
}
