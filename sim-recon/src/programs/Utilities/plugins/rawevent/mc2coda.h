/* mc2coda.h - Header File for mc2coda library */



/* defines */
#define OK    0
#define ERROR -1

#define MAX_PARAM 10

/* Module Types (0-15) and modes */
/* (Updated 6/24/2013 by DL to be consistent with document: */
/* "VME Data Format Standards for JLAB Modules"  */
enum type_id_t {
	TID,               // =0
	FADC250,           // =1
	FADC125,           // =2
	F1TDC32,           // =3
	F1TDC48,           // =4
	JLAB_TS,           // =5
	TD,                // =6
	SSP,               // =7
	JLAB_DISC,         // =8
	MODULE_TYPE_RES1,  // =9
	MODULE_TYPE_RES2,  // =10
	MODULE_TYPE_RES3,  // =11
	MODULE_TYPE_RES4,  // =12
	MODULE_TYPE_RES5,  // =13
	MODULE_TYPE_RES6,  // =14
	MODULE_TYPE_RES7,  // =15
	
	// The following are not defined by the DAQ group (i.e.
	// they don't control the data format so can't encode
	// the module type in it)
	UNKNOWN,           // =16
	VMECPU,            // =17
	CAEN1190,          // =18
	CAEN1290,          // =19
	
	N_MODULE_TYPES     // Make sure this is the last thing in the enum!
};
//#define NOMOD    0
//#define VMECPU   1
//#define TID      2
//#define FADC250  3
#define FADC250_MODE_RW  1   /* Raw Window */
#define FADC250_MODE_WI  2   /* Window Integrated */
#define FADC250_MODE_RP  3   /* Raw Pulse */
#define FADC250_MODE_IP  4   /* Integrated Pulse  */
//#define FADC125  4
#define FADC125_MODE_RW  1
#define FADC125_MODE_WI  2
#define FADC125_MODE_RP  3
#define FADC125_MODE_IP  4

//#define F1TDC32  5
//#define F1TDC48  6
//#define CAEN1190 7
#define CAEN1190_MODE_TM 1    /* Trigger Matching mode */
//#define CAEN1290 8
#define CAEN1290_MODE_TM 1    /* Trigger Matching mode */
//#define JLDISC   10
//#define USERMOD  14


/* Define limits for CODA DAQ Systems */
#define MAX_EXP        10
#define MAX_CRATES    128
#define MAX_SLOTS      21
#define MIN_SLOT        2
#define MAX_SLOT       21

#define MIN_TRIG_TIME  2048    /* minimum 8 microseconds for trigger time */

#define MAX_EVENT_SIZE     1048576    /* 1 Megabyte Default Event buffer size */
#define MAX_HIT_SIZE       256        /* Max Bytes for 1 Hit - eg RAW window mode for FADC250 */
#define MAX_TOTAL_HITS     50000      /* Define maximum number of hits that can be stored */
#define MAX_HITS_PER_SLOT  500        /* Max hits stored per slot (module) */
#define MAX_HITS_PER_CHAN  50         /* Max hits allowed per channel - most modules will less than this */

#define MC2CINIT_NULL      0
#define MC2CINIT_DEFINED   1
#define MC2CINIT_COMPLETE  2


typedef struct coda_crate_map {
	int crate_id;
	int nModules;
	unsigned int moduleMask;
	int module_map[MAX_SLOTS];
	int det_map[MAX_SLOTS];
} CODA_CRATE_MAP;

typedef struct coda_exp_info {
	char expname[64];
	int ncrates;
	int inited;
	int openevents;
	unsigned short rocid[MAX_CRATES];
	CODA_CRATE_MAP *crate[MAX_CRATES];
} CODA_EXP_INFO;


typedef struct coda_event_info {
	uint64_t eventid;
	uint64_t trigtime;
	uint32_t evtype;
	struct coda_exp_info *expid;
	int nhits;
	int hcount[MAX_CRATES][MAX_SLOTS];
	struct coda_hit_info *hits[MAX_CRATES][MAX_SLOTS];
	int maxBytes;
	unsigned int *evbuf;
} CODA_EVENT_INFO;


typedef struct coda_hit_info {
	int hit_id;                  /* Unique identifier for each "hit" */
	int det_id;                  /* Detector ID   0<=val<=4095       */
	int crate_id;                /* Crate ID      0<=val<=127        */
	int slot_id;                 /* Slot ID       2<=val<=21         */
	int chan_id;                 /* Channel number - dependent on module type */
	int module_id;               /* User defined module type         */
	int module_mode;             /* Module Operational Mode          */
	int module_param[MAX_PARAM]; /* Module specific program parameters - times - windows latencies etc...*/
	uint32_t nwords;             /* Number of words in the hdata array */
	uint32_t *hdata;             /* Hit data - dependent on module type and mode */
} CODA_HIT_INFO;



/* Macros */

#define ROC_BANK_OPEN(status, id, nevents) {			    \
StartOfRocBank = dabufp; \
*(++dabufp) = (((status) << 28) | (id) << 16) | 0x1000 | (nevents);\
(dabufp)++;

#define ROC_BANK_CLOSE \
*StartOfRocBank = (uint32_t) (dabufp - StartOfRocBank -1);	\
} \


#define DATA_BANK_OPEN(status, detid, nevents) {			    \
uint32_t *StartOfBank; \
StartOfBank = dabufp; \
*(++dabufp) = (((status) << 28) | (detid) << 16) | 0x0100 | (nevents);\
(dabufp)++;

#define DATA_BANK_CLOSE \
*StartOfBank = (uint32_t) (dabufp - StartOfBank - 1);	\
} \






/* Prototypes */
void mc2codaSetRunNumber(unsigned int run_number);
CODA_EXP_INFO *mc2codaInitExp(int nCrates, const char *name);
int mc2codaSetCrate(CODA_EXP_INFO *expID, int crateid, int nmod, int *modules, int *detid);
CODA_EVENT_INFO *mc2codaOpenEvent(CODA_EXP_INFO *expID, uint64_t eventNum, uint64_t trigTime, unsigned short eventType, int maxSize);
int mc2codaWrite(CODA_EVENT_INFO *eventID, int nHits, struct coda_hit_info *codaHits);
void mc2codaStats(CODA_EVENT_INFO *eventID, int sflag);
unsigned int mc2codaCloseEvent(CODA_EVENT_INFO *eventID);
int mc2codaResetEvent(CODA_EVENT_INFO *eventID, uint64_t eventNum, uint64_t trigTime, unsigned short eventType);
int mc2codaFreeEvent(CODA_EVENT_INFO *eventID);
void mc2codaFree(CODA_EXP_INFO *expID);

