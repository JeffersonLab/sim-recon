/*
 * disData.h
 *
*/

#ifndef disDataH_INCLUDED
#define disDataH_INCLUDED

static const char sccsid_disDataH[] = "@(#)disData.h\t5.7\tCreated 12/11/97 15:27:54, \tcompiled "__DATE__;

#include <rectypes.h>

/********************************************

An event will look like this:

---------------------------------------------
event      header              itape_header_t
           ----------------------------------
           group1              group_header_t
                               data
           ----------------------------------
           group2              group_header_t
                               data
           ----------------------------------
	   ...
           ----------------------------------
           groupN              group_header_t
                               data
---------------------------------------------

Assigned group numbers:

    RANGE           SUBSYSTEM                     OWNER

   0 -   99 - reserved for debugging purposes    C.Olchanski
 100 -  199 - TCYL data structures               C.Olchanski
 200 -  299 - LGD data structures                Rob
 300 -  399 - CSI data structures                ND?
 400 -  499 - MPS modules                        C.Olchanski
 500 -  599 - magnet DVMs                        C.Olchanski

**********************************************************************/

/*
 * Note: This file is used to automatically generate little <--> big
 *       endian data conversion routines. For each data group,
 *       it needs to know the data type and the include file
 *       where the data type is defined.
 *
 * Therefore each data group should be defined using the following
 *       format:
 *
 * #define GROUP_XXX        NUM  ** header.h:type_t - user comments **
 *
 *
 * There exist "special" predefined types and header files that can be used
 *       when there is no "real" data type or header file corresponding
 *       to the data group:
 *
 *   NOCONVERT_t --- a data type that does not need endian conversion,
 *                        for example an array of char.
 *   INTARRAY_t  --- a data type that should be converted as an array
 *                        of 32-bit integers, for example the raw data group.
 *
*/


#define GROUP_DELETED         0  /* groups deleted with data_removeGroup()    disData.h:NOCONVERT_t */
#define GROUP_RAW             1  /* raw data group, see rawData.h, INTARRAY_t */

#define GROUP_TIME0_TDCS     30  /* time zero for this event (time pedestal) (itypes.h:tdc_values_t) */
#define GROUP_TRIGGER_TDCS   40  /* trigger signals timing tdcs (itypes.h:tdc_values_t) */
#define GROUP_TRIGGER_TDCS_COMP 41 /* itypes.h:tdc_values_t - Has only c1,c2,and c3 */
#define GROUP_SCALERS        50  /* scaler_values_t in itypes.h */
#define GROUP_LATCHES        52  /* latch_values_t  in itypes.h */

#define GROUP_MISC_ADCS      54   /* itypes.h:adc_values_t */
#define GROUP_DIBBUK         55   /* dibbuk_data_t, see dibbuk2.h and rawUnpack.c:unpack_dibbuk() */
#define GROUP_EVBV_ENCODERS  56   /* itypes.h:adc_values_t    */
#define GROUP_TEMP_PROBES    57   /* itypes.h:adc_values_t    */
#define GROUP_DIBBUKnfs      58   /* dibbuk2nfs.h:dibbuk2nfs_t dibbuk data read through NFS */

#define GROUP_BOSSUMMARY     61   /* bosSummaryRecord_t, see bosSummary.h */

#define GROUP_VETO_HITS      71   /* adc_hits_t     in itypes.h */
#define GROUP_VETO_BITMAP    72   /* vetoCounters_t in itypes.h */

/*#define GROUP_TCYL          100*/
#define GROUP_TCYL_ADCS     110   /* itypes.h:adc_values_t */
#define GROUP_TCYL_TDCS     111   /* itypes.h:tdc_values_t */
#define GROUP_TCYL_COORS_A  120   /* tcyl.h:tcylHits_t  */
#define GROUP_TCYL_COORS_B  121   /* tcyl.h:tcylHits_t  */
#define GROUP_TCYL_COORS_C  122   /* tcyl.h:tcylHits_t  */
#define GROUP_TCYL_COORS_D  123   /* tcyl.h:tcylHits_t  */
#define GROUP_TCYL_COORS    GROUP_TCYL_COORS_A
#define GROUP_TCYL_PATREC   130   /* tcyl.h:tcylPatrec_hits_t  */
#define GROUP_TCYL_TRACKS   140   /* tcyl.h:tcylTracks_t  */
#define GROUP_TCYL_TRACKCOV 141   /* tcyl.h:tcylTracksCov_t  (covariance matrices for TCYL tracks) */

/*#define GROUP_LGD           200*/
#define GROUP_LGD_MCHITS    209        /* lgdCluster.h:lgd_hits_t  */
#define GROUP_LGD_ADCS      210        /* itypes.h:adc_values_t    */
                                    /* These are defined in lgdCluster.h */
#define GROUP_LGD_HITS      211     /* lgdCluster.h:lgd_hits_t -- this group duplicates GROUP_LGD_ENG below */
#define GROUP_LGD_HITS2     212     /* itypes.h:adc_hits_t        */
#define GROUP_LGD_SHOWER_MC 213     /* LgdShowerLibrary.h:lgd_showerHits_t  */

#define GROUP_LGD_ROCS          240  /* itypes.h:lgd_roc_t */
#define GROUP_LGD_MAM           241  /* itypes.h:lgd_mam_t */

#define GROUP_LGD_CLUSTERS      252  /* lgdCluster.h:lgd_clusters_t    */
#define GROUP_LGD_CLUSTER_HITS  253  /* lgdCluster.h:lgd_hits_t        */
#define GROUP_LGD_ENCODERS      254  /* itypes.h:adc_values_t      */
#define GROUP_LGD_CLUSTER_HITS2 255  /* itypes.h:adc_hits_t        */
#define GROUP_LGD_CLUSTERS_MC   256  /* lgdCluster.h:lgd_clusters_t -- generated by montecarlo */

#define GROUP_NN_CLUSTERS       260  /* lgdCluster.h:lgd_clusters_t  */      
#define GROUP_NN_CLUSTER_HITS   261  /* lgdCluster.h:lgd_hits_t      */


/*#define GROUP_CSI           300*/
#define GROUP_CSI_ADCS      310     /* itypes.h:adc_values_t     */
#define GROUP_CSI_SCALERS   311     /* itypes.h:scaler_values_t  */
#define GROUP_CSI_HITS      312     /* itypes.h:csi_hits_t       */
#define GROUP_CSI_CLUSTERS  313	    /* CsI.h:csi_clusters_t*/
#define GROUP_CSI_PACKED    314     /* CsI.h:csi_packed_adcs_t */
#define GROUP_CSI_TAGS      315     /* tags for the CsI blocks. csi_tags_t from itypes.h */
#define GROUP_CSI_TDCS      316     /* itypes.h:tdc_values_t     */

/* General MPS PWC group */

#define GROUP_PWC_PLANES    401     /* pwc_planes_t, see itypes.h for definitions */
#define GROUP_PWC_HITS      402     /* pwc_hits_t, see itypes.h for definitions */

/* Beam PWCs */

#define GROUP_BEAMPWC_PLANES   411  /* pwc_planes_t, see itypes.h for definitions */
#define GROUP_BEAMPWC_HITS     412  /* pwc_hits_t, see itypes.h for definitions   */

/* Beam Hodoscopes */

#define GROUP_BEAMHODO_TDCS    415  /* tdc_values_t, itypes.h */

#define GROUP_BEAMHODO_PLANES  414  /* pwc_planes_t, see itypes.h for definitions */
#define GROUP_BEAMHODO_HITS    416  /* pwc_hits_t, see itypes.h for definitions */

/* The BEAMSCINT groups are a misnomer- the correct name is 'BEAMHODO' */

#define GROUP_BEAMSCINT_TDCS    GROUP_BEAMHODO_TDCS
#define GROUP_BEAMSCINT_PLANES  GROUP_BEAMHODO_PLANES
#define GROUP_BEAMSCINT_HITS    GROUP_BEAMHODO_HITS

#define GROUP_TPX123_PLANES    421  /* pwc_planes_t, see itypes.h for definitions */
#define GROUP_TPX123_HITS      422  /* pwc_hits_t, see itypes.h for definitions */

#define GROUP_MPSDC_PLANES     431  /* mpsdc_planes_t, see itypes.h for definitions */
#define GROUP_MPSDC_HITS       432  /* mpsdc_hits_t, see itypes.h for definitions */
#define GROUP_MPSDC_COMP       433  /* mpsdc_comp_t see itypes.h */

#define GROUP_TDX4_TDCS        461  /* tdc_values_t see itypes.h */

#define GROUP_COORS            471  /* saved chamber hits, disData.h:INTARRAY_t */

#define GROUP_DVMS             511  /* adc_values_t see itypes.h */

#define GROUP_C9_ADCS          551  /* adc_values_t see itypes.h */
#define GROUP_C9_TDCS          552  /* tdc_values_t see itypes.h */
#define GROUP_C9_HITS          553  /* adc_hits_t see itypes.h   */
#define GROUP_C9_COOKED_HITS   556  /* adc_hits_t, better def'n of hit. 
				       note: fvalue = no. of photons,
				       negative number = don't know. */
#define GROUP_C9_MC_HITS       554  /* adc_hits_t see itypes.h   */
#define GROUP_NOTC9H9          555  /* c9.h: notc9h9_hits_t  */
#define GROUP_H9_ADCS          561  /* adc_values_t see itypes.h */
#define GROUP_H9_TDCS          562  /* tdc_values_t see itypes.h */
#define GROUP_H9_HITS          563  /* adc_hits_t see itypes.h   */
#define GROUP_H9_MC_HITS       564  /* adc_hits_t see itypes.h   */

/* MonteCarlo */

#define GROUP_MC_TRIGGERMASK   601  /* mcTriggerMask.h: mcTriggerMask_t */

/* Analysis */

#define GROUP_MONTECARLO_EVENT 700 /* event_t with SAGE Monte Carlo event, see cstruct.h */

#define GROUP_BEAM_RECORD       710 /* see: mpsX.h BeamRec_t */
/*#define GROUP_BEAM_TRACKCOV     711*/ /* beam track covariance matrix: see: beam.h BeamCov_t */

/*
 * The TRACKS groups below are the reconstructed tracks
 * before the vertex is located
*/

#define GROUP_MONTECARLO_TRACKS 721 /* tracks_t (see tracking.h) */
#define GROUP_TRACKS            722 /* tracks_t (see tracking.h) */
#define GROUP_HTR_TRACKS        723 /* tracks_t (see tracking.h), see the HTR package */

#define GROUP_TrackPositions    724 /* trackPositions_t (see trackPositions.h) */

/* Pattern recognition groups:
     GROUP_PATREC     - output of the pattern recognition program
     GROUP_PATREC_TRK - output of the geometry reconstruction program
     GROUP_PATREC_HTR - output of the HTR geometry package
*/

#define GROUP_MONTECARLO_PATREC 731 /* trackHits_t (see tracking.h) */
#define GROUP_PATREC            732 /* regular pattern recognition data, trackHits_t (see tracking.h) */
#define GROUP_PATREC_TRK        733 /* patrec after fitting,             trackHits_t (see tracking.h) */
#define GROUP_PATREC_HTR        734 /* patrec from HTR,                  trackHits_t (see tracking.h) */
#define GROUP_PATREC_PROJ       735 /* patrec projected tracks,   patrecProjTracks_t (see tracking.h) */
#define GROUP_PATREC_BEAM       736 /* patrec from beam chambers,        trackHits_t (see tracking.h) */

/*
 * The GEO groups are the output of the geometry reconstruction program
 *
 * Note, that the GEO_TRACKS use the vertex as one more measurement
*/

#define GROUP_GEO_MC_TRACKS     741 /* geo_tracks_rec_t (see tracking.h) */
#define GROUP_GEO_TRACKS        742 /* geo_tracks_rec_t (see tracking.h) */

#define GROUP_GEO_MC_VERTICES   751 /* geo_vertices_rec_t (see tracking.h) */
#define GROUP_GEO_VERTICES      752 /* geo_vertices_rec_t (see tracking.h) */

/* Event Summary Records */

#define GROUP_ESR_NPRONG        801 /* esr_nprong_t    (see esr.h) */
#define GROUP_ESR_NPARTICLE	802 /* esr_nparticle_t (see esr.h) */
#define GROUP_ESR_VETO          803 /* esrVeto_t       (see veto.h) */
#define GROUP_NDESR_VETO        GROUP_ESR_VETO


#define GROUP_EVENT_WEIGHTS     804 /* EventWeights_t  (esr.h) used by PWA */
#define GROUP_PWA_AMPLITUDES    805 /* amplitudes_t    (esr.h) used by PWA */
#define GROUP_ESR_NPARTICLE_MC  806 /* esr_nparticle_t (see esr.h) */ 
#define GROUP_ESR_VERTICES      807 /* esr_vertices_t  (see esr.h) */ 

#define GROUP_ESR_COMPRESSED    808 /* esr_compressed_t (see esr.h) */
#define GROUP_ESR_GENERIC       809 /* generic_t (see esr.h) */ 

#define GROUP_ESR_VERTICES_MC   810 /* esr_vertices_t  (see esr.h) */ 

#define GROUP_ESR_BEFORE_SQUAW	811 /* esr_nparticle_t (see esr.h) */
#define GROUP_ESR_NPARTICLE_ALL	812 /* esr_nparticle_t (see esr.h) */

#define GROUP_SQUAW_SVFITS      840 /* svfits_t   (see squaw.h)   */
#define GROUP_SQUAW_SVPULLS     841 /* svpulls_t  (see squaw.h)   */
#define GROUP_SQUAW_SV          842 /* sqfits_t   (see squaw97.h) */
#define GROUP_SQUAW_SVP         843 /* sqpulls_t  (see squaw97.h) */

#define GROUP_SQUAW_MVFITS      850 /* mvfits_t   (see squaw.h)   */
#define GROUP_SQUAW_MVPULLS     851 /* mvpulls_t  (see squaw.h)   */
#define GROUP_SQUAW_MV          852 /* sqfits_t   (see squaw97.h) */
#define GROUP_SQUAW_MVP         853 /* sqpulls_t  (see squaw97.h) */


#ifdef mips
#define TRANSCOMPUTERCODE 1
#endif

#ifndef TRANSCOMPUTERCODE
#define TRANSCOMPUTERCODE 123
#endif

#include <ntypes.h>

typedef struct
{
  uint32   length;       /* record length in bytes      */
  uint32   type;         /* should always be TYPE_ITAPE */
  uint32   transComputerCode;
  uint32   ngroups;
   int32   runNo;
  uint32   spillNo;
   int32   eventNo;      /* 0-BOS, (-1)-EOS, others- events number */
  uint16   eventType;    /* event type (decoded latch word)        */
   int16   trigger;      /* trigger: 0: BOS, 1-6: events, -1: EOS  */
   int32   time;         /* time as returned by time()             */
  uint32   latch;        /* virtual latch word                     */
} itape_header_t;

typedef struct
{
  uint32   length;
  uint32   type;
} group_header_t;


/*
 * 'data_newItape' will fill initialize an empty itape
*/

int data_newItape(void* itape);

/*
 * 'data_getGroup' returns the pointer to the group data or NULL
*/

void* data_getGroup(const void*buf,uint32 group);

unsigned long data_getGroupSize(const void*buf,uint32 group);

/*
 * 'data_getGroupHeader' returns thr pointer to the group header or NULL
*/

group_header_t* data_getGroupHeader(const void*buf,uint32 type);

void* data_addGroup(void* buffer,int bufsize,uint32 group,uint32 length);

group_header_t* data_addGroup_header(void* buffer,int bufsize,
				     uint32 group,uint32 length);

/*
 * 'data_renameGroup' will rename a given group. The return value is:
 * 0: success, non 0: error.
*/

int data_renameGroup(void* event,int oldName,int newName);

/*
 * 'data_removeGroup' returns 0 if group found and removed, 1 otherwise
 *    the space occupied by the removed group is not deallocated and is lost
 *    until 'data_clean' is called to compact the itape
*/

int data_removeGroup(void* buffer,uint32 groupName);


int data_removeAllgroups(void* buffer);  /* returns 1 */
int data_saveGroups(void* buffer,int nsave,int *isave);  /* returns 1 */

/*
 * returns the address of the end of used space inside the buffer
*/

void* data_getAdr(const void*buffer);

/*
 * 'data_clean' removes the unused parts of the itape created by the 'removeGroup' functions
 *
 *    NOTICE: this call moves data inside the itape and therefore all pointers
 *            into the itape (for example, pointers returned by 'getGroup')
 *            become invalid.
*/

int data_clean(void* buffer);

/*
 * Calculates a new CRC for the itape
*/

int data_addCRC(void*buffer);

/*
 * Check the CRC. If it is broken, return 1, otherwise, return 0
*/

int data_checkCRC(const void*buffer);
int data_checkCRC1(const void*buffer,int bufferLength);

/*
 * perform big <--> little endian conversion
*/

int endian_convertItape(void *itape,int itapeLength);

/* tmask() is used to check the trigger word in the itape against a trigger mask.
   The trigger word is an integer from -1 to 6 set by the DAQ code in the master SSP
   according to the trigger latch. It corresponds to the 8 possible triggers
   implemented in the E852 DAQ system.

   To find out which trigger latch values correspond to which trigger, see the DAQinit
   database (default.db, current.db, the currently selected db, etc...).

   Each bit in the trigger mask selects one of the triggers:

       bit 0 (0x01) - trigger -1 - EOS
       bit 1 (0x02) - trigger  0 - BOS
       bit 2 (0x04) - trigger 1
       bit 3 (0x08) - trigger 2
       bit 4 (0x10) - trigger 3
       bit 5 (0x20) - trigger 4
       bit 6 (0x40) - trigger 5
       bit 7 (0x80) - trigger 6
*/

int tmask(unsigned long mask,unsigned long trigger);

/*
 * To extract an event from a spill in shared memory
*/

void* data_getItape(const void*buffer,const void**nextItape);

/*
 * The 'data_listGroups' function will return-
 *    in 'ngroups' - number of groups in the event
 *    in 'groupsList' - a pointer to the array of groups
 *    in 'groupsSize' - a pointer to the array of their lengths
 *
 * All the three arguments can be NULL if the returned information is not needed
 * The arrays should be free'ed after use
*/

int data_listGroups(const void*event,int*ngroups,int **groupsList,int **groupsSize);

#endif
/* end file */
