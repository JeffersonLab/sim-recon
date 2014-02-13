/*
 * itypes.h - unpacked data - types (after the unpacker)
*/

#ifndef ITYPES_INCLUDED
#define ITYPES_INCLUDED

static const char sccsid_itypesH[] = "@(#)itypes.h\t5.2\tCreated 7/27/97 18:59:14, \tcompiled "__DATE__;

#include <ntypes.h>      /* for int32, etc... definitions */


typedef struct
{
  uint32 channel;
  uint32 value;
} adc_value_t;

typedef struct
{
  uint32 channel;
  int32  le;         /* leading edge  */
  int32  te;         /* trailing edge */
} tdc_value_t;

typedef struct
{
  int32 nadc;
  adc_value_t adc[1];
} adc_values_t;

typedef struct
{
  int32 ntdc;
  tdc_value_t tdc[1];
} tdc_values_t;

typedef struct
{
  int32  nscalers;
  uint32 scaler[1];
} scaler_values_t;

typedef struct
{
  int32  nlatches;
  uint32 latch[1];
} latch_values_t;

typedef struct
{
  int32   channel;
  float32 fvalue;
} oldAdc_hit_t;

typedef struct
{
  int32   channel;
  int32   adc;
  float32 fvalue;
} adc_hit_t;

typedef struct
{
  int32 nhits;
  adc_hit_t hit[1];
} adc_hits_t;


/*
  These records go into the groups created by PWC, and drift chamber unpacking code
*/

typedef struct
{
  int32 plane;  /* plane id, same as in the 'space' printout (i.e. counting from one) */
  int32 first;  /* index of the hit in the hits group */
  int32 nhits;  /* number of hits on this plane */
} mps_plane_t;

typedef struct
{
  int32 nplanes;
  mps_plane_t plane[1];
} mps_planes_t;

typedef struct
{
  int32 first;   /* first wire in cluster */
  int32 last;    /* last wire in the cluster -
		  for one-cluster hits first and last will be equal */
} pwc_cluster_t;

typedef struct
{
  int32 nclusters;
  pwc_cluster_t cluster[1];
} pwc_hits_t;

typedef mps_plane_t pwc_plane_t;
typedef mps_planes_t pwc_planes_t;

typedef struct
{
  int32 nwire;      /* struck wire in scintillator    */
  int32 tdc_count;  /* count on the corresponding tdc */
} scint_hit_t;

typedef struct
{
  int32 nhits;
  scint_hit_t hit[1];
} scint_hits_t;

typedef mps_plane_t scint_plane_t;
typedef mps_planes_t scint_planes_t;

/*
  These records constitute groups created by MPS Drift Chambers
  unpacking code

  GROUP_MPSDC_HITS is made of one mpsdc_hits_t structure and 
                 many mpsdc_hit_t structures - one for every hit on every plane.

  GROUP_MPSDC_PLANES is made of one mpsdc_planes_t record and many
                 mpsdc_plane_t record - one for every MPS DC plane

*/

typedef struct
{
  int32 wire;     /* wire number (as given by the MPS DC readout system */
  int32 drift;    /* drift time (in MPS DC ticks) */
} mpsdc_hit_t;


typedef struct
{
  int32 nhits;
  mpsdc_hit_t hit[1];
} mpsdc_hits_t;


typedef mps_plane_t mpsdc_plane_t;
typedef mps_planes_t mpsdc_planes_t;


typedef struct{
  uint8  plane;
  uint8  time;
  uint16 wire;
} mpsdc_hit_comp_t;

typedef struct{
  int32 nhits;
  mpsdc_hit_comp_t hit[1];
} mpsdc_comp_t;


/* CsI data structures */

typedef struct {
  int32   block;
  float32 energy;
} csi_hit_t;

typedef struct {
  int32  nhits;
  int32  HaloLatch;
  uint32 DiscLatch[7];
  csi_hit_t hits[1];
} csi_hits_t;


typedef struct
{
  int32  nhits;
  uint32 tags[1];

} csi_tags_t;

/*
 * Veto counters bit maps
*/

typedef struct
{
  int32 length;         /* size of the structure, including the length word */
  int32 flags;
  uint32 cpva;
  uint32 cpvb;
  uint32 cpvc;
  uint32 dea;
  uint32 ev;
  uint32 bv;

} vetoCounters_t;

/*
 * For LGD trigger processor readout
*/

typedef struct{
  int32 nwords;
  int32 macE;
  int32 macPhi;
  int32 macPy;
  int32 macPx;
  int32 cnt;
  int32 max;
  int32 mult;
  int32 etot;
  int32 massSqrd;
  int32 t;
  int32 input[1];
} mam_indiv_t;

typedef struct{
  int32 nwords;
  int32 data[1];
} lgd_tp_list_t;

typedef lgd_tp_list_t lgd_roc_t;
typedef lgd_tp_list_t mam_list_t;

typedef union{
  mam_indiv_t indiv;
  mam_list_t list;
}lgd_mam_t;

#endif
/* end file */
