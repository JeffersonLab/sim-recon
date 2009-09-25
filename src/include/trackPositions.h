/*
 * trackPositions.h
*/


typedef struct
{
   char name[8]; /* name of the plane - for example "LGD"
		  * the name has to be a miltiple of 4 bytes
		  * to avoid alignment-to-word-boundary problems */
   int itrack;   /* index into the GROUP_TRACKS           */
   vector3_t v;  /* particle position */
   vector3_t p;  /* particle momentum */
   int charge;   /* particle charge   */
   int flags;    /* flags, currently set to 0 */
} trackPosition_t;


typedef struct
{
  int nentries;
  trackPosition_t entry[1];
} trackPositions_t;


void MakeTrackPositions(itape_header_t*event,int eventSize);

/* end file */
