/* mc2coda Library
 Includes module data generating functions */

 #include <DAQ/daq_param_type.h>
 
void GetPedestals(uint32_t *peds, uint32_t Npeds);

#define SUPPRESS_DRIFT_CHAMBER_HITS_OVERFLOW_WARNINGS 1


/* FADC 250 Paramters */
#define FADC250_MAX_CHAN      16
#define FADC250_MAX_HITS       4
#define FADC250_MAX_WINDOW    512
#define FADC250_MAX_LATENCY   2048
#define FADC250_MAX_NSB       512
#define FADC250_MAX_NSA       512
#define FADC250_WINDOW_WIDTH  50

/* FADC 250 Macros */
#define FADC250_BL_HEADER(slot,blnum,cnt) {*dabufp++ =  0x80000000 | (slot << 22) | (FADC250 << 18) | ((blnum&0x3ff) << 8) | cnt; }
#define FADC250_BL_TRAILER(slot,nwords)   {*dabufp++ =  0x88000000 | (slot << 22) | nwords; }

#define FADC250_EV_HEADER(slot,trig)  {*dabufp++ = 0x90000000 | (slot << 22) | (trig&0x3fffff); }
#define FADC250_EV_TS_LOW(timestamp)  {*dabufp++ = 0x98000000 | (timestamp&(0x0000000000ffffffLLU)); }
#define FADC250_EV_TS_HIGH(timestamp) {*dabufp++ = (timestamp&(0x0000ffffff000000LLU)) >> 24; }
#define FADC250_EV_PARAM1(nsb,nsa,pl) {*dabufp++ = (nsb) << 20 | (nsa) << 11 | pl; }

#define FADC250_RW_HEADER(chan,ww) {*dabufp++ = 0xA0000000 | (chan << 23) | ww ; }
#define FADC250_RW_DATA(s1,s2) {*dabufp++ = (s1 << 16) | s2 ; }

#define FADC250_WI_SUM(chan,sum) {*dabufp++ = 0xA8000000 | (chan << 23) | (sum&0x3fffff) ; }

#define FADC250_RP_HEADER(chan,pn,fs) {*dabufp++ = 0xB0000000 | (chan << 23) | (pn << 21) | (fs) ; }
#define FADC250_RP_DATA(s1,s2) {*dabufp++ = (s1 << 16) | s2 ; }

#define FADC250_PI_SUM(chan,pn,sum)   {*dabufp++ = 0xB8000000 | (chan << 23) | (pn << 21) | (sum&0x7ffff) ; }
#define FADC250_PI_TIME(chan,pn,time) {*dabufp++ = 0xC0000000 | (chan << 23) | (pn << 21) | (time&0x7ffff); }
#define FADC250_PI_PED(chan,pn,ped,peak) {*dabufp++ = 0xD0000000 | (chan << 23) | (pn << 21) | ((ped&0x1ff)<<12) | (peak&0xfff); }

#define FADC250_FILLER {*dabufp++ = 0xF8000000; }


int
fadc250_write_data (CODA_EVENT_INFO *event, int roc, int slot, int mode)
{
   
   int ii, jj, chan, hcnt, nwords;
   uint32_t  eventNum;
   uint64_t  timestamp;
   unsigned int *start = dabufp;
   CODA_HIT_INFO *hit;
   CODA_HIT_INFO *chit[MAX_HITS_PER_CHAN];
   
   eventNum  = (event->eventid)&0xffffffff;
   timestamp = (event->trigtime);
   hcnt      =  event->hcount[(roc-1)][(slot-1)];
   
   FADC250_BL_HEADER(slot,eventNum,1);
   FADC250_EV_HEADER(slot,eventNum);
   FADC250_EV_TS_LOW(timestamp);
   
   /*printf("fadc250:_write_data: DEBUG: %d hits available"
    *       " for  slot %d roc %d.\n",hcnt, slot, roc);
    */
	
	// Get pedestal values (possibly randmized)
	uint32_t peds[FADC250_MAX_CHAN+1];
	GetPedestals(peds, FADC250_MAX_CHAN+1);
   
   /*Loop over all channels */
   for (chan=0; chan<FADC250_MAX_CHAN; chan++) {
      
      /* check for all hits for this channel */
      jj=0;
      for (ii=0; ii<hcnt; ii++) {
         hit = (CODA_HIT_INFO *) &event->hits[(roc-1)][(slot-1)][ii];
         /* printf("%d %d %d %d\n",hit->crate_id,hit->slot_id,
          *                        hit->chan_id,hit->module_id);
          */
         if ( (roc == hit->crate_id) && (slot == hit->slot_id) && 
              (chan == hit->chan_id) && (hit->module_id == FADC250) &&
              (hit->module_mode == mode) )
         {
            /* printf("fadc250:_write_data: DEBUG: found hit %d for roc"
             *        " %d slot %d chan %d.\n",jj, roc,slot,chan);
             */
            
            chit[jj] = hit;
            jj++;
            if (jj >= MAX_HITS_PER_CHAN) {
               printf("fadc250_write_data: ERROR: HIT overflow (%d) for "
                      "crate,slot,chan = %d,%d,%d (Bank will be corrupt)\n",
                      jj,roc,slot,chan);
               return 0;
            }
            else if (jj > FADC250_MAX_HITS) {
               printf("fadc250_write_data: WARN: Too many hits (%d) for"
                      " (crate, slot, chan) = %d, %d, %d (truncating)\n",
                      jj,roc,slot,chan);
               jj = FADC250_MAX_HITS;
            }
         }
      }
      /* printf("write hit data %d\n",jj); */
      uint32_t ped = peds[FADC250_MAX_CHAN]; // common pedestal for all channels
      for (ii=0; ii < jj; ii++) {
		 uint32_t peak = ped + chit[ii]->hdata[0]/10; // Assume pulse is 20samples(=80ns) so area is ~0.5*b*h=pi -> h=pi/(b/2) where b=20samples
		 uint32_t sum = chit[ii]->hdata[0] + ped*FADC250_WINDOW_WIDTH;
         FADC250_PI_SUM(chan,ii,sum);
         FADC250_PI_TIME(chan,ii,chit[ii]->hdata[1]);
	 FADC250_PI_PED(chan,ii,(ped + peds[chan]),peak); // measured pedestal has additional stochastic element
      }
   }
   
   nwords = dabufp - start;
   if ((nwords%2) == 0) {
      FADC250_FILLER;
      nwords += 2;
      FADC250_BL_TRAILER(slot,nwords);
   } else {
      nwords += 1;
      FADC250_BL_TRAILER(slot,nwords);
   }
   
   if (nwords%4 != 0) {
      FADC250_FILLER;
      FADC250_FILLER;
      nwords += 2;
   }
   
   return nwords;
}


/* FADC 125 Paramters */
#define FADC125_MAX_CHAN      72
#define FADC125_MAX_HITS       4
#define FADC125_MAX_WINDOW    512
#define FADC125_MAX_LATENCY   2048
#define FADC125_MAX_NSB       512
#define FADC125_MAX_NSA       512
#define FADC125_WINDOW_WIDTH  50

/* FADC 125 Macros */
#define FADC125_BL_HEADER(slot,blnum,cnt) {*dabufp++ =  0x80000000 | (slot << 22) | (FADC125 << 18) | ((blnum&0x3ff) << 8) | cnt; }
#define FADC125_BL_TRAILER(slot,nwords)   {*dabufp++ =  0x88000000 | (slot << 22) | nwords; }

#define FADC125_EV_HEADER(slot,trig)  {*dabufp++ = 0x90000000 | (slot << 22) | (trig&0x3fffff); }
#define FADC125_EV_TS_LOW(timestamp)  {*dabufp++ = 0x98000000 | (timestamp&(0x0000000000ffffffLLU)); }
#define FADC125_EV_TS_HIGH(timestamp) {*dabufp++ = (timestamp&(0x0000ffffff000000LLU)) >> 24; }

#define FADC125_RW_HEADER(chan,ww) {*dabufp++ = 0xA0000000 | (chan << 20) | (ww&0x1fffff) ; }
#define FADC125_RW_DATA(s1,s2) {*dabufp++ = (s1 << 16) | s2 ; }

#define FADC125_WI_SUM(chan,sum) {*dabufp++ = 0xA8000000 | (chan << 20) | (sum&0x1fffff) ; }

#define FADC125_RP_HEADER(chan,pn,fs) {*dabufp++ = 0xB0000000 | (chan << 20) | (pn << 18) | (fs&0x3ffff) ; }
#define FADC125_RP_DATA(s1,s2) {*dabufp++ = (s1 << 16) | s2 ; }

#define FADC125_PI_SUM(chan,sum)   {*dabufp++ = 0xB8000000 | (chan << 20) | (sum&0xfffff) ; }
#define FADC125_PI_TIME(chan,qf,time) {*dabufp++ = 0xC0000000 | (chan << 20) | (qf << 18) | (time&0xffff); }
#define FADC125_PI_PED(chan,pn,ped,peak) {*dabufp++ = 0xD0000000 | ((chan << 23)&0x0f) | (pn << 21) | ((ped&0x1ff)<<12) | (peak&0xfff) ; }

#define FADC125_FILLER {*dabufp++ = 0xF8000000; }


int
fadc125_write_data (CODA_EVENT_INFO *event, int roc, int slot, int mode)
{
   
   int ii, jj, chan, hcnt, nwords;
   uint32_t  eventNum;
   uint64_t  timestamp;
   unsigned int *start = dabufp;
   CODA_HIT_INFO *hit;
   CODA_HIT_INFO *chit[MAX_HITS_PER_CHAN];
   
   eventNum  = (event->eventid)&0xffffffff;
   timestamp = (event->trigtime);
   hcnt      =  event->hcount[(roc-1)][(slot-1)];
   
   /* Global timestamp is in 4ns ticks. The local clock on the FADC 125
    * is half that so change timestamp to 8ns ticks (divide by two).
    */
   timestamp = timestamp >> 1;

	// Get pedestal values (possibly randmized)
	uint32_t peds[FADC125_MAX_CHAN+1];
	GetPedestals(peds, FADC125_MAX_CHAN+1);

   FADC125_BL_HEADER(slot,eventNum,1);
   FADC125_EV_HEADER(slot,eventNum);
   FADC125_EV_TS_LOW(timestamp);
   
   /*  printf("fadc125_write_data: %d hits checking slot %d roc %d.\n",
    *         hcnt, slot, roc);
    */
   
   /*Loop over all channels */
   for (chan=0; chan < FADC125_MAX_CHAN; chan++) {
      
      /* check for all hits for this channel */
      jj=0;
      for (ii=0; ii < hcnt; ii++) {
         hit = (CODA_HIT_INFO *) &event->hits[(roc-1)][(slot-1)][ii];
         /* printf("%d %d %d %d\n",hit->crate_id,hit->slot_id,
          *        hit->chan_id,hit->module_id);
          */
         if ( (roc == hit->crate_id) && (slot == hit->slot_id) &&
              (chan == hit->chan_id) && (hit->module_id == FADC125) &&
              (hit->module_mode == mode) )
         {
            /* printf("fadc125:_write_data: DEBUG: found hit for crate"
             *        " %d slot %d chan %d.\n",roc,slot,chan);
             */
            chit[jj] = hit;
            jj++;
            if (jj >= MAX_HITS_PER_CHAN) {
               printf("fadc125_write_data: ERROR: HIT overflow (%d)"
                      " for (crate, slot, chan) = %d, %d, %d\n",
                      jj,roc,slot,chan);
               printf("fadc125_write_data: ERROR: "
                      "ROC Bank will be corrupted\n");
               return 0;
            }
            else if (jj > FADC125_MAX_HITS) {
#ifndef SUPPRESS_DRIFT_CHAMBER_HITS_OVERFLOW_WARNINGS
               printf("fadc125_write_data: WARN: Too many hits (%d)"
                      " for (crate, slot, chan) = %d, %d, %d\n",
                      jj,roc,slot,chan);
#endif
               jj = FADC125_MAX_HITS;
            }
         }
      }
      /* printf("write hit data %d\n",jj); */
      uint32_t ped = peds[FADC125_MAX_CHAN]; // common pedestal for all channels
      for (ii=0; ii < jj; ii++) {
	      uint32_t peak = ped + chit[ii]->hdata[0]/10; // Assume pulse is 20samples(=160ns) so area is ~0.5*b*h=pi -> h=pi/(b/2) where b=20samples
	      uint32_t sum = chit[ii]->hdata[0] + ped*FADC125_WINDOW_WIDTH;
	      FADC125_PI_SUM(chan, sum);
	      FADC125_PI_TIME(chan,0,chit[ii]->hdata[1]); // always make quality factor 0
	      FADC125_PI_PED(chan,ii,(ped + peds[chan]),peak); // measured pedestal has additional stochastic element
	      // f125 generally doesn't support multiple pulses so only record first one
	      break;
      }
   }
   
   nwords = dabufp - start;
   if ((nwords%2) == 0) {
      FADC125_FILLER;
      nwords += 2;
      FADC125_BL_TRAILER(slot,nwords);
   } else {
      nwords += 1;
      FADC125_BL_TRAILER(slot,nwords);
   }
   
   if (nwords%4 != 0) {
      FADC125_FILLER;
      FADC125_FILLER;
      nwords += 2;
   }
   
   return nwords;
}


/* F1TDC 32 Channel Parameters */
#define F1TDC32_MAX_CHAN      32
#define F1TDC32_MAX_HITS       8
#define F1TDC32_MAX_CHIPS      8


/* F1TDC 32 channel Macros */
#define F1TDC32_BL_HEADER(slot,blnum,cnt) {*dabufp++ =  0x80000000 | (slot << 22) | (F1TDC32 << 18) | ((blnum&0x3ff) << 8) | cnt; }
#define F1TDC32_BL_TRAILER(slot,nwords)   {*dabufp++ =  0x88000000 | (slot << 22) | nwords; }

#define F1TDC32_EV_HEADER(slot,trig)  {*dabufp++ = 0x90000000 | (slot << 22)| (trig&0x3fffff); }
#define F1TDC32_EV_TS_LOW(timestamp)  {*dabufp++ = 0x98000000 | (timestamp&(0x0000000000ffffffLLU)); }
#define F1TDC32_EV_TS_HIGH(timestamp) {*dabufp++ = (timestamp&(0x000000ffff000000LLU)) >> 24; }


// These changed based on the F1TDC_V2_V3_4_29_14.pdf document
//#define F1TDC32_F1_HEADER(cdata,evt,trig,chan)  {*dabufp++ = 0xC0000000 | (cdata << 24) | ((evt&0x1ff) << 16) | (chan)       ; }
//#define F1TDC32_F1_DATA(cdata,chan,time)        {*dabufp++ = 0xB8000000 | (cdata << 24) | (chan << 16)        | (time&0xffff); }
//#define F1TDC32_F1_TRAILER(cdata,evt,trig,chan) {*dabufp++ = 0xA8000000 | (cdata << 24) | ((evt&0x1ff) << 16) | (chan)       ; }
#define F1TDC32_F1_HEADER(cdata,chip,chan_on_chip,trig,trig_time)  {*dabufp++ = 0xC0000000 | ((cdata&0x1F) << 22) | (0 << 23) | (chip << 3) | (chan_on_chip << 0) | ((trig&0x3f) << 16) | ((trig_time&0x1ff) << 7); }
#define F1TDC32_F1_DATA(cdata,chip,chan_on_chip,time)              {*dabufp++ = 0xB8000000 | ((cdata&0x1F) << 22) | (1 << 23) | (chip << 19) | (chan_on_chip << 16) | (time&0xffff); }

#define F1TDC32_FILLER(slot) {*dabufp++ = 0xF8000000 | (slot << 22); }

#define F1TDC32_CHIP_NUM(chan) (chan >> 2)
#define F1TDC32_CHAN_ON_CHIP(chan) ((chan & 0x03) << 1)

int
f1tdc32_write_data (CODA_EVENT_INFO *event, int roc, int slot, int mode)
{
   
   int ii, jj, chan, hcnt, nwords;
   int chip, chan_on_chip;
   uint64_t tsdiv;
   uint32_t ts, cdata;
   uint32_t  eventNum;
   uint64_t  timestamp;
   unsigned int *start = dabufp;
   CODA_HIT_INFO *hit;
   CODA_HIT_INFO *chit[MAX_HITS_PER_CHAN];
   
   eventNum  = (event->eventid)&0xffffffff;
   timestamp = (event->trigtime);
   hcnt      =  event->hcount[(roc-1)][(slot-1)];
   
   /* Set default value for cdata bits - 3 bits - 100b = 0x4
    *  res locked, ouput fifo ok, hit fifo ok 
    */
   cdata = 0x4;
   
   /* Timestamp is in 4 ns ticks. We need to convert to 
    * F1 clocks = 250/8 = 31.25MHz = 32ns/tick 
    */
   tsdiv = (timestamp >> 3);
   ts = tsdiv&0xffffffff;
   
   
   F1TDC32_BL_HEADER(slot,eventNum,1);
   F1TDC32_EV_HEADER(slot,eventNum);
   F1TDC32_EV_TS_LOW(tsdiv);
   
   /* printf("f1tdc32:_write_data: %d hits checking slot"
    *        " %d roc %d.\n",hcnt, slot, roc);
    */
   
   /*Loop over all channels */
   for (chan=0; chan<F1TDC32_MAX_CHAN; chan++) {

      chip = F1TDC32_CHIP_NUM(chan);
      chan_on_chip = F1TDC32_CHAN_ON_CHIP(chan);
      
      /* Check for outputing Chip Header */
      if (chan_on_chip == 0) {
         //F1TDC32_F1_HEADER(cdata, (eventNum&0x3f),(ts&0x1ff),chan);
         F1TDC32_F1_HEADER(cdata,chip, 7,(eventNum&0x3f), (ts&0x1ff));
      }
            
      /* check for all hits for this channel */
      jj=0;
      for (ii=0; ii<hcnt; ii++) {
         hit = (CODA_HIT_INFO *) &event->hits[(roc-1)][(slot-1)][ii];
         /* printf("%d %d %d %d\n",hit->crate_id,hit->slot_id,
          *        hit->chan_id,hit->module_id);
          */
         if ( (roc == hit->crate_id) && (slot == hit->slot_id) &&
              (chan == hit->chan_id)&& (hit->module_id == F1TDC32) &&
              (hit->module_mode == mode) )
         {
            /* printf("f1tdc32:_write_data: DEBUG: found hit for crate %d"
             *       " slot %d chan %d.\n",roc,slot,chan); 
             */
            
            chit[jj] = hit;
            jj++;
            if (jj >= F1TDC32_MAX_HITS) {
               printf("f1tdc32_write_data: ERROR: Too many hits for channel\n");
               jj--;
					break;
            }
         }
      }
      /* printf("write hit data %d\n",jj); */
      for (ii=0; ii < jj; ii++) {
         /* printf("f1tdc32:_write_data: DEBUG: found hit for crate %d"
          *        " slot %d chan %d  (chip %d chan on chip %d).\n",
          *        roc,slot,chan,chip,chan_on_chip); 
          */
          F1TDC32_F1_DATA(cdata,chip,chan_on_chip,chit[ii]->hdata[0]);
      }
      
      /* Check for outputing chip trailer */
#if 0
      if (chan == 31) {
         F1TDC32_F1_TRAILER(cdata,(eventNum&0x3f),(ts&0x1ff),chan);
      }
#endif   
      
   }
   
   nwords = dabufp - start;
   if ((nwords%2) == 0) {
      F1TDC32_FILLER(slot);
      nwords += 2;
      F1TDC32_BL_TRAILER(slot,nwords);
   } else {
      nwords += 1;
      F1TDC32_BL_TRAILER(slot,nwords);
   }
   
   if (nwords%4) {
      F1TDC32_FILLER(slot);
      F1TDC32_FILLER(slot);
      nwords += 2;
   }
   
   return nwords;
}


/* F1TDC 48 Channel Parameters */
#define F1TDC48_MAX_CHAN      48
#define F1TDC48_MAX_HITS       8
#define F1TDC48_MAX_CHIPS      6


/* F1TDC 48 channel Macros */
#define F1TDC48_BL_HEADER(slot,blnum,cnt) {*dabufp++ =  0x80000000 | (slot << 22) | (F1TDC48 << 18) | ((blnum&0x3ff) << 8) | cnt; }
#define F1TDC48_BL_TRAILER(slot,nwords)   {*dabufp++ =  0x88000000 | (slot << 22) | nwords; }

#define F1TDC48_EV_HEADER(slot,trig)  {*dabufp++ = 0x90000000 | (slot << 22) | (trig&0x3fffff); }
#define F1TDC48_EV_TS_LOW(timestamp)  {*dabufp++ = 0x98000000 | (timestamp&(0x0000000000ffffffLLU)); }
#define F1TDC48_EV_TS_HIGH(timestamp) {*dabufp++ = (timestamp&(0x0000ffffff000000LLU)) >> 24; }


// These changed based on the F1TDC_V2_V3_4_29_14.pdf document
//#define F1TDC48_F1_HEADER(cdata,evt,trig,chan)  {*dabufp++ = 0xA0000000 | (cdata << 24) | ((evt&0x1ff) << 16) | (chan); }
//#define F1TDC48_F1_DATA(cdata,chan,time)        {*dabufp++ = 0xB0800000 | (cdata << 24) | (chan << 16)        | (time&0xffff); }
//#define F1TDC48_F1_TRAILER(cdata,evt,trig,chan) {*dabufp++ = 0xA8000000 | (cdata << 24) | ((evt&0x1ff) << 16) | (chan); }
#define F1TDC48_F1_HEADER(cdata,chip,chan_on_chip,trig,trig_time)  {*dabufp++ = 0xC0000000 | ((cdata&0x1F) << 22) | (0 << 23) | (chip << 3) | (chan_on_chip << 0) | ((trig&0x3f) << 16) | ((trig_time&0x1ff) << 7); }
#define F1TDC48_F1_DATA(cdata,chip,chan_on_chip,time)              {*dabufp++ = 0xB8000000 | ((cdata&0x1F) << 22) | (1 << 23) | (chip << 19) | (chan_on_chip << 16) | (time&0xffff); }

#define F1TDC48_FILLER(slot) {*dabufp++ = 0xF8000000 | (slot << 22); }

#define F1TDC48_CHIP_NUM(chan) (chan >> 3)
#define F1TDC48_CHAN_ON_CHIP(chan) (chan & 0x07)


int
f1tdc48_write_data (CODA_EVENT_INFO *event, int roc, int slot, int mode)
{
   
   int ii, jj, chan, hcnt, nwords;
   int chip, chan_on_chip;
   uint64_t tsdiv;
   uint32_t ts, cdata;
   uint32_t  eventNum;
   uint64_t  timestamp;
   unsigned int *start = dabufp;
   CODA_HIT_INFO *hit;
   CODA_HIT_INFO *chit[MAX_HITS_PER_CHAN];
   
   eventNum  = (event->eventid)&0xffffffff;
   timestamp = (event->trigtime);
   hcnt      =  event->hcount[(roc-1)][(slot-1)];

   /* Set default value for cdata bits - 3 bits - 100b = 0x4
    * res locked, ouput fifo ok, hit fifo ok 
    */
   cdata = 0x4;
   
   /* Timestamp is in 4 ns ticks. We need to convert to 
    * F1 clocks = 250/8 = 31.25MHz = 32ns/tick 
    */
   tsdiv = (timestamp >> 3);
   ts = tsdiv&0xffffffff;
   
   
   F1TDC48_BL_HEADER(slot,eventNum,1);
   F1TDC48_EV_HEADER(slot,eventNum);
   F1TDC48_EV_TS_LOW(tsdiv);
   
   /*  printf("f1tdc48_write_data: %d hits checking slot %d roc %d.\n",
    *         hcnt, slot, roc); 
    */
   
   /* Loop over all channels */
   for (chan=0; chan < F1TDC48_MAX_CHAN; chan++) {
      
      chip = F1TDC48_CHIP_NUM(chan);
      chan_on_chip = F1TDC48_CHAN_ON_CHIP(chan);

      /* Check for outputing Chip Header */
      if (chan_on_chip == 0) {
         //F1TDC48_F1_HEADER(cdata, (eventNum&0x3f),(ts&0x1ff),0);
         F1TDC48_F1_HEADER(cdata,chip,7,(eventNum&0x3f), (ts&0x1ff));
      }
      
      /* check for all hits for this channel */
      jj=0;
      for (ii=0; ii < hcnt; ii++) {
         hit = (CODA_HIT_INFO *) &event->hits[(roc-1)][(slot-1)][ii];
         /* printf("%d %d %d %d\n",hit->crate_id,hit->slot_id,
          *        hit->chan_id,hit->module_id);
          */
         if ( (roc == hit->crate_id) && (slot == hit->slot_id) &&
              (chan == hit->chan_id) && (hit->module_id == F1TDC48) &&
              (hit->module_mode == mode) )
         {
            /* printf("f1tdc48_write_data: DEBUG: found hit for crate %d"
             *        " slot %d chan %d.\n",roc,slot,chan);
             */
            
            chit[jj] = hit;
            jj++;
            if (jj >= F1TDC48_MAX_HITS) {
               printf("f1tdc48_write_data: ERROR: Too many hits for channel\n");
               jj--;
					break;
            }
         }
      }
      /* printf("write hit data %d\n",jj); */
      for (ii=0; ii < jj; ii++) {
         F1TDC48_F1_DATA(cdata,chip,chan_on_chip,chit[ii]->hdata[0]);
      }
      
      /* Check for outputing chip trailer */
#if 0
      if (chan == 47) {
         F1TDC48_F1_TRAILER(cdata,(eventNum&0x3f),(ts&0x1ff),chan);
      }
#endif
      
   }
   
   nwords = dabufp - start;
   if ((nwords%2) == 0) {
      F1TDC48_FILLER(slot);
      nwords += 2;
      F1TDC48_BL_TRAILER(slot,nwords);
   } else {
      nwords += 1;
      F1TDC48_BL_TRAILER(slot,nwords);
   }
   
   if (nwords%4 != 0) {
      F1TDC48_FILLER(slot);
      F1TDC48_FILLER(slot);
      nwords += 2;
   }
   
   return nwords;
}


/* CAEN 1290 Hi Res TDC Parameters */
#define CAEN1290_MAX_CHAN      32
#define CAEN1290_MAX_HITS       8
#define CAEN1290_MAX_CHIPS      4


/* CAEN 1290 TDC Macros */
#define CAEN1290_BL_HEADER(slot,evt) {*dabufp++ =  0x40000000 | (evt << 5) | (slot); }
#define CAEN1290_BL_TRAILER(slot,nwords,status)   {*dabufp++ =  0x80000000 | ((status&0x7) << 24) | (nwords << 5) | (slot); }

#define CAEN1290_TDC_HEADER(chip,evt,bid)  {*dabufp++ = 0x08000000 | ((chip&0x3) << 24) | ((evt&0xfff) << 12) | (bid&0xfff); }
#define CAEN1290_TDC_DATA(edge,chan,time)  {*dabufp++ = 0x00000000 | ((edge&1) << 26) | ((chan&0x1f) << 21) | (time&0x1fffff); }
#define CAEN1290_TDC_TRAILER(chip,evt,nwords)  {*dabufp++ = 0x18000000 | ((chip&0x3) << 24) | ((evt&0xfff) << 12) | (nwords&0xfff); }

#define CAEN1290_TDC_ERROR(chip,eflags) { *dabufp++ = 0x20000000 | ((chip&0x3) << 24) | (eflags&0x7fff); }

#define CAEN1290_FILLER {*dabufp++ = 0xC0000000; }


int
caen1290_write_data (CODA_EVENT_INFO *event, int roc, int slot, int mode)
{
   int ii, jj, chan, hcnt, nwords=0, wcnt=0;
   //uint64_t tsdiv;
   uint32_t chip, stat, edge = 0;
   uint32_t  eventNum;
   // uint64_t  timestamp;
   unsigned int *start = dabufp;
   CODA_HIT_INFO *hit;
   CODA_HIT_INFO *chit[MAX_HITS_PER_CHAN];
   
   eventNum  = (event->eventid)&0xffffffff;
   // timestamp = (event->trigtime);
   hcnt      =  event->hcount[(roc-1)][(slot-1)];
   
   
   /* Set Status to 0 for now */
   stat = 0;
   
   
   CAEN1290_BL_HEADER(slot,eventNum);
   
   /*  printf("caen1290_write_data: %d hits checking slot %d"
    *         " roc %d.\n",hcnt, slot, roc);
    */
   
   /* Loop over all channels */
   chip = 0;
   for (chan=0; chan < CAEN1290_MAX_CHAN; chan++) {
      
      /* Check for outputing Chip Header */
      if ((chan%8) == 0) {
         CAEN1290_TDC_HEADER(chip, eventNum, 0);
         wcnt=0;
         chip++;
      }
      
      /* check for all hits for this channel */
      jj=0;
      for (ii=0; ii < hcnt; ii++) {
         hit = (CODA_HIT_INFO *) &event->hits[(roc-1)][(slot-1)][ii];
         /* printf("%d %d %d %d\n",hit->crate_id,hit->slot_id,
          *        hit->chan_id,hit->module_id);
          */
         if ( (roc == hit->crate_id) && (slot == hit->slot_id) &&
              (chan == hit->chan_id) && (hit->module_id == CAEN1290) &&
              (hit->module_mode == mode) )
         {
            /*printf("caen1290_write_data: DEBUG: found hit for crate %d"
             *       " slot %d chan %d.\n",roc,slot,chan);
             */
            
            chit[jj] = hit;
            jj++;
            wcnt++;
            if (jj >= CAEN1290_MAX_HITS) {
               printf("caen1290_write_data: ERROR: Too many hits for channel\n");
               return 0;
            }
         }
      }
      /* printf("write hit data %d\n",jj); */
      for (ii=0; ii < jj; ii++) {
         CAEN1290_TDC_DATA(edge,chan,chit[ii]->hdata[0]);
      }
      
      /* Check for outputing chip trailer */
      if (((chan+1)%8) == 0) {
         CAEN1290_TDC_TRAILER(chip, eventNum, wcnt);
      }
      
   }
   
   CAEN1290_BL_TRAILER(slot,nwords,stat);
   
   nwords = dabufp - start;
   
   while ((nwords%4) != 0) {
      CAEN1290_FILLER;
      nwords++;
   }
   
   return nwords;
}


#define MAX_CONFIG_PARS 20

void WriteDAQconfigBank(CODA_CRATE_MAP *crate, int roc)
{

	// Get type of digitization modules in this crate. Assume all are the same
	enum type_id_t modtype = UNKNOWN;
	int i;
	for(i=0; i<MAX_SLOTS; i++){
		int mymodtype = crate->module_map[i];
		switch(mymodtype){
			case FADC250:
			case FADC125:
			case F1TDC32:
			case F1TDC48:
			case CAEN1190:
			case CAEN1290:
				modtype = mymodtype;
				break;
			default:
				break;
		}
		if(modtype != UNKNOWN) break;
		
		if(mymodtype>TID && mymodtype<TD){
			modtype = crate->module_map[i];
			break;
		}
	}
	
	// Write parameters based on module type
	uint32_t Npar = 0;
	uint32_t partype[MAX_CONFIG_PARS];
	uint32_t parvalue[MAX_CONFIG_PARS];
	switch(modtype){
		case FADC250:
			partype[Npar] = kPARAM250_NSA_NSB;   parvalue[Npar++] = FADC250_WINDOW_WIDTH;
			//partype[Npar] = kPARAM250_NPED;      parvalue[Npar++] = 4; // hardwired in 1st gen. firmware
			partype[Npar] = kPARAM250_NPED;      parvalue[Npar++] = 1; // we save only one pedestal channel
			break;
		case FADC125:
			partype[Npar] = kPARAM125_NSA_NSB;   parvalue[Npar++] = FADC125_WINDOW_WIDTH;
			partype[Npar] = kPARAM125_WINWIDTH;  parvalue[Npar++] = FADC125_WINDOW_WIDTH;  // this is currently unused in offline code?
			//partype[Npar] = kPARAM125_NPED;      parvalue[Npar++] = 4; // hardwired in 1st gen. firmware
			partype[Npar] = kPARAM125_NPED;      parvalue[Npar++] = 1; // we save only one pedestal channel
			break;
		case F1TDC32:
		case F1TDC48:
			partype[Npar] = kPARAMF1_REFCNT;     parvalue[Npar++] = 0; // Need to get a real number here
			partype[Npar] = kPARAMF1_HSDIV;      parvalue[Npar++] = 0; // Need to get a real number here
			break;
		case CAEN1290:
			partype[Npar] = kPARAMCAEN1290_WINOFFSET; parvalue[Npar++] = 0; // This is just a placeholder
		default:
			break;
	}

	// If no parameters to write, then don't write bank at all
	if(Npar == 0) return;

	// Open data bank of type DAQConfig=0x55
	DATA_BANK_OPEN(0,0x55,1);

	// write header word
	uint32_t slot_mask = crate->moduleMask;
	*dabufp++ = (Npar<<24) | slot_mask;

	// Write all parameters
	uint32_t ii;
	for(ii=0; ii<Npar; ii++){
		*dabufp++ = (partype[ii]<<16) | (parvalue[ii]);
	}
	
	// Close data bank (see note above regarding DATA_BANK_CLOSE)
	DATA_BANK_CLOSE;
}


