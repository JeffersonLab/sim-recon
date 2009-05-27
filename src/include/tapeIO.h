/*
 * tapeIO.h
 *
*/

#ifndef TAPEIOH_INCLUDED

#define TAPEIO_OK       0
#define TAPEIO_ERROR  (-1)
#define TAPEIO_EOT    (-2)   /* end of data (End Of Tape) */

#define TAPEIO_STREAM   1    /* file descriptor is not a tape device      */
#define TAPEIO_FIXED    2    /* it's a tape device in fixed block mode    */
#define TAPEIO_VARIABLE 3    /* it's a tape device in variable block mode */
#define TAPEIO_PIPE     4    /* pipe or TCP socket */

int tape_devType(int fd); /* returns one of the TAPEIO_xxxx or (-1) if error */

int tape_write(int fd,const void*buffer,int len);
int tape_read(int fd,void*buffer,int len);
int tape_flush(int fd);
int tape_writeFM(int fd);

int tape_getBlockSize(int fd,unsigned long*minSize,unsigned long*maxSize,unsigned long*recSize);
int tape_getState(int fd);              /* returns bits defined in mtio.h: MT_BOT, etc... */
unsigned long tape_getPosition(int fd); /* returns current tape block number (does not work for jag tapes) */

int tape_setBlockSize(int fd,unsigned long newBlockSize);

int tape_rewind(int fd);
int tape_unload(int fd);
int tape_findFM(int fd,int count);  /* count > 0 seeks forward, count < 0 seeks backwards */

extern int tape_enablePipeIO;  /* enable the special handing of PIPEs */

#endif
/* end file */
