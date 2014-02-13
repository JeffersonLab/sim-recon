/*
 * dataIO.h
 *
*/

#ifndef _DATAIO_INCLUDED
#define _DATAIO_INCLUDED

static const char sccsidDataIOH[] = "@(#)dataIO.h\t5.2\tCreated 7/27/97 18:54:32, \tcompiled "__DATE__;

#define DATAIO_EOF      0               /* file mark or end of tape              */
#define DATAIO_OK       1               /* no errors                             */
#define DATAIO_EOT      2               /* end of tape                           */
#define DATAIO_ERROR  (-1)              /* read/write error (message->stderr)    */
#define DATAIO_BADCRC (-2)              /* CRC mismatch error: data is damaged   */

int data_write(int fd,const void* event);         /* automatically regenerates the CRC word        */

int data_read(int fd,void* buffer,int bufsize);
int data_read_alloc(int fd,void **buffer);         /* memory pointed by *buffer should be free()ed  */

int data_flush(int fd);   /* Flush internal buffers: should always be called before a close()       */

int data_writeFM(int fd);                          /* Write a file mark on tape                     */

int data_findFM(int fd,int count); /* find a File Mark. count>0 seeks forward, count<0 seeks backwards */

/*
 * data_SeekTape: Position the data tape at the requested run number.
 *
 * Input arguments:     fd - file descriptor of the tape device (as returned by 'open' or 'fileno')
 *                   runNo - requested run number
 *
 * Return value:    0  if run was found.
 *                (-1) if run was not found
 *
 * Output arguments:  tapeNo - if a tape label is encountered while
 *                             seeking the requested run,
 *                             'tapeNo' will be set to the run number from it.
 *                             Otherwise the value is left unchanged.
*/

int data_SeekTape(int fd,int runNo,int *tapeNo);

#endif
/* endfile */
