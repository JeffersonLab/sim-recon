/*
 * This sowftware was developed at Jefferson National Laboratory.
 * (c) 2017.
 */

#include "reader.h"
#include "hipoexceptions.h"
#include "record.h"

#include <cstdlib>
/**
 * HIPO namespace is used for the classes that read write
 * files and records.
 */
namespace hipo {
 /**
  * The constructor for reader, printWarning routine
  * will printout a warning message if the library
  * was not compiled with compression libraries LZ4 or GZIP
  */
  reader::reader(){
    printWarning();
  }
  /**
   * Default destructor. Does nothing
   */
  reader::~reader(){ }
  /**
   * Open file, if file stream is open, it is closed first.
   * At open time verification of file structure is performed.
   * If the signature does not match EVIO/HIPO template, the
   * file will be closed and warning message is printed.
   */
  void reader::open(const char *filename){

    if(inputStream.is_open()==true){
      inputStream.close();
    }

    inputStream.open(filename, std::ios::binary);
    if(inputStream.is_open()==false){
      printf("[ERROR] something went wrong with openning file : %s\n",
            filename);
      return;
    }
    readHeader();
    bool status = verifyFile();
    if(status==false){
      inputStream.close();
    }
    readRecordIndex();
}

/**
 * Reads the file header. The endiannes is determined for bytes
 * swap. The header structure will be filled with file parameters.
 */
 void reader::readHeader(){

    headerBuffer.resize(80);
    inputStream.read(&headerBuffer[0],80);
    header.uniqueid     = *(reinterpret_cast<int *>(&headerBuffer[0]));
    header.filenumber   = *(reinterpret_cast<int *>(&headerBuffer[4]));
    header.headerLength = *(reinterpret_cast<int *>(&headerBuffer[8]));
    header.recordCount  = *(reinterpret_cast<int *>(&headerBuffer[12]));

    header.indexArrayLength = *(reinterpret_cast<int *>(&headerBuffer[16]));
    int word_8 = *(reinterpret_cast<int *>(&headerBuffer[20]));

    header.userHeaderLength = *(reinterpret_cast<int *>(&headerBuffer[24]));
    header.magicNumber = *(reinterpret_cast<int *>(&headerBuffer[28]));
    header.userRegister = *(reinterpret_cast<long *>(&headerBuffer[32]));

    // If magic word is reversed, then the file was written in BIG_ENDIAN
    // format, the bytes have to be swapped
    if(header.magicNumber==0x0001dac0){
       printf(" THIS FILE IS BIG ENDIAN: SWAPPING BYTES\n");
       header.uniqueid  = __builtin_bswap32(header.uniqueid);
       header.filenumber = __builtin_bswap32(header.filenumber);
       header.headerLength = __builtin_bswap32(header.headerLength);
       header.recordCount  = __builtin_bswap32(header.recordCount);
       header.userHeaderLength = __builtin_bswap32(header.userHeaderLength);
       header.indexArrayLength = __builtin_bswap32(header.indexArrayLength);
       word_8 = __builtin_bswap32(word_8);
       header.userRegister = __builtin_bswap64(header.userRegister);
    }

    header.version = word_8&0x000000FF;
    header.bitInfo = (word_8>>8)&0x00FFFFFF;
    header.firstRecordPosition = 4*header.headerLength + header.userHeaderLength;
    //int *signature = reinterpret_cast<int *>(&headerBuffer[0]);
    //printf("signature = %X\n",(unsigned int) *signature);
    //std::cout << "signature = " << std::ios::hex << (*signature) << '\n';
}
/**
 * Verify if the file has a proper format, magic word is checked
 * to and endianness is determined. Byte swap is performed if neccessary.
 */
bool      reader::verifyFile(){
 return true;
}
/**
 * Checks to determine if the file is open.
*/
bool  reader::isOpen(){
  return inputStream.is_open();
}
/**
 * Reads records indicies, it hopes through file Reading
 * only header for each records and fills a vector with
 * record header descriptors, identifying the position,
 * and number of events in each record.
 * If it encounters mistake it will preserve all recovered
 * record information.
 */
void  reader::readRecordIndex(){

    inputStream.seekg(0,std::ios::end);
    long hipoFileSize = inputStream.tellg();
    long positionOffset = header.firstRecordPosition;
    recordIndex.clear();
    inputStream.seekg(positionOffset,std::ios::beg);
    std::vector<char> recheader(80);
    int icounter   = 0;
    int eventCount = 0;
    while(positionOffset+56<hipoFileSize){

        inputStream.read( (char *) &recheader[0],56);
        recordIndex_t recIndex;

        recIndex.recordPosition  = positionOffset;
        recIndex.recordLength    = *(reinterpret_cast<int *>(&recheader[0]));
        recIndex.recordEvents    = *(reinterpret_cast<int *>(&recheader[12]));
        int compressWord         = *(reinterpret_cast<int *>(&recheader[36]));
        int version              = *(reinterpret_cast<int *>(&recheader[20]));
        int magic_number         = *(reinterpret_cast<int *>(&recheader[28]));

        recIndex.recordDataLengthCompressed = compressWord&0x0FFFFFFF;
        //recIndex.compressionType            = (compressWord&0xF0000000)>>28;

        if(magic_number==0x0001dac0){
          recIndex.recordLength = __builtin_bswap32(recIndex.recordLength);
          recIndex.recordEvents = __builtin_bswap32(recIndex.recordEvents);
          compressWord          = __builtin_bswap32(compressWord);
          version               = __builtin_bswap32(version);
        }
        version = (version&0x000000FF);
        if(version!=6) {
          printf(" version error : %d\n",version ); break;
        }
        if(magic_number!=0xc0da0100&&magic_number!=0x0001dac0) {
          printf("magic number error : %X\n", (unsigned int) magic_number);
          break;
        }

        //printf("version = %d , magic number = %X\n",version,(unsigned int) magic_number);
        positionOffset += recIndex.recordLength*4;
        inputStream.seekg(positionOffset,std::ios::beg);
        recordIndex.push_back(recIndex);
        //positionOffset +=
        icounter++;
        /*printf(" next position (%4d) : %16ld file size %ld events = %d\n",
          icounter,positionOffset,hipoFileSize, recIndex.recordEvents);*/
    }
    printf("total records = %d  index array = %d\n",icounter, (unsigned int) recordIndex.size());
}

void  reader::readRecord(hipo::record &record, int index)
{
    long position = recordIndex[index].recordPosition;
    record.readRecord(inputStream,position,0);
}
void  reader::readRecord(int index){
   hipo::record rec;
   long position = recordIndex[index].recordPosition;
   rec.readRecord(inputStream,position,0);
}
/**
 * Print the file information on the screen.
 */
void reader::showInfo(){
    printf("FILE: \n");
    printf(" %18s : %X\n","uniqueid",(unsigned int) header.uniqueid);
    printf(" %18s : %d\n","file #", header.filenumber);
    printf(" %18s : %d\n","header length", header.headerLength);
    printf(" %18s : %d\n","record count", header.recordCount);
    printf(" %18s : %d\n","index length", header.indexArrayLength);
    printf(" %18s : %d\n","version", header.version);
    printf(" %18s : %X\n","bit info", (unsigned int) header.bitInfo);
    printf(" %18s : %d\n","user header", header.userHeaderLength);
    printf(" %18s : %X\n","magic number", (unsigned int) header.magicNumber);
    printf(" %18s : %ld\n","first record", header.firstRecordPosition);
    if(recordIndex.size()<1){
      printf(" there are no records in the file : %d\n", inputStream.is_open());
      return;
    }
    long  recordPosition = 1000;//recordIndex[recordIndex.size()-1].recordPosition;
    float sizePos = recordPosition/1024.0/1024.0;
    printf("-------------------------------------\n");
    //printf("     signature : %X\n", (unsigned int) getSignature());
    //printf(" header Length : %d bytes\n", (unsigned int) getHeaderLength());
    printf("   file Length : %.2f MB\n", sizePos);
    printf("-------------------------------------\n");
}

int   reader::getRecordCount(){
   return recordIndex.size();
}
/**
 * Print warning if the library was not compiled with LZ4 library.
 * When this message appears, the compressed files will be unreadable.
 */
void reader::printWarning(){
    #ifndef __LZ4__
      std::cout << "******************************************************" << std::endl;
      std::cout << "* WARNING:                                           *" << std::endl;
      std::cout << "*   This library war compiled without LZ4 support.   *" << std::endl;
      std::cout << "*   Reading and writing compressed buffers will not  *" << std::endl;
      std::cout << "*   work. However un-compressed file I/O will work.  *" << std::endl;
      std::cout << "******************************************************" << std::endl;
    #endif
  }
}
