/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   hiporecord.h
 * Author: gavalian
 *
 * Created on April 11, 2017, 4:47 PM
 */

#ifndef HIPORECORD_H
#define HIPORECORD_H

#include <iostream>
#include <vector>
#include <string>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

namespace hipo {

    typedef struct {
        int signatureString; // 1) identifier string is HREC (int = 0x43455248
        int recordLength; // 2) TOTAL Length of the RECORD, includes INDEX array
        int recordDataLength; // 3) Length of the DATA uncompressed
        int recordDataLengthCompressed; // 4) compressed length of the DATA buffer
        int numberOfEvents ; // 5) number of event, data buckets in DATA buffer
        int headerLength ; // 6) Length of the buffer represengin HEADER for the record
        int indexDataLength ; // 7) Length of the index buffer (in bytes)
        int userHeaderLength; // user header length in bytes
        int userHeaderLengthPadding; // the padding added to user header Length
        int bitInfo;
        int compressionType;
        int compressedLengthPadding;
        int dataEndianness;
    } recordHeader_t;

    class data {
      private:
        const char  *data_ptr;
        int          data_size;
        int          data_endianness;
        int          data_offset;

      public:
        data(){ data_ptr = NULL; data_size = 0;}
        ~data(){ }

        void setDataPtr(const char *__ptr){ data_ptr = __ptr;}
        void setDataSize(int __size){ data_size = __size;}
        void setDataOffset(int __offset) { data_offset = __offset;}
        void setDataEndianness(int __endianness) { data_endianness = __endianness;}

        const uint32_t   *getEvioPtr(){ return reinterpret_cast<const uint32_t *>(data_ptr);}
        int         getEvioSize(){ return (int) data_size/4 ;}
        const char *getDataPtr(){ return data_ptr;}
        int         getDataSize(){ return data_size;}
        int         getDataEndianness(){ return data_endianness;}
        int         getDataOffset(){ return data_offset;}
    };

    class record {

      private:

        //std::vector< std::vector<char> > eventBuffer;
        std::vector<char>  recordHeaderBuffer;
        recordHeader_t     recordHeader;

        std::vector<char>  recordBuffer;

        char *getUncompressed(const char *data, int dataLength, int dataLengthUncompressed);
        int   getUncompressed(const char *data, char *dest, int dataLength, int dataLengthUncompressed);
        void  showBuffer(const char *data, int wrapping, int maxsize);

    public:

        record();
        ~record();

        void  readRecord(std::ifstream &stream, long position, int dataOffset);
        int   getEventCount();
        void  readEvent( std::vector<char> &vec, int index);
        void  getData(   hipo::data &data, int index);
    };
}
#endif /* HIPORECORD_H */
