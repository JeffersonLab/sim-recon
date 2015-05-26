/*! \file xstream/base64.h
 *
 * \brief C++ streambuf interface to encode/decode base64 data
 */

#ifndef __XSTREAM_BASE64_H
#define __XSTREAM_BASE64_H

#include <xstream/config.h>

#include <xstream/common.h>
#include <streambuf>

namespace xstream{
/*!
 * \brief base64 encoding/decoding objects
 */
namespace base64{


/*!
 * \brief Base64 encode stream class
 *
 * encodes data to base64 and writes it to another streambuf
 *
 * \todo no \c xsputn support yet and \c flush doesn't flush partial data,
 * for example suppose you wrote 2 bytes of a 3 byte chunk, it's possible to know
 * the first 2 characters of the encoded data, but these are not flushed
 *
 */
class ostreambuf: public xstream::ostreambuf
{
    private:

        std::streambuf* _sb;
        char delim;
        unsigned int delim_w;
        unsigned int col;

        char buf[3]; /*!< buffer to store non encoded data */

        /*!
         * \brief flush as much data as possible (overloaded from streambuf)
         *
         * */
        int sync();

        /*!
         * \brief write a character that surpasses buffer end (overloaded from streambuf)
         * 
         */
        int overflow(int c);

#if 0
        //XXX implement later
        /*!
         * \brief write an entire buffer (overloaded from streambuf)
         *
         */
        std::streamsize xsputn(const char *buffer, std::streamsize n);
#endif

        /*!
         * \brief reset input buffer
         *
         */
        void reset();

        /*!
         * \brief Takes care of inserting delimiters every \c delim_w characters
         */

        int write(const char* buf, size_t len);
        
    public:
        /*!
         * \brief construct using a streambuf
         *
         * \param sb streambuf where to write the encoded data
         * \param width width in bytes of lines (=0 for unlimited lines)
         * \param delimiter char to delimit the lines
         *
         * the same parameters need to be used when decoding or it will throw an exception indicating invalid data
         * \c width=76 and \c delim=newline are default because they are the values in the \c RFC
         */
        ostreambuf(std::streambuf* sb, unsigned int width=76, char delimiter='\n');

        /*!
         * \brief closes the base64 stream
         *
         */
        ~ostreambuf();

};

/*!
 * \brief Base64 decode stream class
 *
 * decodes data encoded in base64 and makes it available for reading.
 * Only 4 bytes buffered at a time, and this is not configurable
 *
 * \todo xsgetn support
 * 
 */

class istreambuf: public std::streambuf{
    private:
        std::streambuf* _sb; /*!< stream to read from    */
        bool end; /*!<signals if stream has reached the end */
        char delim;
        unsigned int delim_w;
        unsigned int col;
        char buf[3]; /*!< buffer to store decoded data */

        /*!
         * \brief requests that input buffer be reloaded (overloaded from streambuf)
         */
        int underflow();

#if 0 //implemente later
        /*!
         * \brief reads \c n characters to \c buffer (overloaded from streambuf)
         *
         */

        std::streamsize xsgetn(char *buffer, std::streamsize n);
#endif

    public:
        /*!
         * \brief construct using a streambuf
         *
         * \param sb streambuf where to read the encoded data from
         * \param width width in bytes of lines (=0 for unlimited lines)
         * \param delimiter char to delimit the lines
         *
         * the parameters need to be the same specified at encoding or it will throw an exception indicating invalid data
         * \c width=76 and \c delim=newline are default because they are the values in the \c RFC
         *
         * if data is not delimited exactly at \c width characters with a newline, decoding will throw an exception.
         * RFC specifies that whitespace should be ignored but this library is very strict when decoding, interpreting spurious whitespace as an error.
         *
         */
        istreambuf(std::streambuf* sb, unsigned int width=76, char delimiter='\n' );

        /*!
         * \brief closes the base64 stream
         *
         */
        ~istreambuf();
};

/*!
 * \example b64_encode.cpp
 *
 * shows how to use the xstream::base64::ostreambuf to encode an ordinary ostream into base64
 *
 * it can encode from standard input to standard output or from file to file
 *
 */

/*!
 * \example b64_decode.cpp
 *
 * shows how to use the xstream::base64::istreambuf to decode base64 encoded data from an istream with encoded data
 *
 * it can decode from standard input to standard output or from file to file
 *
 */

}//namespace base64
}//namespace xstream

#endif
