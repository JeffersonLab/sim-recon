/*! \file xstream/digest.h
 *
 * \brief C++ objects to calculate digests of data
 */

#ifndef __XSTREAM_DIGEST_H
#define __XSTREAM_DIGEST_H

#include <xstream/config.h>
#include <xstream/common.h>

#include <streambuf>

#if HAVE_INTTYPES_H
#    include<inttypes.h>
#else
#    if HAVE_STDINT_H
#        include <stdint.h>
#    endif
#    error "I need inttypes.h or stdint.h to exist"    
#endif

#include <iosfwd>


namespace xstream{
/*!
 * \brief digest objects
 */
namespace digest{

/*
 * \brief streambuf class for stream digests
 *
 */
class stream: public xstream::ostreambuf
{
    protected: // changed to protected to avoif clang compiler error when
               // sync() is called from class "common" below. 10/23/2013 DL
        /*!
         * \brief update digest with as much data as possible (overloaded from streambuf)
         *
         * */
        int sync();
    private:   // change back to private (see note above) 10/23/2013 DL
        /*!
         * \brief write a character that surpasses buffer end (overloaded from streambuf)
         * 
         */
        int overflow(int c);

        /*!
         * \brief add an entire buffer to digest calculation (overloaded from streambuf)
         *
         */
        std::streamsize xsputn(const char *buffer, std::streamsize n);

    protected:
        xstream::buffer buf; /*!<buffer data to calculate digest */
        uint64_t length; /*!< number of bytes read so far */
        
        /*!
         * \brief default constructor
         *
         * allocates the buffer
         *
         * \parameter len length of buffer
         *
         * */
        stream(size_t len);

        /*!
         * \brief updates the digest
         * must be inplemented by classes that implement this interface
         *
         */
        virtual void calculate_digest()=0;

        /*!
         * \brief resets digest to it's initial value
         *
         */
        virtual void reset_digest()=0;

        /*! destructor
         *
         * frees the buffer
         *
         * */

        ~stream();
};

/*!
 * \brief Digest base class
 *
 * general interface for digest functions
 *
 */

template <typename digest_type>
class common: public stream{
    public:

        /*
         * \brief constructor
         *
         * \param s size of the buffer
         *
         */
        common(size_t s)
            :stream(s)
        {}

        /*!
         * \brief return the digest value
         *
         */
        virtual digest_type digest()=0;

        /*!
         * \brief resets the digest calculation
         *
         * returns the digest of data so far and for future calculations only considers data entered from now on
         *
         */

        digest_type reset()
        {
            sync();
            digest_type d = digest();
            reset_digest();
            return d;
        }
};
        

#if HAVE_LIBZ

/*!
 * \internal
 * \brief zlib's digest classes base
 *
 * zlib implements \c adler32 and \c crc32 digests
 *
 */

class z_common : public common<unsigned long int> {
    protected:
        unsigned long int _digest; /*!< digest value */

        virtual void reset_digest();

        z_common();

    public:
        virtual unsigned long int digest();
};

/*!
 * \brief adler32 digest class
 *
 */

class adler32 : public z_common {
    private:
        void calculate_digest();
};

/*!
 * \brief crc32 digest class
 *
 */

class crc32 : public z_common {
    private:
        void calculate_digest();
};

/*!
 * \brief base class for digest that work on fixed sized chunks of data
 */

class block_stream: public stream{
    private:
        const size_t chunk_size; /*!< size of the chunk used to calculate the digest */

        /*!
         * \brief update digest with as much data as possible (overloaded from streambuf)
         *
         * */
        int sync();

        /*!
         * \brief write a character that surpasses buffer end (overloaded from streambuf)
         * 
         */
        int overflow(int c);

        /*!
         * \brief add an entire buffer to digest calculation (overloaded from streambuf)
         *
         */
        std::streamsize xsputn(const char *buffer, std::streamsize n);

    public:

        /*!
         * \brief constructor
         *
         * \param chunk size in bytes of the chunk \c digest uses to update it's value 
         * \param blocks number of blocks of size chunk that fit in the buffer size
         *
         * */
        block_stream(size_t chunk, unsigned int blocks);
        
};


/*!
 * \brief md5 digest class
 *
 * return digest as 4 32bit words ( ABCD as defined in the RFC 1321, 128bits total) un a simple wrapping structure
 *
 */

class md5 : public block_stream {

    public:

        /*
         * \brief result of an \c md5 digest
         *
         */

        struct result {
            uint32_t a;
            uint32_t b;
            uint32_t c;
            uint32_t d;

            bool operator==(const result& r) const{
                return (a==r.a && b==r.b && c==r.c && d==r.d);
            }
        };

        /*!
         * \brief returns the digest of the data as an \c result structure
         * 
         */
        result digest();

        /*!
         * \brief resets digest value
         *
         * returns the digest of data so far and for future calculations only considers data entered from now on
         */

        result reset();
    
        /*!
         * \brief default constructor
         *
         */
        md5();

    private:
        result result;
        virtual void reset_digest();

        virtual void calculate_digest();

            
};
            
#endif //have zlib

/*
 * \brief dump md5 value to a stream 
 */

std::ostream& operator<<(std::ostream& o, const struct md5::result& m);

}//namespace digest
}//namespace xstream


/*!
 * \example adler.cpp
 *
 * Calculates adler32 checksum of data
 *
 */

/*!
 * \example crc.cpp
 *
 * Calculates crc32 checksum of data
 *
 */

/*!
 * \example md5.cpp
 *
 * Calculates md5 checksum of data
 *
 */
#endif
