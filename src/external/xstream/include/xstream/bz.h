/*! \file xstream/bz.h
 *
 * \brief C++ streambuf interface to read and write bzip2 streams
 */

#ifndef __XSTREAM_BZ_H
#define __XSTREAM_BZ_H

#include <xstream/config.h>

#if HAVE_LIBBZ2

#include <xstream/common.h>
#include <streambuf>

namespace xstream{
/*!
 * \brief bzip2 compression/decompression objects
 *
 */
namespace bz{

// forward declaration to avoid including any bzlib headers
struct pimpl;

/*!
 * \brief flush methods for compress
 *
 */
enum flush_kind {
    no_sync,    /*!< flush the minimum possible data */
    full_sync,  /*!< writes current "compression block"
    corruption at a previous position, data can be decompressed from this point onward */
    finish_sync /*!< write all data so that the stream can be closed */
};

/*!
 * \brief "private" class to factor some common code related with aquisition and release of bzlib objects
 */
class common: public xstream::common_buffer {
    protected:
        pimpl* z_strm; /*!< bzlib stream "object" */

        std::streampos block_start;
        std::streamoff block_offset;
        pthread_mutex_t *streambuf_mutex;

        /*!
         * \brief construct using a streambuf
         */
        common(std::streambuf* sb);

    public:
        unsigned long int output_count() const;

        /*!
         * \brief number of bytes of input so far
         *
         */
        unsigned long int input_count() const;

        /*!
         * \brief deallocates the buffers
         *
         */
        ~common();

        std::streampos get_block_start() {
            return block_start;
        }
        std::streamoff get_block_offset() {
            return block_offset;
        }
        pthread_mutex_t *get_streambuf_mutex() {
            return streambuf_mutex;
        }
        void set_streambuf_mutex(pthread_mutex_t *mutex) {
            streambuf_mutex = mutex;
        }
};


/*!
 * \brief ouput Bzip2 stream class
 *
 * Just a thin wrapper for streambuf to compress data
 * doesn't use bzlib's fstream like interface, so no unnecessary buffer layer added
 *
 */
class ostreambuf: public common, public xstream::ostreambuf {
    private:
        int level; /*!< compression level */

        /*!
         * \brief inspect bzlib error status and raise exception in case of error
         *
         */
        void raise_error(int err);

        void init(void);

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

        /*!
         * \brief write an entire buffer (overloaded from streambuf)
         *
         */
        std::streamsize xsputn(const char *buffer, std::streamsize n);
        
        /*!
         * \brief fine tuned flushing of stream
         *
         * \param f if determines the kind of flush to do.
         * if full_sync closes the current compression block
         *
         */
        int flush(flush_kind f=no_sync, const char *appendbuf=0, int appendsize=0);

    public:
        /*!
         * \brief construct using a streambuf to write to
         */
        ostreambuf(std::streambuf* sb);

        /*! \brief construct specifying the compression level
         *
         * \param sb streambuf to use
         * \param level compression level 
         *
         * \note level should be between 1(worst compression) and 9 (best compression)
         */
        ostreambuf(std::streambuf* sb, int level);

        /*!
         * \brief closes the bzlib stream
         *
         */
        ~ostreambuf();

        std::streambuf *get_streambuf() {
            return _sb;
        }
};

/*!
 * \brief input Bzip2 stream class
 *
 * Just a thin wrapper for streambuf to decompress bzip2 data streams
 * doesn't use bzlib's fstream like interface, so no unnecessary buffer layer added
 *
 */

class istreambuf: public common, public std::streambuf {
    private:
        
        bool end; /*!<signals if stream has reached the end */

        std::streamsize block_size;
        std::streampos block_next;
        std::streamoff new_block_start;
        unsigned int new_block_offset;
        typedef struct {
            int len;
            char buf[64];
        } leftovers_buf;
        leftovers_buf *leftovers;

        /*!
         * \brief inspect bzlib error status and raise exception in case of error
         *
         */
        void raise_error(int err);

        /*!
         * \brief requests that input buffer be reloaded (overloaded from streambuf)
         */
        int underflow();

        /*!
         * \brief encapsulates call to BZ2_bzDecompress and does error handling
         * */
        void decompress();

        /*!
         * \brief reads data and decompresses it
         */
        void read_decompress();

        /*!
         * \brief reads \c n characters to \c buffer (overloaded from streambuf)
         *
         */
        std::streamsize xsgetn(char *buffer, std::streamsize n);

    public:
        /*!
         * \brief  using a streambuf to read from
         *
         */
        istreambuf(std::streambuf* sb, int* left=0, unsigned int left_size=0);

        /*!
         * \brief closes the bzlib stream
         *
         */
        ~istreambuf();

        std::streambuf *get_streambuf() {
            return _sb;
        }
        std::streamsize get_block_size() {
            return block_size;
        }
        void set_new_position(std::streamoff start, unsigned int offset) {
           new_block_start = start;
           new_block_offset = offset;
        }
};

}//namespace bz
}//namespace xstream

/*!
 * \example bz_decompress.cpp
 *
 * shows how to use the xstream::bz::istreambuf to decompress an ordinary istream with data compressed in bzlib or gzip format
 *
 * it can decompress from standard input to standard output or from file to file
 *
 */

/*!
 * \example bz_compress.cpp
 *
 * shows how to use the xstream::bz::ostreambuf to compress data into bzlib format
 *
 * it can decompress from standard input to standard output or from file to file
 *
 */

#endif //have libbz2

#endif
