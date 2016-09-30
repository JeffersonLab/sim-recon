#include <xstream/config.h>
#include <fstream>

#if HAVE_LIBZ

#include <algorithm>
#include <string.h>
#include <string>
#include <cstring>

#include <xstream/z.h>
#include <xstream/except/z.h>
#include <stdexcept>

#include <stdio.h>
#include <zlib.h>
#include <arpa/inet.h>

#include <cassert>

#include "debug.h"

#define COMPRESSION_BLOCK_SIZE 32000

// The following two macros must always occur in pairs within a single
// block of code, otherwise it will not even compile. This is done on
// purpose, to reduce the risk of blunders with deadlocks. Please take
// the lock, do the operation, and then release the lock as quickly as
// possible. If your function needs to return between the MUTEX_LOCK and
// MUTEX_UNLOCK statements, use MUTEX_ESCAPE before the return statement.

#define MUTEX_LOCK \
   { \
      if (streambuf_mutex != 0) \
         pthread_mutex_lock(streambuf_mutex); \
      pthread_mutex_t *mutex_saved = streambuf_mutex; \
      streambuf_mutex = 0;

#define MUTEX_UNLOCK \
      streambuf_mutex = mutex_saved; \
      if (streambuf_mutex != 0) \
         pthread_mutex_unlock(streambuf_mutex); \
   }

#define MUTEX_ESCAPE \
      streambuf_mutex = mutex_saved; \
      if (streambuf_mutex != 0) \
         pthread_mutex_unlock(streambuf_mutex);

namespace xstream {
namespace z {

// define the standard header for a zlib stream that can be used
// to prime the inflate engine to start at an arbitrary block in
// an input stream

    static int z_header_length = 2;
    static unsigned char z_header[2] = {0x78, 0x9c};

    struct pimpl: public z_stream {};
    
    static const int eof = std::streambuf::traits_type::eof();

    static inline int flush_macro(const flush_kind f) {
        switch (f) {
            case no_sync:
                return Z_NO_FLUSH;
            case normal_sync:
                return Z_SYNC_FLUSH;
            case full_sync:
                return Z_FULL_FLUSH;
            case finish_sync:
                return Z_FINISH;
            case block_sync:
#ifndef Z_BLOCK
# define Z_BLOCK 5
#endif
                return Z_BLOCK;
            default:
                //should probably throw
                return Z_NO_FLUSH;
        }
    }

    const char* error_str(int err) {
        switch(err) {
            case Z_MEM_ERROR:
                return "out of memory";
            case Z_VERSION_ERROR:
                return "zlib version mismatch";
            case Z_DATA_ERROR:
                return "invalid or incomplete data";
            case Z_STREAM_ERROR:
                return "stream error";
            case Z_NEED_DICT:
                return "need dictionary";
            case Z_STREAM_END:
                return "stream end";
            case Z_BUF_ERROR:
                return "buffer error";
        }

        return "unknown error";
    }
    

    common::common(std::streambuf * sb)
    : xstream::common_buffer(sb), z_strm(0), block_start(0), block_offset(0),
      streambuf_mutex(0)
    {
        LOG("z::common");    

        z_strm = new pimpl;

        //initialize zlib structure
        z_strm->zalloc = Z_NULL;
        z_strm->zfree = Z_NULL;
        z_strm->opaque = Z_NULL;
        //buffers
        z_strm->avail_out = out.size;
        z_strm->next_out = reinterpret_cast < Bytef* >(out.buf);

        z_strm->avail_in = 0;
        z_strm->next_in = reinterpret_cast < Bytef* >(in.buf);

    }

    void common::grow_out (unsigned int factor) {

        const size_t taken = out.size - z_strm->avail_out;

        out.grow(factor);

        z_strm->next_out = reinterpret_cast < Bytef* >(out.buf + taken);
        z_strm->avail_out = out.size - taken;
    }

    unsigned long int common::input_count() const {
        return z_strm->total_in;
    }

    unsigned long int common::output_count() const {
        return z_strm->total_out;
    }

    unsigned long int common::checksum() const {
        return z_strm->adler;
    }

    common::~common() {
        LOG("z::~common");
        delete z_strm;
    }

    ostreambuf::ostreambuf (std::streambuf * sb)
    : common(sb), level(Z_DEFAULT_COMPRESSION) {
        LOG("z::ostreambuf without compression level");
        block_start = _sb->pubseekoff(0, std::ios_base::cur, std::ios_base::out);
        init();
    }

    ostreambuf::ostreambuf(std::streambuf *sb, int l)
    : common(sb), level (l) {
        LOG ("z::ostreambuf with compression level " << l);
        block_start = _sb->pubseekoff(0, std::ios_base::cur, std::ios_base::out);
        init();
    }

    void ostreambuf::raise_error(int err) {
        std::string what = error_str(err);

        LOG("z::ostreambuf::raise_error (" << err << ") = " << what);

        if (what.size() > 0) {
            throw compress_error(this, what);
        } else {
            throw compress_error(this);
        }
    }

    void ostreambuf::init() {
        LOG ("z::ostreambuf::init");

        if (Z_DEFAULT_COMPRESSION == level || (level <= 9 && level >= 1)) {
            int cret =::deflateInit(z_strm, level);
            if (Z_OK != cret) {
                LOG ("z::ostreambuf::init: error creating zstream " << cret);
                //XXX exception ins constructor
                raise_error(cret);
            }
            //initialize streambuf interface functions
            setp(in.buf, in.buf + in.size);
        } else {
            char str[256];
            sprintf(str, "invalid compression level %d", level);
            throw std::domain_error(str);
        }
    }

    ostreambuf::~ostreambuf() {
        LOG ("z::ostreambuf::~ostreambuf");
        //sync (write remaining data)
        flush(finish_sync);

        //sync underlying streambuf
        MUTEX_LOCK
        _sb->pubsync();
        MUTEX_UNLOCK

        if (0 != z_strm) {
            //XXX should I throw an exception in case of error?
            //remember this is a destructor
            //I should definitely LOG something
            int cret = ::deflateEnd(z_strm);
            if (Z_OK != cret){
                LOG("z::~ostreambuf error dealocating zstream");
            }
        }
    }

    int ostreambuf::sync () {
      LOG ("z::ostreambuf::sync");
      int ret;
      MUTEX_LOCK
      ret = flush(finish_sync);
      _sb->pubsync();
      MUTEX_UNLOCK
      return ret;
    }


    int ostreambuf::overflow(int c) {
        LOG ("z::ostreambuf::overflow(" << c << ")\t available=" << (available ()) << "\tEOF=" << eof);
        if (eof == c) {
            LOG ("\tEOF");
            flush(no_sync);
            return eof;
        } else {
            if (0 == available ()) {
                LOG ("\t have to flush :[]");
                flush(no_sync);
            }
            *pptr () = static_cast < char >(c);
            pbump (1);
        }
        return c;
    }

    std::streamsize ostreambuf::xsputn (const char *buffer, std::streamsize n) {
        LOG ("z::ostreambuf::xsputn(" << buffer << "," << n << ")");

        return flush(no_sync, buffer, n);
    }

    int ostreambuf::flush(flush_kind f, const char *appendbuf, int appendsize) {
        LOG ("z::ostreambuf::flush(" << f << ")");
        std::streamsize in_s = taken ();
        LOG ("\tinput_size=" << in_s);

        //set up compression engine input feed
        int written;
        if (in_s > 0) {
           z_strm->next_in = reinterpret_cast < Bytef* >(pbase());
           z_strm->avail_in = in_s;
           written = in_s;
        } else if (appendsize > 0) {
           z_strm->next_in = (Bytef*)appendbuf;
           z_strm->avail_in = appendsize;
           written = appendsize;
           appendsize = 0;
        } else {
           z_strm->next_in = reinterpret_cast < Bytef* >(pbase());
           z_strm->avail_in = 0;
           written = 0;
        }
        block_offset += written;
        if (block_offset > (std::streamoff)COMPRESSION_BLOCK_SIZE) {
            f = (f == no_sync)? finish_sync : f;
        }

        if (z_strm->avail_in + z_strm->total_in == 0)
           return 0;

        bool redo = false;
        bool reinit_deflator = false;

        do {
            int cret;
            redo = false;

            do {
                cret = ::deflate(z_strm, flush_macro(f));

                if (finish_sync == f && Z_OK == cret) {
                    grow_out();
                    continue;
                } else {
                    break;
                }
            } while (1);

            //error handling
            if (f == finish_sync) {
                if (Z_STREAM_END == cret) {
                    redo = false;
                    reinit_deflator = true;
                }
                else {
                    //serious error, throw exception
                    LOG ("\terror :" << cret);
                }
            } else if (Z_OK != cret) {
                LOG ("\terror deflating " << cret);
                //XXX throw exception here
                raise_error(cret);
            }

            if (f == finish_sync) { // only completed streams can be written
                std::streamsize count = out.size - z_strm->avail_out;
                if (count > 0) { // ignore empty blocks
                    LOG ("\twriting " << count << " bytes");
                    int size = htonl(count);
                    MUTEX_LOCK
                    const std::streamsize wrote = _sb->sputn((char*)&size, 4) +
                                                  _sb->sputn(out.buf, count);
                    if (wrote != count + 4) {
                        MUTEX_ESCAPE
                        LOG("\terror writting, only wrote " << wrote 
                            << " but asked for " << count);
                        raise_error(Z_STREAM_ERROR);
                    }
                    block_start = _sb->pubseekoff(0, std::ios_base::cur,
                                                     std::ios_base::out);
                    block_offset = 0;
                    MUTEX_UNLOCK
                }
                z_strm->next_out = reinterpret_cast < Bytef* >(out.buf);
                z_strm->avail_out = out.size;
            }

            if (0 == z_strm->avail_out) { // && 0 != z_strm->avail_in)
                LOG("\tavail_out=0 => redo");
                redo = true;
            }

            if (!redo && appendbuf && appendsize > 0) {
                z_strm->next_in = (Bytef*)appendbuf;
                z_strm->avail_in = appendsize;
                written += appendsize;
                appendsize = 0;
                redo = true;
            }
        }
        while (redo);
        assert (0 == z_strm->avail_in);

        if (reinit_deflator) {
            int cret;
            cret = ::deflateEnd(z_strm);
            if (Z_OK != cret) {
                LOG("\tERROR: deflateEnd returned " << cret);
                raise_error(cret);
            }
            z_strm->zalloc = Z_NULL;
            z_strm->zfree = Z_NULL;
            z_strm->opaque = Z_NULL;
            z_strm->avail_out = out.size;
            z_strm->next_out = reinterpret_cast < Bytef* >(out.buf);
            z_strm->avail_in = 0;
            z_strm->next_in = reinterpret_cast < Bytef* >(in.buf);
            cret =::deflateInit(z_strm, level);
            if (Z_OK != cret) {
                LOG("\tERROR: deflateInit returned " << cret);
                raise_error(cret);
            }
        }

        //reset buffer
        setp(in.buf, in.buf + in.size);
        return written;
    }

    /////////////////////
    // istream follows //
    /////////////////////

    istreambuf::istreambuf (std::streambuf *sb, int *left, unsigned int left_size)
    : common(sb), end(false), block_size(0), block_next(0), 
      new_block_start(0), new_block_offset(0),
      leftovers(0)
    {
        LOG ("z::istreambuf");

        memset(z_strm, 0, sizeof(*z_strm));
        int cret = ::inflateInit(z_strm);

        if (Z_OK != cret) {
            LOG ("\terror creating zstream " << cret);
            //XXX throw exception here
            raise_error(cret);
        }
        //initialize streambuf interface functions
        //first call will call uflow and this will set the buffer accordingly
        //no buffering
        setg(out.buf, out.buf, out.buf);
        block_start = _sb->pubseekoff(0, std::ios_base::cur, std::ios_base::in);

        if (left_size >= sizeof(leftovers_buf)) {
            leftovers = (leftovers_buf*)left;
        }
        else {
            LOG("\terror - insufficient space for leftovers buffer");
            raise_error(cret);
        }
    }

    void istreambuf::raise_error(int err) {
        std::string what = error_str(err);

        LOG("z::istreambuf::raise_error (" << err << ") = " << what);

        if (what.size() > 0) {
            throw decompress_error(this, what);
        } else {
            throw decompress_error(this);
        }
    }

    int istreambuf::underflow() {
        LOG("z:istreambuf::underflow");

        if (end) {
            LOG("\tend of stream (EOF)");
            //signal the stream has reached it's end
            return eof;
        }

        if (new_block_start > 0 || new_block_offset > 0) {
            if (block_start != new_block_start ||
                block_offset > new_block_offset ||
                block_size == 0)
            {
                z_strm->next_out = reinterpret_cast < Bytef* >(out.buf);
                z_strm->avail_out = 0;
                read_inflate();
            }
            while (block_offset < new_block_offset) {
                z_strm->next_out = reinterpret_cast < Bytef* >(out.buf);
                z_strm->avail_out = new_block_offset - block_offset;
                if (z_strm->avail_out > out.size) {
                    z_strm->avail_out = out.size;
                }
                block_offset += z_strm->avail_out;
                inflate();
            }
            new_block_start = 0;
            new_block_offset = 0;
        }

        z_strm->avail_out = out.size;
        z_strm->next_out = reinterpret_cast < Bytef* >(out.buf);

        if (0 < z_strm->avail_in) {
            LOG("\tdata in queue, inflating");
            inflate();
        }
        while (!end && z_strm->avail_out > 0) {
            read_inflate();
        }
        if (end && z_strm->avail_out > 0) {
            LOG("\tend of stream (EOF)");
            //signal the stream has reached it's end
            return eof;
        }
            
        //set streambuf pointers
        setg(out.buf, out.buf, reinterpret_cast <char*> (z_strm->next_out) );

        return int(out.buf[0]);
    }
    
    //read to buffer in place (apart from data already buffered)
    std::streamsize istreambuf::xsgetn(char *buffer, std::streamsize n) {
        LOG("z::istreambuf::xsgetn (" << n << ")");

        if (new_block_start > 0 || new_block_offset > 0) {
            if (block_start != new_block_start ||
                block_offset > new_block_offset ||
                block_size == 0)
            {
                z_strm->next_out = reinterpret_cast < Bytef* >(out.buf);
                z_strm->avail_out = 0;
                read_inflate();
                setg(out.buf, out.buf, out.buf);
            }
            else
            {
                std::streamsize available = egptr() - gptr();
                int waste = new_block_offset - block_offset;
                waste = (available < waste)? available : waste;
                if (waste > 0) {
                    gbump(waste);
                    block_offset += waste;
                }
            }
            while (block_offset < new_block_offset) {
                z_strm->next_out = reinterpret_cast < Bytef* >(out.buf);
                z_strm->avail_out = new_block_offset - block_offset;
                if (z_strm->avail_out > out.size) {
                    z_strm->avail_out = out.size;
                }
                block_offset += z_strm->avail_out;
                inflate();
            }
            new_block_start = 0;
            new_block_offset = 0;
        }

        //try to satisfy request from buffered input
        std::streamsize available = egptr() - gptr();
        int read = (available >= n)? n : available;
        if (read) {
            std::copy(gptr(), gptr() + read, buffer);
            gbump(read);
            block_offset += read;
        }

        //inflate the rest directly into the user's buffer
        if (read < n) {
            if (end) {
                LOG("\tend of stream (EOF)");
                //signal the stream has reached it's end
                return eof;
            }

            z_strm->next_out = reinterpret_cast < Bytef* >(buffer) + read;
            z_strm->avail_out = n - read;

            if (0 < z_strm->avail_in) {
                inflate();
            }
            while (!end && z_strm->avail_out > 0) {
                read_inflate();
            }
            if (end && z_strm->avail_out > 0) {
                LOG("\tend of stream (EOF)");
                //signal the stream has reached it's end
                return eof;
            }
            block_offset += n - read;
        }
        return n;
    }

    void istreambuf::read_inflate( const flush_kind f) {
        LOG("z::istreambuf::read_inflate " << f);
        bool reinit_inflator = false;
        int read;
        if (block_size < 0) { // stream has no blocksize markers
            MUTEX_LOCK
            if (new_block_start > 0) {
               _sb->pubseekoff(new_block_start, std::ios_base::beg,
                                                std::ios_base::in);
               new_block_start = 0;
               leftovers->len = 0;
               block_next = 0;
               end = false;
            }
            read = _sb->sgetn(in.buf, in.size);
            MUTEX_UNLOCK
        }
        else { // look for prefixed blocksize: leading byte = 0
            MUTEX_LOCK
            if (new_block_start > 0) {
               _sb->pubseekoff(new_block_start, std::ios_base::beg,
                                                std::ios_base::in);
               block_start = new_block_start;
               new_block_start = 0;
               leftovers->len = 0;
               block_next = 0;
               end = false;
            }
            else {
               block_start = _sb->pubseekoff(0, std::ios_base::cur,
                                                std::ios_base::in);
               block_start -= leftovers->len;
            }
            reinit_inflator = (block_next != block_start);
            read = leftovers->len;
            if (read < 4) {
                read += _sb->sgetn(leftovers->buf + read, 4 - read);
                if (read != 4) {
                    end = true;
                    MUTEX_ESCAPE
                    return;
                }
            }
            if (leftovers->buf[0] == 0) { // z blocks have prefixed blocksize
                int *size = (int*)leftovers->buf;
                block_size = ntohl(*size);
                read -= 4;
                if (reinit_inflator && read > 0) {
                    std::memcpy(in.buf, leftovers->buf + 4, read);
                    read += _sb->sgetn(in.buf + read, block_size - read);
                }
                else {
                    read = _sb->sgetn(in.buf, block_size - read);
                }
                leftovers->len = _sb->sgetn(leftovers->buf, 8);
                if (leftovers->len > 4) {
                    std::memcpy(in.buf + read, leftovers->buf + 4,
                                               leftovers->len - 4);
                    read += leftovers->len - 4;
                }
            }
            else { // z blocks are jammed together, no blocksize available
                read += _sb->sgetn(in.buf + read, in.size - read);
                leftovers->len = 0;
                block_size = -1;
            }
            block_next = _sb->pubseekoff(0, std::ios_base::cur,
                                            std::ios_base::in);
            block_next -= leftovers->len;
            MUTEX_UNLOCK
        }
        LOG("\tread " << read << " bytes");
        block_offset = 0;

        if (0 == read) {
            end = true;
            return;
        }

        const char* head = (const char*)z_header;
        if (reinit_inflator) {
            int cret = ::inflateEnd(z_strm);
            if (Z_OK != cret) {
                LOG ("\terror terminating zstream " << cret);
                //XXX throw exception here
                raise_error(cret);
            }
            z_strm->zalloc = Z_NULL;
            z_strm->zfree = Z_NULL;
            z_strm->opaque = Z_NULL;
            cret = ::inflateInit(z_strm);
            if (Z_OK != cret) {
                LOG ("\terror initializing zstream " << cret);
                //XXX throw exception here
                raise_error(cret);
            }
            if (strncmp(head, in.buf, z_header_length) != 0) {
                z_strm->avail_in = z_header_length;
                z_strm->next_in = reinterpret_cast < Bytef* >(z_header);
                inflate(f); // inject the z stream header
            }
        }
        z_strm->next_in = reinterpret_cast < Bytef* >(in.buf);
        z_strm->avail_in = read;
        inflate(f);
    }

    void istreambuf::inflate(const flush_kind f) {
        LOG("z::istreambuf::inflate " << f);

        int cret = ::inflate(z_strm, flush_macro(f));

        if (Z_STREAM_END == cret) {
            z_strm->avail_in = 0;
            block_next = 0;
        }
        else if (cret == Z_DATA_ERROR && z_strm->avail_in == 0) {
            // Ignore CRC errors at the end of stream because we may not have
            // started inflating at the beginning. We can rely on the CRC
            // checks that are present within each compressed block anyway.
            end = true;
        }
        else if (Z_OK != cret) {
            printf("z input stream crapping out, cret is %d\n", cret);
            LOG("\terror inflating: " << cret);
            //XXX throw exception
            raise_error(cret);
            //can try to salvage some more data with inflateSync (on some cases)
        }
    }

    istreambuf::~istreambuf() {
        LOG("z::~istreambuf");
        if (0 != z_strm) {
            //XXX should I throw an exception in case of error?
            //remember this is a destructor
            ::inflateEnd(z_strm);
        }
    }

}//namespace z
}//namespace xstream

#endif    //zlib
