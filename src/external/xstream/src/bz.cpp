#include <xstream/config.h>
#include <fstream>

#if HAVE_LIBBZ2

#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <cstring>
#include <algorithm>
#include <cassert>

#include <xstream/bz.h>
#include <xstream/except/bz.h>

#include <bzlib.h>
#include <arpa/inet.h>

#include "debug.h"

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
namespace bz {

// define a set of compressed stream headers that can be used
// to remove arbitrary bit-alignment offsets in an input stream

    static const int bz_header_length[8] = {35,36,34,44,40,38,40,43};
    static const unsigned char bz_header[8][48] = {
           { 0x42, 0x5a, 0x68, 0x39, 0x31, 0x41, 0x59, 0x26,
             0x53, 0x59, 0x86, 0xad, 0x3f, 0xf0, 0x00, 0x00,
             0x02, 0x48, 0x00, 0x04, 0x00, 0x20, 0x00, 0x20,
             0x00, 0x21, 0x00, 0x82, 0x0b, 0x31, 0x41, 0x59,
             0x26, 0x53, 0x59 },
           { 0x42, 0x5a, 0x68, 0x39, 0x31, 0x41, 0x59, 0x26,
             0x53, 0x59, 0x7b, 0x68, 0x3e, 0x4a, 0x00, 0x00,
             0x01, 0x88, 0x00, 0x0f, 0xc0, 0x20, 0x00, 0x21,
             0x80, 0x0c, 0x01, 0x37, 0xa4, 0xbb, 0x98, 0xa0,
             0xac, 0x93, 0x29, 0xac },
           { 0x42, 0x5a, 0x68, 0x39, 0x31, 0x41, 0x59, 0x26,
             0x53, 0x59, 0x83, 0x69, 0xfc, 0x04, 0x00, 0x00,
             0x01, 0x88, 0x00, 0x30, 0x00, 0x20, 0x00, 0x30,
             0x80, 0x2a, 0x69, 0x11, 0xcc, 0x50, 0x56, 0x49,
             0x94, 0xd6 },
           { 0x42, 0x5a, 0x68, 0x39, 0x31, 0x41, 0x59, 0x26,
             0x53, 0x59, 0xd5, 0x0d, 0x1c, 0x28, 0x00, 0x00,
             0x01, 0x51, 0x80, 0x00, 0x10, 0x01, 0x01, 0x80,
             0x02, 0x01, 0x80, 0x20, 0x00, 0x31, 0x0c, 0x01,
             0x06, 0x9b, 0x47, 0xe7, 0x38, 0x04, 0xa6, 0x28,
             0x2b, 0x24, 0xca, 0x6b },
           { 0x42, 0x5a, 0x68, 0x39, 0x31, 0x41, 0x59, 0x26,
             0x53, 0x59, 0x44, 0x1f, 0x23, 0x2f, 0x00, 0x00,
             0x02, 0x11, 0x00, 0x00, 0x00, 0xa5, 0xa2, 0xa0,
             0x00, 0x22, 0x06, 0x9a, 0x7a, 0x10, 0xc0, 0x8e,
             0x10, 0xd8, 0x43, 0x14, 0x15, 0x92, 0x65, 0x35 },
           { 0x42, 0x5a, 0x68, 0x39, 0x31, 0x41, 0x59, 0x26,
             0x53, 0x59, 0x83, 0x90, 0x27, 0xae, 0x00, 0x00,
             0x01, 0x81, 0x80, 0x0c, 0x00, 0x14, 0x20, 0x20,
             0x00, 0x21, 0x86, 0x81, 0x9a, 0x09, 0x4d, 0xa8,
             0xb9, 0x8a, 0x0a, 0xc9, 0x32, 0x9a },
           { 0x42, 0x5a, 0x68, 0x39, 0x31, 0x41, 0x59, 0x26,
             0x53, 0x59, 0x8d, 0x4f, 0x1e, 0x72, 0x00, 0x00,
             0x04, 0x41, 0x80, 0x40, 0x00, 0x00, 0x20, 0x14,
             0x60, 0x20, 0x00, 0x30, 0xc0, 0x08, 0x63, 0x45,
             0x84, 0x0f, 0xb8, 0xc5, 0x05, 0x64, 0x99, 0x4d },
           { 0x42, 0x5a, 0x68, 0x39, 0x31, 0x41, 0x59, 0x26,
             0x53, 0x59, 0x64, 0xf2, 0x3a, 0xdc, 0x00, 0x00,
             0x00, 0x9e, 0x00, 0x04, 0x00, 0x30, 0x00, 0x02,
             0x08, 0x08, 0x80, 0x20, 0x00, 0x31, 0x0c, 0x01,
             0x06, 0x99, 0xa4, 0xe4, 0x4c, 0x0a, 0x62, 0x82,
             0xb2, 0x4c, 0xa6 }
    };

    static const int eof = std::streambuf::traits_type::eof();

    struct pimpl: public bz_stream {};

    static inline int flush_macro(const flush_kind f) {
        switch (f) {
            case no_sync:
                return BZ_RUN;
            case full_sync:
                return BZ_FLUSH;
            case finish_sync:
                return BZ_FINISH;
            default:
                //should probably throw
                return BZ_RUN;
        }
    }

    common::common(std::streambuf *sb)
    : xstream::common_buffer(sb), z_strm(0), block_start(0), block_offset(0),
      streambuf_mutex(0)
    {
        LOG("bz::common");    

        z_strm = new pimpl;

        //initialize bzlib structure
        z_strm->bzalloc = NULL;
        z_strm->bzfree = NULL;
        z_strm->opaque = NULL;
        //buffers
        z_strm->avail_out = out.size;
        z_strm->next_out = out.buf;

        z_strm->avail_in = 0;
        z_strm->next_in = in.buf;
    }

    unsigned long int common::input_count() const {
        return ((uint64_t)(z_strm->total_in_hi32)<< 32) + (uint64_t)(z_strm->total_in_lo32);
    }

    unsigned long int common::output_count() const {
        return ((uint64_t)(z_strm->total_out_hi32)<< 32) + (uint64_t)(z_strm->total_out_lo32);
    }

    common::~common() {
        LOG("bz::~common");
        delete z_strm;
    }


    //default compression 9
    ostreambuf::ostreambuf(std::streambuf * sb)
    : common(sb), level(9) {
        LOG("bz::ostreambuf without compression level");
        block_start = _sb->pubseekoff(0, std::ios_base::cur, std::ios_base::out);
        init ();
    }

    ostreambuf::ostreambuf (std::streambuf * sb, int l)
    : common(sb), level(l) {
        LOG("bz::ostreambuf with compression level " << l);
        block_start = _sb->pubseekoff(0, std::ios_base::cur, std::ios_base::out);
        init ();
    }

    const char* error_str(int err) {
        switch(err) {
            case BZ_MEM_ERROR:
                return "out of memory";
            case BZ_CONFIG_ERROR:
                return "bzlib badly configured (bad sizes of int32 (4), int16 (2) or char (1), check and recompile)";
            case BZ_PARAM_ERROR:
                return "invalid parameter, possibly invalid compression level";
            case BZ_SEQUENCE_ERROR:
                return "bad sequence (this means xstream is buggy)";
            case BZ_DATA_ERROR:
                return "invalid or incomplete data (crc failed)";
            case BZ_DATA_ERROR_MAGIC:
                return "magic bytes not found in stream";
            case BZ_IO_ERROR:
                return "io error";
            case BZ_UNEXPECTED_EOF:
                return "premature end of data";
            case BZ_OUTBUFF_FULL:
                return "output buffer full";
        }
        
        return "unknown error";
    }

    void ostreambuf::raise_error(int err) {
        std::string what = error_str(err);

        LOG("bz::ostreambuf::raise_error (" << err << ") = " << what);

        if (what.size() > 0) {
            throw compress_error(this,what);
        } else {
            throw compress_error(this);
        }
    }


    void ostreambuf::init() {
        LOG("bz::ostreambuf::init");
        int cret =::BZ2_bzCompressInit(
            z_strm,
            level, 
            0, //verbosity
            30 //workFactor (default value) controls when to switch to the fallback algorithm
        );

        if (BZ_OK != cret) {
            LOG("bz::ostreambuf::init: error creating bz2stream " << cret);
            raise_error(cret);
        }
        //initialize streambuf interface functions
        setp(in.buf, in.buf + in.size);
    }

    ostreambuf::~ostreambuf() {
        LOG("bz::ostreambuf::~ostreambuf");
        //fullsync (write remaining data)
        flush(finish_sync);
 
        //sync underlying streambuf
        MUTEX_LOCK
        _sb->pubsync();
        MUTEX_UNLOCK

        if (0 != z_strm) {
            //XXX should I throw an exception in case of error?
            //remember this is a destructor
            int cret = ::BZ2_bzCompressEnd(z_strm);
            if (BZ_OK != cret) {
                LOG("\tERROR: BZ2_bzCompressEnd returned " << cret);
            }
        }
    }

    int ostreambuf::sync () {
        LOG("bz::ostreambuf::sync");
        int ret;
        MUTEX_LOCK
        ret = flush(finish_sync);
        _sb->pubsync();
        MUTEX_UNLOCK
        return ret;
    }


    int ostreambuf::overflow (int c) {
        LOG("bz::ostreambuf::overflow(" << c << ")\t available=" << (available ()) << "\tEOF=" << eof);
        if (eof == c) {
            LOG("\tEOF");
            flush(no_sync);
            return eof;
        } else {
            if (0 == available()) {
                LOG("\t have to flush :[]");
                flush(no_sync);
            }
            *pptr() = static_cast < char >(c);
            pbump(1);
        }
        return c;
    }

    std::streamsize ostreambuf::xsputn(const char *buffer, std::streamsize n) {
        LOG("bz::ostreambuf::xsputn(" << buffer << "," << n << ")");

        return flush(no_sync, buffer, n);
    }

    int ostreambuf::flush(flush_kind f, const char *appendbuf, int appendsize) {
        LOG("bz::ostreambuf::flush(" << f << ")");
        std::streamsize in_s = taken();
        LOG("\tinput_size=" << in_s);

        //set up compression engine input feed
        int written;
        if (in_s > 0) {
            z_strm->next_in = pbase();
            z_strm->avail_in = in_s;
            written = in_s;
        } else if (appendsize > 0) {
            z_strm->next_in = (char*)appendbuf;
            z_strm->avail_in = appendsize;
            written = appendsize;
            appendsize = 0;
        } else {
            z_strm->next_in = pbase();
            z_strm->avail_in = 0;
            written = 0;
        }
        block_offset += written;
        if (block_offset > (std::streamoff)level * 100000) {
            f = (f == no_sync)? finish_sync : f;
        }

        if (z_strm->avail_in + 
            z_strm->total_in_lo32 + z_strm->total_in_hi32 == 0)
        {
           return 0;
        }

        bool redo = false;
        bool reinit_compressor = false;

        do {
            int cret;
            redo = false;
            bool error = false;

            cret = ::BZ2_bzCompress(z_strm, flush_macro(f));

            //error handling
            if (finish_sync == f) {
                if (BZ_STREAM_END == cret) {
                    redo = false;
                    reinit_compressor = true;
                }
                else if (BZ_FINISH_OK == cret) {
                    redo = true;
                }
                else {
                    //serious error, throw exception
                    LOG("\terror in finish:" << cret);
                    error = true;
                }
            }
            else if (full_sync == f) {
                if (BZ_FLUSH_OK == cret) {
                    LOG("\tanother go at sync");
                    redo = true;
                }
                else if (BZ_RUN_OK == cret) {
                    LOG("\tsync ok");
                    redo = false;
                }
                else {
                    LOG("\terror in sync: " << cret);
                    error = true;
                }
            }
            else if (no_sync == f) {
                if (BZ_RUN_OK != cret) {
                    LOG("\terror compressing " << cret);
                    error = true;
                }
            }
            else {
                LOG("\tERROR: unknown flush mode " << flush_macro(f));
                throw general_error();
                error = true;
            }

            if (error) {
                raise_error(cret);
            }

            if (f == finish_sync) { // only complete streams can be written
                std::streamsize count = out.size - z_strm->avail_out;
                if (count > 0) {  // ignore empty blocks
                    LOG("\twriting " << count << " bytes");
                    int size = htonl(count);
                    MUTEX_LOCK
                    const std::streamsize wrote = _sb->sputn((char*)&size, 4) +
                                                  _sb->sputn(out.buf, count);
                    if (wrote != count + 4) {
                        MUTEX_ESCAPE
                        LOG("\terror writting, only wrote " << wrote 
                            << " but asked for " << count);
                        raise_error(BZ_IO_ERROR);
                    }
                    block_start = _sb->pubseekoff(0, std::ios_base::cur,
                                                     std::ios_base::out);
                    block_offset = 0;
                    MUTEX_UNLOCK
                }
                z_strm->next_out = out.buf;
                z_strm->avail_out = out.size;
            }

            if ((0 == z_strm->avail_out) && (0 != z_strm->avail_in)) {
                LOG("\tavail_out=0 => redo");
                redo = true;
            }

            if (!redo && appendbuf && appendsize > 0) {
                z_strm->next_in = (char*)appendbuf;
                z_strm->avail_in = appendsize;
                written += appendsize;
                appendsize = 0;
                redo = true;
            }
        } while (redo);
        assert (0 == z_strm->avail_in);

        if (reinit_compressor) {
            int cret;
            cret = ::BZ2_bzCompressEnd(z_strm);
            if (BZ_OK != cret) {
                LOG("\tERROR: BZ2_bzCompressEnd returned " << cret);
                raise_error(cret);
            }
            z_strm->bzalloc = NULL;
            z_strm->bzfree = NULL;
            z_strm->opaque = NULL;
            z_strm->avail_out = out.size;
            z_strm->next_out = out.buf;
            z_strm->avail_in = 0;
            z_strm->next_in = in.buf;
            cret =::BZ2_bzCompressInit(z_strm, level, 0, 30);
            if (BZ_OK != cret) {
                LOG("\tERROR: BZ2_bzCompressInit returned " << cret);
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

    istreambuf::istreambuf(std::streambuf *sb, int *left, unsigned int left_size)
    : common(sb), end(false), block_size(0), block_next(0), 
      new_block_start(0), new_block_offset(0),
      leftovers(0)
    {
        LOG("bz::istreambuf");
        int cret =::BZ2_bzDecompressInit(z_strm,
            0, //verbosity
            0  //no small memory
        );

        if (BZ_OK != cret) {
            LOG("\terror creating bz2stream " << cret);
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

    void istreambuf::raise_error(int err){
        std::string what = error_str(err);

        LOG("bz::istreambuf::raise_error (" << err << ") = " << what);

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
                z_strm->next_out = out.buf;
                z_strm->avail_out = 0;
                read_decompress();
            }
            while (block_offset < new_block_offset) {
                z_strm->next_out = out.buf;
                z_strm->avail_out = new_block_offset - block_offset;
                if (z_strm->avail_out > out.size) {
                    z_strm->avail_out = out.size;
                }
                block_offset += z_strm->avail_out;
                decompress();
            }
            new_block_start = 0;
            new_block_offset = 0;
        }

        z_strm->avail_out = out.size;
        z_strm->next_out = out.buf;
        if (0 < z_strm->avail_in) {
            LOG("\tdata in queue, inflating");
            decompress();
        }
        while (!end && z_strm->avail_out > 0) {
            read_decompress();
        }
        if (end && z_strm->avail_out > 0) {
            LOG("\tend of stream (EOF)");
            //signal the stream has reached it's end
            return eof;
        }
            
        //set streambuf pointers
        setg(out.buf, out.buf, z_strm->next_out);

        return int(out.buf[0]);
    }
    
    //read to buffer in place (apart from data already buffered)
    std::streamsize istreambuf::xsgetn(char *buffer, std::streamsize n) {
        LOG("bz::istreambuf::xsgetn (" << n << ")");

        if (new_block_start > 0 || new_block_offset > 0) {
            if (block_start != new_block_start ||
                block_offset > new_block_offset ||
                block_size == 0)
            {
                z_strm->next_out = out.buf;
                z_strm->avail_out = 0;
                read_decompress();
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
                z_strm->next_out = out.buf;
                z_strm->avail_out = new_block_offset - block_offset;
                if (z_strm->avail_out > out.size) {
                    z_strm->avail_out = out.size;
                }
                block_offset += z_strm->avail_out;
                decompress();
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

            z_strm->next_out = buffer + read;
            z_strm->avail_out = n - read;

            if (0 < z_strm->avail_in) {
                decompress();
            }
            while (!end && z_strm->avail_out > 0) {
                read_decompress();
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

    void istreambuf::read_decompress() {
        LOG("bz::istreambuf::read_decompress ");
        bool reinit_decompressor = false;
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
            reinit_decompressor = (block_next != block_start);
            read = leftovers->len;
            if (read < 4) {
                read += _sb->sgetn(leftovers->buf + read, 4 - read);
                if (read != 4) {
                    end = true;
                    MUTEX_ESCAPE
                    return;
                }
            }
            if (leftovers->buf[0] == 0) { // bz2 blocks have prefixed blocksize
                int *size = (int*)leftovers->buf;
                block_size = ntohl(*size);
                read -= 4;
                if (reinit_decompressor && read > 0) {
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
            else { // bz2 blocks are jammed together, no blocksize available
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

        // We want to be able to start decompression at an arbitrary position
        // in the input stream. This is possible with bzip2 streams, but there
        // is a problem that the compressed blocks are arbitrary numbers of 
        // bits long and they are catenated one after another in a bit stream
        // without any reference to byte boundaries. This makes it difficult 
        // to jump into the middle of a stream and start the decompressor 
        // since it expects a stream header followed by the first block that
        // happens to start on a byte boundary. To make this work, I splice
        // an artificial stream header followed by a dummy compressed block
        // onto the beginning of the input stream, where the length of the
        // dummy block in bits is chosen so that it abutts without padding 
        // the next block in the input stream. I have prepared 8 dummy blocks
        // so there should be one to match the alignment of any input block.
        // To match them, I look for the first place in the input stream
        // with a byte string that matches the last 5 bytes in one of my
        // dummy headers, and then I inject the dummy header into the 
        // decompressor ahead of the actual data. The dummy blocks are all
        // contrived to decompress to an 8-byte string, so throwing away the
        // first 8 bytes out of the decompressor, it is primed to decompress
        // the remaining stream without any need for bit-shifting the input.

        const char* head = (const char*)bz_header;
        if (reinit_decompressor) {
            // reinitialize bzlib structure
            int saved_buflen = z_strm->avail_out;
            char *saved_buffer = z_strm->next_out;
            int cret = ::BZ2_bzDecompressEnd(z_strm);
            if (BZ_OK != cret) {
                LOG("\tERROR: BZ2_bzDecompressEnd returned " << cret);
            }
            z_strm->bzalloc = NULL;
            z_strm->bzfree = NULL;
            z_strm->opaque = NULL;
            cret = ::BZ2_bzDecompressInit(z_strm, 0, 0);
            if (BZ_OK != cret) {
                LOG("\terror creating bz2stream " << cret);
                raise_error(cret);
            }
            if (strncmp(head, in.buf, 8) != 0) {
                int hdr;
                int splice;
                int match = 10;
                for (hdr = 0; hdr < 8; ++hdr) {
                    splice = bz_header_length[hdr] - 5;
                    const char* shead = (const char*)bz_header[hdr];
                    for (match = 1; match < 10; ++match) {
                        if (strncmp(&in.buf[match], &shead[splice], 5) == 0)
                            break;
                    }
                    if (match < 10)
                        break;
                }
                if (hdr > 7) {
                    LOG("\tbz2 stream format error on input");
                    raise_error(BZ_DATA_ERROR_MAGIC);
                }
                char dummy_buffer[10];
                z_strm->next_out = dummy_buffer;
                z_strm->avail_out = 8;
                if (hdr == 3)
                    z_strm->avail_out = 9;
                while (match > 0)
                    in.buf[--match] = (const char)bz_header[hdr][--splice];
                z_strm->avail_in = splice;
                z_strm->next_in = (char*)bz_header[hdr];
                decompress(); // waste the first 8 bytes
            }
            z_strm->avail_out = saved_buflen;
            z_strm->next_out = saved_buffer;
        }
        z_strm->next_in = in.buf;
        z_strm->avail_in = read;
        decompress();
    }

    void istreambuf::decompress() {
        LOG("bz::istreambuf::decompress ");

        int cret = ::BZ2_bzDecompress(z_strm);

        const int* head = (const int*)bz_header;
        int* buf = (int*)z_strm->next_in;

        if (BZ_STREAM_END == cret) {
            z_strm->avail_in = 0;
            block_next = 0;
        }
        else if (cret == BZ_DATA_ERROR && z_strm->avail_in == 0) {
            // Ignore CRC errors at the end of stream because we may not have
            // started decompressing at the beginning. We can rely on the CRC
            // checks that are present within each compressed block anyway.
            end = true;
        }
        else if (cret == BZ_DATA_ERROR && z_strm->avail_in == 4 && *buf == *head)
        {
            // Decompressor may complain that the first 4 bytes of the next
            // input block were appended to the previous block, if it had
            // a stream terminus. This is not the true end of the stream,
            // and not an error, just ignore it and reset block_next so
            // that the next buffer is treated as a new stream.
            z_strm->avail_in = 0;
            block_next = 0;
            end = false;
        }
        else if (BZ_OK != cret) {
            printf("bz2 input stream crapping out, cret is %d\n", cret);
            LOG("\terror decompressing: " << cret);
            raise_error(cret);
        }
    }

    istreambuf::~istreambuf() {
        LOG("bz::~istreambuf");
        if (0 != z_strm) {
            int cret = ::BZ2_bzDecompressEnd(z_strm);
            if (BZ_OK != cret) {
                LOG("\tERROR: BZ2_bzDecompressEnd returned " << cret);
            }
        }
    }

}//namespace bz
}//namespace xstream

#endif //bzlib
