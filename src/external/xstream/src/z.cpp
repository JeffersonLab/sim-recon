#include <xstream/config.h>

#if HAVE_LIBZ

#include <algorithm>
#include <string.h>
#include <string>

#include <xstream/z.h>
#include <xstream/except/z.h>
#include <stdexcept>

#include <stdio.h>
#include <zlib.h>
#include <arpa/inet.h>

#include <cassert>

#include "debug.h"

#define COMPRESSION_BLOCK_SIZE 32000

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
    : xstream::common_buffer(sb), z_strm(0), block_offset(0)
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
        init();
    }

    ostreambuf::ostreambuf(std::streambuf *sb, int l)
    : common(sb), level (l) {
        LOG ("z::ostreambuf with compression level " << l);
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
        _sb->pubsync();

        if (0 != z_strm) {
            //XXX should I throw an exception in case of error?
            //remember this is a destructor
            //I should definitly LOG something
            int cret = ::deflateEnd(z_strm);
            if (Z_OK != cret){
                LOG("z::~ostreambuf error dealocating zstream");
            }
        }
    }

    int ostreambuf::sync () {
      LOG ("z::ostreambuf::sync");
      int ret = flush(normal_sync);
      _sb->pubsync();
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
            f = (f == no_sync)? full_sync : f;
        }

        bool redo = false;

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
                if (Z_STREAM_END != cret) {
                    //serious error, throw exception
                    LOG ("\terror :" << cret);
                }
            } else if (Z_OK != cret) {
                LOG ("\terror deflating " << cret);
                //XXX throw exception here
                raise_error(cret);
            }

            if (f != no_sync) {
                std::streamsize count = out.size - z_strm->avail_out;
                LOG ("\twriting " << count << " bytes");
                int size = htonl(count);
                const std::streamsize wrote = _sb->sputn((char*)&size, 4) +
                                              _sb->sputn(out.buf, count);
                if (wrote != count + 4) {
                    LOG("\terror writting, only wrote " << wrote 
                        << " but asked for " << count);
                    raise_error(Z_STREAM_ERROR);
                }

                // reset output
                z_strm->next_out = reinterpret_cast < Bytef* >(out.buf);
                z_strm->avail_out = out.size;
                block_offset = 0;
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
        //reset buffer
        setp(in.buf, in.buf + in.size);
        return written;
    }

    /////////////////////
    // istream follows //
    /////////////////////

    istreambuf::istreambuf (std::streambuf * sb)
    : common(sb), end(false), block_size(0) {
        LOG ("z::istreambuf");

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
        int read;
        int block_pending = 0;
        if (block_size < 0) { // stream has no blocksize markers
            read = _sb->sgetn(in.buf, in.size);
        }
        else { // look for prefixed blocksize: leading byte = 0
            read = _sb->sgetn(in.buf, 4);
            if (in.buf[0] == 0) { // z blocks have prefixed blocksize
                int *size = (int*)in.buf;
                block_pending = ntohl(*size);
                read = _sb->sgetn(in.buf, block_pending);
            }
            else { // z blocks are jammed together, no blocksize available
                read += _sb->sgetn(in.buf + 4, in.size - 4);
                block_size = -1;
            }
        }
        LOG("\tread " << read << " bytes");
        block_offset = 0;

        if (0 == read) {
            LOG("\tpremature end of stream");
            //XXX throw exception
            raise_error(Z_DATA_ERROR);
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

        const char* head = (const char*)z_header;
        if (block_size == 0 && strncmp(head, in.buf, z_header_length) != 0) {
            z_strm->avail_in = z_header_length;
            z_strm->next_in = reinterpret_cast < Bytef* >(z_header);
            inflate(f); // inject the z stream header
        }
        block_size = read;
        z_strm->next_in = reinterpret_cast < Bytef* >(in.buf);
        z_strm->avail_in = read;
        inflate(f);
    }

    void istreambuf::inflate(const flush_kind f) {
        LOG("z::istreambuf::inflate " << f);

        int cret = ::inflate(z_strm, flush_macro(f));

        if (Z_STREAM_END == cret) {
            end = true;
        }
        else if (cret == Z_DATA_ERROR && z_strm->avail_in == 0) {
            // Ignore CRC errors at the end of stream because we may not have
            // started decompressing at the beginning. We can rely on the CRC
            // checks that are present within each compressed block anyway.
            cret = Z_OK;
            end = true;
        }
        else if (Z_OK != cret) {
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
