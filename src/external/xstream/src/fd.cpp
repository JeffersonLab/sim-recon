#include <xstream/config.h>
#include <xstream/posix.h>
#include <xstream/fd.h>

#include <algorithm>
#include <streambuf>

#include "debug.h"

#if HAVE_FD

namespace xstream
{
    namespace fd
    {
/*
        streambuf::streambuf(xstream::posix::fd _fd):fd(_fd)
        {
            LOG("fd::streambuf (posix::fd)");
        }
*/
        static const int eof = std::streambuf::traits_type::eof();
        static const size_t buflen = 4000 * 1024;

        streambuf::streambuf(int f, bool c)
        : xstream::posix::fd(f, c), rbuf(buflen), wbuf(buflen)
        {
            LOG("fd::streambuf (int fd)");

            //next read/write calls will call uflow overflow
        }

        void streambuf::flush_write() {
            LOG("fd::streambuf::flush");
            write(wbuf.buf, taken());
            reset_write();
        }

        int streambuf::sync() {
            LOG("fd::streambuf::sync");

            //write remaining data
            flush_write();

            //physical sync data
            xstream::posix::fd::sync();
            //if there's an error an exception is thrown
            //and it doesn't get here
            return 0;
        }

        int streambuf::overflow(int c) {
            LOG("fd::streambuf::overflow (" << c << ")");

            if (eof == c) {
                return eof;
            }

            flush_write();
            *pptr() = static_cast < char >(c);
            pbump (1);
            return c;
        }

        std::streamsize streambuf::xsputn(const char* buffer, std::streamsize n)
        {
            flush_write();
            std::streamsize count = write(buffer, n);
            return count;
        }

        int streambuf::underflow() {
            LOG("fd::streambuf::underflow");
            
            std::streamsize nread = read(rbuf.buf,rbuf.size);

            setg(rbuf.buf, rbuf.buf, rbuf.buf + nread);

            if (0 == nread) {
                return eof;
            } else {
                return static_cast<int>(*rbuf.buf);
            }
        }

        std::streamsize streambuf::xsgetn(char *buffer, std::streamsize n) {
            LOG("fd::streambuf::xsgetn " << n);
            std::streamsize av = available();

            std::streamsize nread=0;

            char* beg = gptr();
            char* end = egptr();

            bool need_read=true;

            if (av >= n) {
                end = beg + n;
                need_read = false;
            }
            
            //copy buffered data to suplied buffer
            std::copy(beg, end, buffer);
            nread = end - beg;
            setg(rbuf.buf, end, rbuf.buf + rbuf.size);
        
            if (need_read) {
                nread += read(buffer + av, n - av);
            }
            LOG("\tread " << nread);
            return nread;
        }

        void streambuf::reset_write()
        {
            setp(wbuf.buf, wbuf.buf + wbuf.size);
        }

        streambuf::~streambuf() {
            LOG("fd::~streambuf");

            if (taken() > 0) {
                sync();
            }
        }
    }    //namespace fd
}    //namespace xstream

#endif
