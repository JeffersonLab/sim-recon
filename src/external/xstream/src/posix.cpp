#include <xstream/config.h>
#include <xstream/posix.h>
#include <xstream/except/posix.h>
#include <iosfwd>

#include <string>
#include <ctime>

//for read write, etc
#include <unistd.h>

//for errno
#include <errno.h>

//for strerror and strerror_r
#include <cstring>

#include "debug.h"

/*
 * I don't prefix time localtime_r and strftime with std:: because lot's of implementaions are broken and it's an unnecessary purism
 *
 */

namespace xstream {
namespace posix{

    date_format::date_format(const std::string& f):format(f){};

    std::string date_format::now() const {
        LOG("posix::date_format::now");
        ::time_t t;
        ::time(&t);
        
        struct std::tm tm;

        //XXX check error and throw exception
        ::localtime_r(&t,&tm);

        //first size of the allocated buffer
        //if needed it will be reallocated
        size_t len=512;
        std::string ret;

        do {
            LOG("\tlen=" << len);
            char buf[len];
            size_t cret=::strftime(buf,len,format.c_str(),&tm);

            if (0 != cret && cret <= len) {
                buf[len - 1] = '\0'; //just in case
                ret = buf;
            }
            len *= 2;
        } while (0 == ret.size());

        LOG("\tdate=" << ret);
        return ret;
    }

    void check_return(int code, const std::string& call) {
        LOG("posix::check_return " << call << " => " << code);
        if (-1 == code) {
            //XXX please try to use strerror_r instead
            const std::string desc(strerror(errno));
            LOG("\tthrowing " << errno << " => " << desc);
            throw general_error(call, errno, desc);
        }
    }

#if HAVE_FD
    fd::fd(int f, bool c)
    : fdn(f), dest_close(c)
    {
        LOG("posix::fd (" << f << "," << c << ")");
    }
    
    std::streamsize fd::read(char* buf, std::streamsize len)
    {
        LOG("posix::fd::read " << len);
        
        ssize_t count;
        
        do {
            count = ::read(fdn, buf, len);
        } while (-1 == count && EINTR == errno);

        check_return(count, "read");

        return count;
    }
    
    std::string fd::read(std::streamsize len)
    {
        char buf[len];
        ssize_t count = read(buf, len);

        return std::string(buf, buf + count);
    }
    
    std::streamsize fd::write(const char* buf, std::streamsize len)
    {
        LOG("posix::fd::write " << len);
        
        ssize_t count;
        
        do {
            count = ::write(fdn, buf, len);
        } while (-1 == count && EINTR == errno);

        check_return(count, "write");

        return count;
    }
    
    void fd::sync()
    {
        LOG("posix::fd::sync");
        int cret = fsync(fdn);
        check_return(cret, "fsync");
    }

    fd::~fd()
    {
        LOG("posix::fd::~fd");
        if (dest_close) {
            LOG("\tclosing");
            int cret = ::close(fdn);
            check_return(cret, "close");
        }
    }
#endif

}//namespace posix
}//namespace xstream
