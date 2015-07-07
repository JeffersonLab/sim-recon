#include <xstream/dater.h>

#include <string>
#include <streambuf>

#include "debug.h"

namespace xstream {

    dater::dater(std::streambuf* sb, const std::string& f, char sep, bool wn)
        : _sb(sb), date(f), separator(sep), write_next(wn)
    {
        LOG("dater::dater format=" << f);
    }

    int dater::sync() {
        LOG("dater::sync");
        return _sb->pubsync();
    }

    int dater::overflow(int c) {
        LOG("dater::overflow (" << c << ")");
        if (write_next) {
            write_date();
        }

        //XXX error handling could improve here

        _sb->sputc(c);
        
        if (separator == c) {
            write_next = true;
        }
        return c;
    }

    std::streamsize dater::xsputn(const char *buffer, std::streamsize n) {
        LOG("dater::xsputn " << buffer);
        std::streamsize written=0;

        const char* end_buffer = buffer + n;
        const char* beg;
        const char* end;

        std::string d;
        for (beg = buffer, end = buffer; end < end_buffer; ++end) {
            if (write_next) {
                write_date();
            }

            if (separator == *end) {
                //XXX check return code
                std::streamsize w = end - beg + 1;
                _sb->sputn(beg, w);
                written += w;
                beg = end + 1;
                write_next = true;
            }
        }
        return written;
    }

    void dater::write_date() {
        LOG("dater::write_date");
        //XXX check return code
        std::string d = date.now();
        _sb->sputn(d.c_str(), d.size());
        write_next = false;
    }

    dater::~dater() {
        LOG("dater::~dater");
        _sb->pubsync();
    }
}//namespace xstream
