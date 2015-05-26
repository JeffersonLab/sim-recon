#include <xstream/common.h>
#include <algorithm>

#include "debug.h"

namespace xstream {

    static const size_t buffer_size     = 4000 * 1024;
    static const size_t out_buffer_size = static_cast < size_t > (buffer_size * 1.02 + 12);    //see manual

    buffer::buffer(size_t s)
    : buf(0), size(0)
    {
        LOG("buffer::buffer ("<<s<<")");
        resize(s);
    }

    void buffer::grow(unsigned int f)
    {
        LOG("buffer::grow " << f);
        if (f < 1) {
            LOG("\tERROR: just tried to grow to a smaller size");
            return;
        }
        const size_t new_s = size * f;
        char* new_b = new char[new_s];

        std::copy(buf, buf + size, new_b);
        delete[] buf;

        size = new_s;
        buf = new_b;
    }

    void buffer::resize(size_t s)
    {
        LOG("buffer::resize " << s);
        if (buf) {
            LOG("\tdeleting buf");
            delete[] buf;
        }
        size = s;
        buf = new char[size];
    }

    buffer::~buffer()
    {
        LOG("buffer::~buffer");    
        delete[] buf;
    }

    common_buffer::common_buffer(std::streambuf * sb)
    :_sb (sb), in(buffer_size), out(out_buffer_size)
    {
        LOG("common_buffer");    
    }

    common_buffer::~common_buffer() {
            LOG("~common_buffer");    
    }

}//namespace xstream
