#include <xstream/config.h>
#include <xstream/base64.h>
#include <xstream/except/base64.h>

#include "debug.h"

namespace xstream
{
    namespace base64
    {

        static const char dictionary[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

        static const int eof = std::streambuf::traits_type::eof();

        ostreambuf::ostreambuf(std::streambuf* sb, unsigned int w, char c)
        : _sb(sb), delim(c), delim_w(w), col(0) {
            LOG("base64::ostreambuf (streambuf*)");
            reset();
        }

        void ostreambuf::reset() {
            LOG("base64::ostreambuf::reset");
            setp(buf,buf+sizeof(buf)/sizeof(char));
        }

        //no flush of base64 data, just flushes the streambuf it writes to
        int ostreambuf::sync() {
            LOG("base64::ostreambuf::sync");
            return _sb->pubsync();
        }

        static inline void encode(const char* in, char* out) {
            LOG("base64::encode");
            static const unsigned int hi6 = ((1 << 6) -1 ) << 2;
            static const unsigned int lo2 = (1 << 2) - 1;
            static const unsigned int hi4 = ((1 << 4) - 1) << 4;
            static const unsigned int lo4 = (1 << 4) - 1;
            static const unsigned int hi2 = ((1 << 2) - 1) << 6;
            static const unsigned int lo6 = (1 << 6) - 1;

            out[0] = dictionary[(hi6 & in[0]) >> 2];
            out[1] = dictionary[(lo2 & in[0]) << 4 | (hi4 & in[1]) >> 4];
            out[2] = dictionary[(lo4 & in[1]) << 2 | (hi2 & in[2]) >> 6];
            out[3] = dictionary[lo6 & in[2]];
        }

        int ostreambuf::write(const char* buf, size_t len) {
            unsigned int rcol = delim_w - col;
            int ret = 0;

            //XXX all these sputn calls need to be checked
            if (delim_w > 0 && rcol <= len) {
                LOG("\t" << rcol << " columns to padding");
                ret += _sb->sputn(buf, rcol);
                ret += (_sb->sputc(delim) == delim)? 1 : 0;
                ret += _sb->sputn(buf + rcol, len - rcol);
                col = len - rcol;
                LOG("\tcol = " << col);
            }
            else {
                ret += _sb->sputn(buf,len);
                col += len;
            }

            reset();

            return ret;
        }

        int ostreambuf::overflow(int c) {
            LOG("base64::ostreambuf::overflow (" << c << ")\n\t[" << *buf << *(buf+1) << *(buf+2) << "]");

            char enc[4];                //encoded buffer
            encode(buf, enc);

            write(enc, 4);

            *pptr() = static_cast < char >(c);
            pbump(1);

            return c;
        }

        ostreambuf::~ostreambuf() {
            LOG("base64::~ostreambuf");
            int av = available();

            if (3 == av) {
                //no need to do anything
            }
            else {
                char enc[4];
                for(int i=0; i < av; ++i) {
                    //pad to zero
                    buf[2-i] = '\0';
                }
                encode(buf,enc);

                if (1 <= av) {
                    enc[3] = '=';
                }
                if (2 <= av) {
                    enc[2] = '=';
                }

                //XXX need to check return code
                write(enc,4);
            }

            _sb->pubsync();
        }

        //should search in constant time
        static inline char index(char c) {
            static const int window = 'Z'- 'A' + 1 ;

            //if c is not in dict throw exception
            //could be a lot cleaner if no debugging was needed

            int d = c - 'A';
            int ret = -1;

            if (d >= 0 && d < window) {
                ret = d;
            }
            else {
                d = c - 'a';
                if (d >= 0 && d < window) {
                    ret = d + window;
                }
                else {
                    d = c - '0';
                    if (d >= 0 && d < 10) {
                        ret = d + 2 * window;
                    }
                    else {
                        if ('+' == c) {
                            ret = 62;
                        }
                        else {
                            if ('/' == c) {
                                ret = 63;
                            }
                        }
                    }
                }
            }

            if (-1 == ret) {
                LOG("base64::index (" << c << ") [unknown]");
                throw(decode_error(std::string("character '") + c + "' not part of base64 alphabet"));
            }
            else {
                return ret;
            }
        }

        static inline void decode(const char* in, char* out) {
            LOG("base64::decode");
            static const unsigned char hi2 = ((1 << 2) - 1) << 4;
            static const unsigned char lo4 = ((1 << 4) - 1);
            static const unsigned char hi4 = ((1 << 4) - 1) << 2;
            static const unsigned char lo2 = ((1 << 2) - 1);

            char _in[4];

            for (int i=0; i < 4; ++i) {
                _in[i] = index(in[i]);
            }

            char c;
            c = _in[0] << 2;
            c |= (_in[1] & hi2) >> 4;
            out[0] = c;

            c = (_in[1] & lo4) << 4;
            c |= (_in[2] & hi4) >> 2;
            out[1] = c;

            c = (_in[2] & lo2) << 6;
            c |= _in[3];
            out[2] = c;
        }

        istreambuf::istreambuf(std::streambuf* sb, unsigned int d_w, char d)
        : _sb(sb), end(false), delim(d), delim_w(d_w), col(0) {
            LOG("base64::istreambuf (streambuf*)");
            setg(buf, buf, buf);
        }

        int istreambuf::underflow() {
            LOG("base64::istreambuf::underflow");

            if (end) {
                LOG("\tattempt to read from an ended stream");
                return eof;
            }
            else {
                char enc[4 + 1];
                int av = (delim_w - col);
                if (delim_w > 0 && av <= 4) {
                    LOG("\texpecting delimiter at " << av);
                    int ret = _sb->sgetn(enc, 5);
                    LOG("\tenc=" << enc);
                    if (5 != ret) {
                        LOG("\tread " << ret << " but wanted " << 5);
                        //XXX throw exception?
                        end = true;
                        return eof;
                    }
                    if (delim != enc[av]) {
                        LOG("\texpected delimiter " << delim << " but got " << enc[av]);
                        end = true;
                        throw decode_error("expected delimiter missing");
                        return eof;
                    }
                    else {
                        col = 4 - av;
                        for (int i=av; i < 4; ++i) {
                            enc[i] = enc[i + 1];
                        }
                    }
                }
                else {
                    int ret = _sb->sgetn(enc, 4);
                    col += 4;

                    if (4 != ret) {
                        //XXX maybe I should throw an exception
                        end = true;
                        return eof;
                    }
                }
                int len=3;
                if ('=' == enc[3]) {
                    enc[3] = 'A';    //guard
                    if ('=' == enc[2]) {
                        len = 1;
                        enc[2] = 'A';    //guard
                    } else {
                        len = 2;
                    }
                }
                if (3 != len) {
                    end = true;
                }

                decode(enc, buf);
                setg(buf, buf, buf + len);
                return 0;
            }
        }

        istreambuf::~istreambuf() {
            LOG("base64::~istreambuf");
        }

    }//namespace base64
}//namespace xstream
