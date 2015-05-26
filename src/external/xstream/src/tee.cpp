#include <xstream/tee.h>

#include <iostream>
#include <map>

#include "debug.h"

namespace xstream
{
    namespace tee
    {

        static const int eof = std::streambuf::traits_type::eof();

        void ostreambuf::add( std::streambuf* sb) {
            LOG("tee::ostreambuf::add (streambuf*)");
            destinations.insert(sb);
        }

        void ostreambuf::add( std::ostream& os) {
            LOG("tee::ostreambuf::add (ostream&)");
            destinations.insert(os.rdbuf());
        }

        void ostreambuf::remove( std::streambuf* sb) {
            LOG("tee::ostreambuf::remove (streambuf*)");
            destinations.erase(sb);
        }

        void ostreambuf::remove( std::ostream& os) {
            LOG("tee::ostreambuf::remove (ostream&)");
            destinations.erase(os.rdbuf());
        }

        int ostreambuf::sync() {
            LOG("tee::ostreambuf::sync");
            if (0 == destinations.size()) {
                return eof;
            }

            std::map<std::streambuf*, int> return_code;

            std::set<std::streambuf*>::iterator it = destinations.begin();
            for (; it != destinations.end(); it++) {
                int ret = (*it)->pubsync();
                return_code[*it] = ret;
            }

            //remove streambufs with errors
            //second go

            it = destinations.begin();
            for (;it != destinations.end(); it++) {
                if (eof == return_code[*it]) {
                    destinations.erase(it);
                }
            }

            if (0 == destinations.size()) {
                return eof;
            }

            return 0;
        }

        int ostreambuf::overflow(int c) {
            LOG("tee::ostreambuf::overflow " << c);

            if (0 == destinations.size()) {
                return eof;
            }

            std::map<std::streambuf*, int> return_code;

            std::set<std::streambuf*>::iterator it = destinations.begin();
            for (;it != destinations.end(); it++) {
                int ret = (*it)->sputc(c);
                return_code[*it] = ret;
            }

            //remove streambufs with errors
            //second go

            it = destinations.begin();
            for (;it != destinations.end(); it++) {
                if (eof==return_code[*it]) {
                    destinations.erase(it);
                }
            }

            if (0 == destinations.size()) {
                return eof;
             }

            return c;
        }

        std::streamsize ostreambuf::xsputn(const char *buffer, std::streamsize n) {
            LOG("tee::ostreambuf::xsputn " << buffer);

            if (0 == destinations.size()) {
                return eof;
            }

            std::map<std::streambuf*, int> return_code;

            std::set<std::streambuf*>::iterator it = destinations.begin();
            for (;it != destinations.end(); it++) {
                int ret = (*it)->sputn(buffer,n);
                return_code[*it] = ret;
            }

            //remove streambufs with errors
            //second go

            it = destinations.begin();
            for (;it != destinations.end(); it++) {
                if (0 > return_code[*it]) {
                    destinations.erase(it);
                }
            }

            if (0 == destinations.size()) {
                return eof;
            }

            return n;
        }

        ostreambuf::~ostreambuf() {
            //no need to delete anything
        }

    }                                    //namespace fork
}                                        //namespace xstream
