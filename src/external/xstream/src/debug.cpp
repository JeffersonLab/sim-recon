#include "debug.h"
#include <fstream>
#include <iomanip>
#include <cstdlib>

int AGlobalVariableToPreventWarningAboutNoSymbols;

#if ENABLE_LOGGING


namespace xstream{

std::ostream* debug = &std::clog;

#if 0
    class log {
        log()
        {
            char* fname = getenv("XSTREAM_LOGFILE");
            if (0 != fname) {
                debug = new std::fstream(fname, std::ios::out | std::ios::trunc);
                (*debug) << std::ios::unitbuf;
            }
        }
    };

    log _log_();

#endif

} //namespace xstream


#endif
