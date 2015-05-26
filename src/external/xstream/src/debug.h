/*! \file debug.h
 * 
 * \brief debugging/logging support
 *
 */


#ifndef XSTREAM_DEBUG_H
#define XSTREAM_DEBUG_H

#include <xstream/config.h>

#if ENABLE_LOGGING

#include <iostream>

namespace xstream{

extern std::ostream* debug;

} // namespace xstream

/*! \brief logging macro using an ostream
 */

#define LOG(d) do {(*debug) << "[" << (__FILE__) << ": " << __LINE__ << "]\t" << d << std::endl;} while (0)

#else

#define LOG(s)

#endif

#endif
