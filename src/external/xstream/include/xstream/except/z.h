/*! \file xstream/except/z.h
 *
 * \brief exceptions related to zlib usage xstream::z namespace
 *
 */

#ifndef __XSTREAM_EXCEPT_Z_H
#define __XSTREAM_EXCEPT_Z_H

#include <xstream/config.h>

#include <string>
#include <xstream/except.h>
#include <xstream/z.h>

namespace xstream{
    namespace z{

/*!
 * \brief errors in zlib usage
 *
 */
class general_error: public xstream::fatal_error
{
    public:
        general_error(
                const std::string& w="generic error in zlib stream"
            )
            :xstream::fatal_error(w){};
        virtual std::string module() const
        {
            return (xstream::fatal_error::module()+"::zlib");
        }
};

/*!
 * \brief general zlib compression errors
 *
 */


class compress_error: public general_error
{
    public:
        /*!
         * \brief ostreambuf that caused the exception 
         *
         * */
        xstream::z::ostreambuf* stream;
        compress_error(
                xstream::z::ostreambuf* p,
                const std::string& w
            )
            :general_error(w),stream(p)
            {};

        compress_error(xstream::z::ostreambuf* p)
            :general_error(),stream(p)
            {};

        virtual std::string module() const
        {
            return (general_error::module()+"::compress");
        }
};


/*!
 * \brief general zlib decompression errors
 *
 */

class decompress_error: public general_error {
    public:
        /*!
         * \brief istreambuf that caused the exception 
         *
         * */
        xstream::z::istreambuf* stream;
        decompress_error(
                xstream::z::istreambuf* p,
                const std::string& w
            )
            :general_error(w),stream(p){};

        decompress_error(xstream::z::istreambuf* p)
            :general_error(),stream(p){};

        virtual std::string module() const{
            return (general_error::module()+"::decompress");
        }
};

}//namespace z
}//namespace xstream

#endif
