/*! \file xstream/except/bz.h
 *
 * \brief exceptions related to bzlib usage xstream::bz namespace
 *
 */

#ifndef __XSTREAM_EXCEPT_BZ_H
#define __XSTREAM_EXCEPT_BZ_H

#include <xstream/config.h>

#include <string>
#include <xstream/except.h>
#include <xstream/bz.h>

namespace xstream{
    namespace bz{

/*!
 * \brief errors in bzlib usage
 *
 */
class general_error: public xstream::fatal_error
{
    public:
        general_error(
                const std::string& w="generic error in bzlib stream"
            )
            :xstream::fatal_error(w)
            {};
        virtual std::string module() const
        {
            return (xstream::fatal_error::module()+"::bzlib");
        }
};

/*!
 * \brief general bzlib compression errors
 *
 */


class compress_error: public general_error
{
    public:
        /*!
         * \brief ostreambuf that caused the exception 
         *
         * */
        xstream::bz::ostreambuf* stream;
        compress_error(
                xstream::bz::ostreambuf* p,
                const std::string& w
            )
            :general_error(w),stream(p){};

        compress_error(xstream::bz::ostreambuf* p)
            :general_error(),stream(p){};

        virtual std::string module() const{
            return (general_error::module()+"::compress");
        }
};


/*!
 * \brief general bzlib decompression errors
 *
 */

class decompress_error: public general_error {
    public:
        /*!
         * \brief istreambuf that caused the exception 
         *
         * */
        xstream::bz::istreambuf* stream;
        decompress_error(
                xstream::bz::istreambuf* p,
                const std::string& w
            )
            :general_error(w),stream(p)
            {};

        decompress_error(xstream::bz::istreambuf* p)
            :general_error(),stream(p)
            {};

        virtual std::string module() const
        {
            return (general_error::module()+"::decompress");
        }
};

}//namespace bz
}//namespace xstream

#endif
