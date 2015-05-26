/*! \file xstream/except/base64.h
 *
 * \brief exception related to base64 encode/decode, xstream::base64 namespace
 *
 */

#ifndef __XSTREAM_EXCEPT_BASE64_H
#define __XSTREAM_EXCEPT_BASE64_H

#include <xstream/config.h>

#include <string>
#include <xstream/except.h>
#include <xstream/z.h>

namespace xstream{
    namespace base64{

/*!
 * \brief errors in base64 usage
 *
 */
class general_error: public xstream::fatal_error
{
    public:
        general_error(
                const std::string& w="generic error in base64"
            )
            :xstream::fatal_error(w){};
        virtual std::string module() const{
            return (xstream::fatal_error::module()+"::base64");
        }
};

/*!
 * \brief general base64 encoding errors
 *
 */


class encode_error: public general_error
{
    public:
        encode_error(const std::string& w):general_error(w){};

        virtual std::string module() const{
            return (general_error::module()+"::encode");
        }
};


/*!
 * \brief general base64 decoding errors
 *
 */

class decode_error: public general_error
{
    public:
        decode_error(const std::string& w):general_error(w){};

        virtual std::string module() const{
            return (general_error::module()+"::decode");
        }
};

}//namespace base64
}//namespace xstream

#endif
