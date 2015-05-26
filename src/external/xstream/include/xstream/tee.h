/*! \file xstream/tee.h
 *
 * \brief C++ streambuf to fork output
 * data written to it is also written to several other streambufs
 *
 */

#ifndef __XSTREAM_TEE_H
#define __XSTREAM_TEE_H

#include <xstream/config.h>

#include <xstream/common.h>
#include <streambuf>
#include <set>

namespace xstream{

/*!
 * \brief  tee related ostreambuf objects
 *
 */
namespace tee{

/*!
 * \brief fork output stream class
 *
 * when data is written to it, it writes it to several other streambufs
 * error handling is as follows:
 *  \arg if writting to a streambuf returns an error, that streambuf is removed from the internal list and no more data is written to it
 *  \arg when there is no streambuf left to write to \c eof is returned
 *
 *  \note doesn't catch any exceptions, so let's the io library decide if they should be propagated
 *
 */
class ostreambuf: public xstream::ostreambuf
{
    private:

        std::set<std::streambuf*> destinations; /*!< set of streambufs to write to */

        /*!
         * \brief flush as much data as possible (overloaded from streambuf)
         *
         * */
        int sync();

        /*!
         * \brief write a character that supasses buffer end (overloaded from streambuf)
         * 
         */
        int overflow(int c);

        /*!
         * \brief write an entire buffer (overloaded from streambuf)
         *
         */
        std::streamsize xsputn(const char *buffer, std::streamsize n);
        
    public:
        /*!
         * \brief construct \c NOP object
         */
        ostreambuf()
        {};

        /*!
         * \brief add an output streembuf to write to
         */
        void add(std::streambuf* sb);

        /*!
         * \brief remove a streambuf to write to
         *
         */
        void remove(std::streambuf* sb);

        /*!
         * \brief add an output ostream to write to
         * 
         * higher level syntatic sugar for the streambuf version
         */
        void add(std::ostream& os);

        /*!
         * \brief remove an output ostream to write to
         * 
         * higher level syntatic sugar for the streambuf version
         */
        void remove(std::ostream& os);

        /*!
         * \brief closes the streambuf stream
         *
         */
        ~ostreambuf();

};

 /*!
  * \example n_tee.cpp
  *
  * shows how to write the content of standard input to several files using the xstream::tee::ostreambuf class
  *
  */


}//namespace fork
}//namespace xstream

#endif
