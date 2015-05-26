/*! \file xstream/dater.h
 *
 * \brief C++ objects to date output lines
 */

#ifndef __XSTREAM_DATER_H
#define __XSTREAM_DATER_H

#include <xstream/config.h>

#include <xstream/common.h>
#include <xstream/posix.h>

#include <streambuf>
#include <string>

namespace xstream{

/*!
 * \brief filter that add a time stamp on the begining of each line
 *
 */
class dater: public xstream::ostreambuf
{
    private:
        std::streambuf* _sb; /*!< streambuf were data is written to */
        posix::date_format date;
        char separator; /*!< data is written after this character */
        bool write_next;
        
        /*!
         * \brief writes the current date to the streambuf
         *
         */
        void write_date();
        
        /*!
         * \brief tries to write as much data as possible (overloaded from streambuf)
         *
         * */
        int sync();

        /*!
         * \brief write a character that surpasses buffer end (overloaded from streambuf)
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
         * \brief construct using another streambuf
         *
         * \param sb streambuf to write data to
         * \param format format string for \c strftime defaults to [%c] the preferred date and time representation for the current locale
         * \param separator character that signal a change of line, defaults to newline
         * \param write_first if true date is written right before the first write. default is true
         *
         */

        dater(std::streambuf* sb, const std::string& format="[%c] ", char separator='\n', bool write_first=true);
        
        /*! destructor
         *
         * */
        ~dater();
};

}//namespace xstream

#endif
