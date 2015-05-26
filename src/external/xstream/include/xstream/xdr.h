/*! \file xstream/xdr.h
 *
 * \brief C++ iostream like interface to read and write xdr streams
 * 
 * \arg endian safe
 * \arg needs floats and doubles stores in IEEE format in hardware
 * (most architecture work this way, \c SPARC is the only exception I know)
 *
 */


#ifndef __XSTREAM_XDR
#define __XSTREAM_XDR

#include <xstream/config.h>

#include <vector>
#include <string>

//for pair definition
#include <utility>

#include <streambuf>
#include <iostream>
#include <istream>
#include <ostream>

namespace{
    
#if HAVE_INTTYPES_H
#    include<inttypes.h>
#else
#    if HAVE_STDINT_H
#        include <stdint.h>
#    endif
#    error "I need inttypes.h or stdint.h to exist"    
#endif

}

/*!
 * \brief Library top namespace
 */
namespace xstream{
/*!
 * \brief xdr related objects
 */
namespace xdr{

using std::string;
using std::streambuf;
using std::pair;
using std::vector;

/*!
 * \brief Output xdr stream class
 *
 * Just a thin wrapper for streambuf to serialize all the datatypes
 *
 */

class ostream
{
    private:
        streambuf *_sb;
    public:
        /*!
         * \brief construct using a streambuf
         * 
         */
        ostream(streambuf* sb):_sb(sb){};

        /*!
         * \brief construct using an ostream
         * 
         */

        ostream(const std::ostream &os):_sb(os.rdbuf()){};

        ostream& operator<<(int32_t v);
        //ostream& operator<<(int v);
        ostream& operator<<(uint32_t v);
        //ostream& operator<<(unsigned int v);
        ostream& operator<<(int64_t v);
        //ostream& operator<<(const long int v);
        ostream& operator<<(uint64_t v);
        //ostream& operator<<(unsigned long v);
        ostream& operator<<(float v);
        ostream& operator<<(double v);
        ostream& operator<<(const string &v);

        /*!
         * \brief Serializes STL pair containers to xdr
         * 
         */
        template <typename A, typename B>
            ostream& operator<<(const pair<A,B> &p){
                (*this)<<(p.first);
                (*this)<<(p.second);
                return *this;
            }

        /*!
         * \brief Serializes STL vector containers to xdr
         * 
         */
        template <typename T>
            ostream& operator<<( const vector<T> &t){

            uint32_t sz = static_cast<uint32_t>( t.size() );
            (*this)<<sz;
            typename vector<T>::const_iterator cit = t.begin();
            const typename vector<T>::const_iterator cend = t.end();
            while(cit!=cend){
                (*this)<<*(cit++);
            }
            return (*this);
        }
};

/*!
 * \brief Input xdr stream class
 *
 * Just a thin wrapper for streambuf to deserialize all the datatypes
 *
 */
class istream
{
    private:
        streambuf *_sb;
        
    public:
        /*!
         * \brief construct using a streambuf
         * 
         */
        istream(streambuf* sb):_sb(sb){};

        /*!
         * \brief construct using an istream
         * 
         */
        istream(const std::istream &os):_sb(os.rdbuf()){};

        istream& operator>>(int32_t &v);
        //istream& operator>>(int &v);
        istream& operator>>(uint32_t &v);
        //istream& operator>>(unsigned int &v);
        istream& operator>>(int64_t &v);
        //istream& operator>>(long int &v);
        istream& operator>>(uint64_t &v);
        //istream& operator>>(unsigned long int &v);
        istream& operator>>(float &v);
        istream& operator>>(double &v);
        istream& operator>>(string &v);

        /*!
         * \brief Deserializes STL pair containers from xdr
         * 
         */
        template <typename A, typename B>
            istream& operator>>(pair<A,B> &p){
                (*this)>>(p.first);
                (*this)>>(p.second);
                return *this;
            }

        /*!
         * \brief Deserializes STL vector containers from xdr
         * 
         */
        template <typename T>
            istream& operator>>( vector<T> &t){
            uint32_t sz;
            T val;

            //maybe I could allocate all the necessary entries here
            
            (*this)>>sz;

            std::clog<<"VECTOR SIZE: "<<sz<<std::endl;

            while(sz--!=0){
                (*this)>>val;
                std::clog<<"read: "<<val<<std::endl;
                t.push_back(val);
            }
            return (*this);
        }
};

/*!
 * \example xdr_in.cpp
 *
 * shows how to use the xstream::xdr::istream to restore a structure from \c XDR data
 *
 * reads from standard input
 *
 */

/*!
 * \example xdr_out.cpp
 *
 * shows how to use the xstream::xdr::ostream to serialize a structure to \c XDR data
 *
 * writes to standard output
 *
 */

}//namespace xdr
}//namespace xstream

#endif
