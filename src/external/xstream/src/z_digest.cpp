#include <xstream/config.h>
#include <xstream/digest.h>

#include "debug.h"

#if HAVE_LIBZ

//for adler32 and crc32
#include <zlib.h>


namespace xstream{
namespace digest{

    ////////////////////////////////
    /////// Adler32 digest /////////
    ////////////////////////////////
    
    static const size_t buffer_size = 4000 * 1024;
    typedef unsigned long int uli;
    
    z_common::z_common()
    : common<uli>(buffer_size), _digest(0)
    {
        LOG("digest::z_common");
    }

    unsigned long int z_common::digest()
    {
        LOG("digest::z_common::digest");
        return _digest;
    }

    void z_common::reset_digest()
    {
        _digest=0;
    }

    void adler32::calculate_digest()
    {
        LOG("digest::adler32::calculate_digest");
        _digest = ::adler32(_digest, reinterpret_cast<Bytef*>(pbase()), taken());
        LOG("\tdigest = " << _digest);
    }

    //////////////////////////////
    /////// CRC32 digest /////////
    //////////////////////////////
    

    void crc32::calculate_digest()
    {
        LOG("digest::crc32::calculate_digest");
        _digest = ::crc32(_digest, reinterpret_cast<Bytef*>(pbase()), taken());
        LOG("\tdigest = " << _digest);
    }

    }//namespace digest
}//namespace xstream


#endif
