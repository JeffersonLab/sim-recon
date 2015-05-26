/* RFC 1321 compliante MD% implementation
 *
 * little-endian optimizations but works on all byte orderings
 *
 */

#include <xstream/config.h>
#include <xstream/digest.h>

#include <iosfwd>
#include <iomanip>
#include <iostream>

#include "debug.h"

//see md5_t.pl to understand T values

#include "md5_t.h"

//initial values of digest

#define A0 0x01234567
#define B0 0x89abcdef
#define C0 0xfedcba98
#define D0 0x76543210

//auxiliary functions

static inline uint32_t rotate_left(uint32_t x, unsigned int n)
{
    return ((x << n) | (x >> (32 - n)));
}

static inline uint32_t F(uint32_t x, uint32_t y, uint32_t z)
{
    return ((x & y) | ((~x) & z));
}

static inline uint32_t G(uint32_t x, uint32_t y, uint32_t z)
{
    return ((x & z) | (y & (~z)));
}

static inline uint32_t H(uint32_t x, uint32_t y, uint32_t z)
{
    return (x ^ y ^ z);
}

static inline uint32_t I(uint32_t x, uint32_t y, uint32_t z)
{
    return (y ^ (x | (~z)));
}

static const size_t block_size = 16 * /*4 byte words*/ 4;

namespace xstream{
namespace digest{

    //this is where the actual digest work is made
    //buf should be 32*16 bytes long
    static void process_chunk(
        uint32_t& AA,
        uint32_t& BB,
        uint32_t& CC,
        uint32_t& DD,
        char* _buf)
    {
        LOG("md5::process_chunk (A,B,C,D) = (" << AA << "," << BB << "," << CC << "," << DD << ")");

        //hope register is at least tolerated by most compilers
        register uint32_t A=AA;
        register uint32_t B=BB;
        register uint32_t C=CC;
        register uint32_t D=DD;
        
        //needed for endian agnostic copy bellow
        unsigned char* buf = reinterpret_cast<unsigned char*>(_buf);

#if ARCH_LITTLE_ENDIAN
        //on some systems this can be tricky due to allignment issues
        //some compiler flags may solve that
        //if it doesn' work ok, try configure without ARCH_LITTLE_ENDIAN
        uint32_t* wbuf = reinterpret_cast<uint32_t*>(buf);
#else
        uint32_t wbuf[block_size/4]; //word buffer
        //we have to copy the words by a byte ordering agnostic procedure

        for (unsigned int i=0, j=0; i < block_size; i += 4, ++j) {
            //LOG("\t(i,j) = ("<<i<<","<<j<<")");
            wbuf[j] = buf[i] + (buf[i+1] << 8) + (buf[i+2] << 16) + (buf[i+3] << 24);
        }
#endif

        //now that we have the words we process the word buffer

        //wish I could remove these defines and use templates instead
        
#define BRA(q,w,e,r,k,s,Ti) \
        q = w + ( rotate_left( q + FUN(w,e,r) + wbuf[k] + Ti, s)  );

#define FUN F

        BRA(A, B, C, D, 0, 7, T1);
        BRA(D, A, B, C, 1, 12, T2);
        BRA(C, D, A, B, 2, 17, T3);
        BRA(B, C, D, A, 3, 22, T4);
        BRA(A, B, C, D, 4, 7, T5);
        BRA(D, A, B, C, 5, 12, T6);
        BRA(C, D, A, B, 6, 17, T7);
        BRA(B, C, D, A, 7, 22, T8);
        BRA(A, B, C, D, 8, 7, T9);
        BRA(D, A, B, C, 9, 12, T10);
        BRA(C, D, A, B, 10, 17, T11);
        BRA(B, C, D, A, 11, 22, T12);
        BRA(A, B, C, D, 12, 7, T13);
        BRA(D, A, B, C, 13, 12, T14);
        BRA(C, D, A, B, 14, 17, T15);
        BRA(B, C, D, A, 15, 22, T16);

#undef FUN
#define FUN G

        BRA(A, B, C, D, 1, 5, T17);
        BRA(D, A, B, C, 6, 9, T18);
        BRA(C, D, A, B, 11, 14, T19);
        BRA(B, C, D, A, 0, 20, T20);
        BRA(A, B, C, D, 5, 5, T21);
        BRA(D, A, B, C, 10, 9, T22);
        BRA(C, D, A, B, 15, 14, T23);
        BRA(B, C, D, A, 4, 20, T24);
        BRA(A, B, C, D, 9, 5, T25);
        BRA(D, A, B, C, 14, 9, T26);
        BRA(C, D, A, B, 3, 14, T27);
        BRA(B, C, D, A, 8, 20, T28);
        BRA(A, B, C, D, 13, 5, T29);
        BRA(D, A, B, C, 2, 9, T30);
        BRA(C, D, A, B, 7, 14, T31);
        BRA(B, C, D, A, 12, 20, T32);

#undef FUN
#define FUN H

        BRA(A, B, C, D, 5, 4, T33);
        BRA(D, A, B, C, 8, 11, T34);
        BRA(C, D, A, B, 11, 16, T35);
        BRA(B, C, D, A, 14, 23, T36);
        BRA(A, B, C, D, 1, 4, T37);
        BRA(D, A, B, C, 4, 11, T38);
        BRA(C, D, A, B, 7, 16, T39);
        BRA(B, C, D, A, 10, 23, T40);
        BRA(A, B, C, D, 13, 4, T41);
        BRA(D, A, B, C, 0, 11, T42);
        BRA(C, D, A, B, 3, 16, T43);
        BRA(B, C, D, A, 6, 23, T44);
        BRA(A, B, C, D, 9, 4, T45);
        BRA(D, A, B, C, 12, 11, T46);
        BRA(C, D, A, B, 15, 16, T47);
        BRA(B, C, D, A, 2, 23, T48);

#undef FUN
#define FUN I

        BRA(A, B, C, D, 0, 6, T49);
        BRA(D, A, B, C, 7, 10, T50);
        BRA(C, D, A, B, 14, 15, T51);
        BRA(B, C, D, A, 5, 21, T52);
        BRA(A, B, C, D, 12, 6, T53);
        BRA(D, A, B, C, 3, 10, T54);
        BRA(C, D, A, B, 10, 15, T55);
        BRA(B, C, D, A, 1, 21, T56);
        BRA(A, B, C, D, 8, 6, T57);
        BRA(D, A, B, C, 15, 10, T58);
        BRA(C, D, A, B, 6, 15, T59);
        BRA(B, C, D, A, 13, 21, T60);
        BRA(A, B, C, D, 4, 6, T61);
        BRA(D, A, B, C, 11, 10, T62);
        BRA(C, D, A, B, 2, 15, T63);
        BRA(B, C, D, A, 9, 21, T64);

        //end of "brackets"

        LOG("\tchunk values (A,B,C,D) = (" << A << "," << B << "," << C << "," << D << ")");

        AA += A;
        BB += B;
        CC += C;
        DD += D;
    }

    md5::md5()
    : block_stream(block_size,4)
    {
        LOG("digest::md5");
        reset_digest();
    }

    void md5::reset_digest()
    {
        LOG("digest::md5::reset_digest");
        result.a = 0x67452301;
        result.b = 0xefcdab89;
        result.c = 0x98badcfe;
        result.d = 0x10325476;
    }

    void md5::calculate_digest()
    {
        LOG("digest::md5::calculate_digest");
        process_chunk(
            result.a,
            result.b,
            result.c,
            result.d,
            pbase()
        );
    }

    struct md5::result md5::digest()
    {
        LOG("digest::md5::digest");
        //now for the real deal;
        pubsync();
        
        //cache current digest
        struct result d = result;
        
        const unsigned long int t = taken();
        const unsigned long int l = length + t;
        const char* orig = pbase();

        //I could make this a litle more efficient, but since this only occurs at the end, maybe it's ok
        char b[block_size * 2];

        std::copy(orig, orig + t, b);

        //pad data
        std::fill(b + t, b + 2 * block_size, '\0');
        b[t] = 128; //add one bit after data

        LOG("\ttaken = " << t << "\tlen=" << l);

        //write size of data
        unsigned long int ll = l * 8; //length in bits not bytes
        char* const end = (b + block_size - 8) + (t >= block_size - 8 ? block_size : 0);
        LOG("\tend-b = " << (end - b));
        
        for (int i=0; i < 8; ++i){
            end[i] = ll & ((1 << 8) - 1);
            ll >>= 8;
        }

        LOG("\tprocessing chunks");
        for(char* ptr=b; ptr < end; ptr += block_size) {
            LOG("\tend-ptr = " << (end - ptr));
            process_chunk(d.a, d.b, d.c, d.d, ptr);
        }
        return d;
    }

    struct md5::result md5::reset()
    {
        LOG("digest::md5::reset");
        struct result r = digest();
        reset_digest();
        return r;
    }

    static void print_hex(std::ostream& o, uint32_t n)
    {
        for (int i=0; i < 4; ++i) {
            o.width(2);
            o.fill('0');
            o << (n & ((1 << 8) - 1));
            n >>= 8;
        }
    }

    std::ostream& operator<<(std::ostream& o, const struct md5::result& r)
    {
        std::ios::fmtflags orig = o.flags();
        o.setf(std::ios::hex, std::ios::basefield);

        print_hex(o, r.a);
        print_hex(o, r.b);
        print_hex(o, r.c);
        print_hex(o, r.d);
        
        o.flags(orig);

        return o;
    }

    }//namespace digest
}//namespace xstream
