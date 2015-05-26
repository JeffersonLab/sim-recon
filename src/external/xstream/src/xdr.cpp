#include <xstream/xdr.h>


/*
 * Minimal xdr implementation
 * should be endian agnostic
 *
 * XXX need to check for EOF, and have some kind of failbit
 * XXX throw exceptions
 */

namespace xstream {
    namespace xdr{

// OUTPUT
    
ostream& ostream::operator<<(const string &s) {
    static const char* cpad = "\0\0\0";
    size_t len =  s.size();
    // mudar isto XXX
    (*this) << (static_cast<uint32_t>(len));
    _sb->sputn(s.c_str(), len);
    //mod 4 via bitmask (only works for powers of 2)
    size_t pad = len & 3;
    if (pad > 0){
        _sb->sputn(cpad, 4 - pad);
    }
    return *this;
}

ostream& ostream::operator<<(uint32_t _v) {
    uint32_t v = _v;
    unsigned char c[4];
    for (int i=0; i < 4; i++, v >>= 8) {
        c[i] = static_cast<unsigned char>(v & 0xff);
    }
    for (int i=0; i < 4; i++) {        // RFC1832 mandates msb...lsb order
        _sb->sputc(c[3 - i]);
    }
    return *this;
}

/* ostream& ostream::operator<<(unsigned int v) {
    uint32_t n = *(reinterpret_cast<const uint32_t*>(&v));
    return (*this << n);
}
*/

ostream& ostream::operator<<(int32_t v) {
    uint32_t n = *(reinterpret_cast<const uint32_t*>(&v));
    return (*this << n);
}

/*
ostream& ostream::operator<<(int v) {
    uint32_t n = *(reinterpret_cast<const uint32_t*>(&v));
    return (*this << n);
}
*/

ostream& ostream::operator<<(int64_t v) {
    uint64_t n = *(reinterpret_cast<const uint64_t*>(&v));
    return (*this << n);
}

/*
ostream& ostream::operator<<(const long int v) {
    uint64_t n = *(reinterpret_cast<const uint64_t*>(&v));
    return (*this << n);
}
*/

ostream& ostream::operator<<(uint64_t _v) {
    uint64_t v=_v;
    unsigned char c[8];
    for (int i=0; i < 8; i++, v >>= 8) {
        c[i] = static_cast<unsigned char>(v & 0xff);
    }
    for (int i=0; i < 8; i++) {        // RFC1832 mandates msb...lsb order
        _sb->sputc(c[7 - i]);
    }
    return *this;
}

/*
ostream& ostream::operator<<(unsigned long int v) {
    uint64_t n = *(reinterpret_cast<const uint64_t*>(&v));
    return (*this << n);
}
*/

/*
ostream& ostream::operator<<(float v) {
    uint32_t n = *(reinterpret_cast<const uint32_t*>(&v));
    return (*this << n);
}
*/

// assume floats on this platform are IEEE-754 binary32
ostream& ostream::operator<<(float v) {
    union jack {
        float binary32;
        uint32_t ui32;
    } _v;
    _v.binary32 = v;
    return (*this << _v.ui32);
}

/*
ostream& ostream::operator<<(double v) {
    uint64_t n = *(reinterpret_cast<const uint64_t*>(&v));
    return (*this << n);
}
*/

// assume doubles on this platform are IEEE-754 binary32
ostream& ostream::operator<<(double v) {
    union jack {
        float binary64;
        uint32_t ui64;
    } _v;
    _v.binary64 = v;
    return (*this << _v.ui64);
}

// INPUT

istream& istream::operator>>(string &s) {
    uint32_t len;
    (*this) >> len;

    char line[len];
    _sb->sgetn(line, len);
    s = string(line, line + len);
    
    //mod 4 via bitmask (only works for powers of 2)
    size_t pad = (4 - len) & 3;
    //ignore padding zeros
    char dummy[pad];
    _sb->sgetn(dummy, pad);
    return *this;
}

istream& istream::operator>>(uint32_t &v) {
    v=0;
    for (int i=0; i < 32; i += 8) {
        uint32_t c = static_cast<uint32_t>(_sb->sbumpc());
        v <<= 8;
        v |= c;
    }
    return *this;
}

/*
istream& ostream::operator>>(unsigned int &v) {
    uint32_t _v;
    (*this) >> _v;
    v=reinterpret_cast<unsigned int>(_v);
    return (*this);
}
*/

istream& istream::operator>>(int32_t &v) {
    uint32_t _v;
    (*this) >> _v;
    v = static_cast<int32_t>(_v);
    return (*this);
}

/*
istream& istream::operator>>(int &v) {
    uint32_t _v;
    (*this) >> _v;
    v=reinterpret_cast<int>(_v);
    return (*this);
}
*/

istream& istream::operator>>(uint64_t &v) {
    v = 0;
    for (int i=0; i < 64; i += 8) {
        uint64_t c = static_cast<uint64_t>(_sb->sbumpc());
        v <<= 8;
        v |= c;
    }
    return *this;
}

/*
istream& istream::operator>>(unsigned long int &v) {
    uint32_t _v;
    (*this) >> _v;
    v = reinterpret_cast<unsigned long int>(_v);
    return (*this);
}
*/

istream& istream::operator>>(int64_t &v) {
    uint64_t _v;
    (*this) >> _v;
    v = static_cast<int64_t>(_v);
    return (*this);
}

/*
istream& istream::operator>>(long int &v) {
    uint32_t _v;
    (*this) >> _v;
    v=reinterpret_cast<long int>(_v);
    return (*this);
}
*/

istream& istream::operator>>(float &v) {
    uint32_t n;
    (*this) >> n;
    float* vp = reinterpret_cast<float*>(&n);
    v = *vp;
    return *this;
}

istream& istream::operator>>(double &v) {
    uint64_t n;
    (*this) >> n;
    double* vp = reinterpret_cast<double*>(&n);
    v = *vp;
    return *this;
}

}//namespace xdr
}//namespace xstream
