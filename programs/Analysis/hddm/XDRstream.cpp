//
// Copyright (C) 1994-1995 Harry Danilevsky and Andrew Renalds
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish and distribute copies of the 
// Software, and to permit persons to whom the Software is furnished to do so,
// subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESSED OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// HARRY DANILEVSKY AND ANDREW RENALDS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//

#include "XDRstream.hpp"
#include <stdio.h>
#include <string.h>
#include <generic.h>
#include <limits.h>

static int const XDRBUFSIZE=1024;

XDRstream::~XDRstream() {}


// Input stream

XDRistream::XDRistream(char *buf, size_t size)
{
  iostr = &i;
  xdrmem_create(&xdr, buf, size, XDR_DECODE);
}

XDRistream::~XDRistream()
{
  xdr_destroy(&xdr);
}

#define XdrGetOp(type, xdrtype)						      \
XDRistream&								      \
XDRistream::operator>>(type& c)					      \
{									      \
    if (!name2(xdr_,xdrtype(&xdr, &c)))					      \
      clear( rdstate() | ios::failbit );				      \
    return *this;							      \
}									      \


XdrGetOp(char,char)
XdrGetOp(double,double)
XdrGetOp(float,float)
XdrGetOp(int,int)
XdrGetOp(long,long)
XdrGetOp(short,short)

XdrGetOp(unsigned int,u_int)
XdrGetOp(unsigned long,u_long)
XdrGetOp(unsigned short,u_short)


// operator >>(char *) reads null-terminated string from incoming stream
// Since xdr stores string size I don't need a get(char *, size_t) method.
// If s is not NULL, it is assumed to point to buffer of sufficient size.
// Otherwise memory will be allocated with new and it will be caller's
// responsibilty to delete s;     
     
     
XDRistream&
XDRistream::operator >> (char*& s)
{
  char *tmps = NULL;
  
  if (!xdr_wrapstring(&xdr, &tmps))
    clear( rdstate() | ios::failbit );
  else
    {
      if (!s)
	s = new char[strlen(tmps)+1];
      strcpy(s, tmps);
      free(tmps);
    }
  return *this;
}

XDRistream&
XDRistream::operator >> (u_char*& s)
{
  return XDRistream::operator >> ((char*&) s);
}

// read N bytes from input stream (which must be written out by XDRostream::write()).
// If input array is larger than N bytes, extra bytes will be ignored.
// If p is NULL, memory will be allocated with new and it will be caller's
// responsibilty to delete p;  

XDRistream&
XDRistream::read(char*& p, int N)
{
  char *tmps = NULL;
  u_int size, maxsize = UINT_MAX;
  
  if (!xdr_bytes(&xdr, &tmps, &size, maxsize))
    clear( rdstate() | ios::failbit );
  else if (tmps)
    {
      if (!p)
	p = new char[size];
      memcpy (p, tmps, N);
      free(tmps);
    }
  return *this;
}

XDRistream&
XDRistream::read(u_char*& p, int N)
{
  return XDRistream::read((char*&) p, N);
}

// Output Stream

XDRostream::XDRostream(streambuf *sb): o(sb)
{
  iostr = &o;
  buf = new char[size=XDRBUFSIZE];
  xdrmem_create(&xdr, buf, size, XDR_ENCODE);
}

XDRostream::~XDRostream()
{
  xdr_destroy(&xdr);
  delete buf;
}

#define XdrPutOp(type,xdrtype)				 		      \
XDRostream&								      \
XDRostream::operator << (type c)					      \
{									      \
  if (name2(xdr_,xdrtype(&xdr, &c)))					      \
    o.write(buf, xdr_getpos(&xdr));					      \
  else									      \
    clear( rdstate() | ios::failbit );					      \
  xdr_setpos(&xdr, 0);							      \
  return *this;								      \
}

XdrPutOp(char,char)
XdrPutOp(double,double)
XdrPutOp(int,int)
XdrPutOp(unsigned int,u_int)
XdrPutOp(long,long)
XdrPutOp(unsigned long,u_long)
XdrPutOp(short,short)
XdrPutOp(unsigned short,u_short)

// Special version for floats to avoid problem with conversion
// between floats and doubles     
     
XDRostream&
XDRostream::operator<<(float c)
{
  float temp = c;
  if (xdr_float(&xdr, &temp))
    o.write(buf, xdr_getpos(&xdr));
  else    
    clear( rdstate() | ios::failbit );
  xdr_setpos(&xdr, 0);
  return *this;
}

// For C strings we may have to reallocate xdr memory buffer.

XDRostream&
XDRostream::operator << (const char* s)
{
  if (s)
    {
      // realloc buf is necessary
      // Note that xdr internally prepends strings with integer length,
      // that's why extra 4 bytes_per_unit are aXDRed
      int newsize = RNDUP(strlen(s))+BYTES_PER_XDR_UNIT;
      if (size < newsize)
	{
	  cout << "XDRostream::<<(const char *) reallocing buf. old size: "
	       << size << "\tnewsize: " << newsize << endl;
	  delete buf;
	  xdr_destroy(&xdr);
	  buf = new char[size=newsize];
	  xdrmem_create(&xdr, buf, size, XDR_ENCODE);
	}
    if (xdr_string(&xdr, (char**)&s, strlen(s)))
      o.write(buf, xdr_getpos(&xdr));
    else
      clear( rdstate() | ios::failbit );
    }
  xdr_setpos(&xdr, 0);
  return *this;
}

XDRostream&
XDRostream::put(char c)
{
  return (*this << c);
}


XDRostream&
XDRostream::write(const char* p, int N)
{
  if (p)
    {
      /* realloc buf is necessary */

      int newsize = RNDUP(N)+BYTES_PER_XDR_UNIT;

      if (size < newsize)
	{
	  cout << "XDRostream::write reallocing buf. old size: "
	       << size << "\tnewsize: " << newsize << endl;
	  delete buf;
	  xdr_destroy(&xdr);
	  buf = new char[size=newsize];
	  xdrmem_create(&xdr, buf, size, XDR_ENCODE);
	}
      if (xdr_bytes(&xdr, (char **) &p, (u_int *)&N, N))
	o.write(buf, xdr_getpos(&xdr));
      else
	clear( rdstate() | ios::failbit );
    }
  xdr_setpos(&xdr, 0);
  return *this;
}

XDRostream&
XDRostream::write(const u_char *p, int N)
{
  return XDRostream::write((const char *)p, N);
}
