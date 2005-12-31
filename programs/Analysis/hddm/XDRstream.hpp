#ifndef _XDRSTREAM_
#define _XDRSTREAM_
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

#include <stdio.h>
#include <stdlib.h>
#include <stream.h>
#include <strstream.h>
#include <rpc/types.h>
#include <rpc/xdr.h>

class XDRstream
{
public:
  virtual int eof() 	{ return iostr->eof(); }
  virtual int fail() 	{ return iostr->fail(); }
  virtual int bad()	{ return iostr->bad(); }
  virtual int good()	{ return iostr->good(); }
  virtual int rdstate() { return iostr->rdstate(); }     
  virtual void clear(int v=0) { iostr->clear(v); }
  virtual ~XDRstream() = 0;
protected:
  XDRstream(): iostr(NULL) {}
  ios *iostr;
};


class XDRistream: public XDRstream
{
public:
  XDRistream(char *, size_t size);
  ~XDRistream();
  XDRistream&	get(char& c);
  XDRistream&	read(char*& s, int n);
  XDRistream&	read(u_char*& s, int n);

  XDRistream&	operator>>(char*&  );
  XDRistream&	operator>>(u_char*&);
  XDRistream&	operator>>(char&  );
  XDRistream&	operator>>(double&);
  XDRistream&	operator>>(float& );
  XDRistream&	operator>>(int&   );
  XDRistream&	operator>>(long&  );
  XDRistream&	operator>>(short& );
  XDRistream&	operator>>(unsigned char& );
  XDRistream&	operator>>(unsigned int&  );
  XDRistream&	operator>>(unsigned long& );
  XDRistream&	operator>>(unsigned short&);
private:  
  XDRistream& operator=(const XDRistream&) { return *this;}
protected:
  XDR xdr;
  char *buf;
  strstream i;			// this is used only for set/get status
};

class XDRostream: public XDRstream
{
public:
  XDRostream(streambuf *);
  ~XDRostream();

  XDRostream&	operator<<(const char*);
  XDRostream&	operator<<(char);
  XDRostream&	operator<<(u_char);
  XDRostream&	operator<<(double);
  XDRostream&	operator<<(float);
  XDRostream&	operator<<(int);
  XDRostream&	operator<<(u_int);
  XDRostream&	operator<<(long);
  XDRostream&	operator<<(u_long);
  XDRostream&	operator<<(short);
  XDRostream&	operator<<(u_short);
  
  XDRostream&	put(char c);

  XDRostream&	write(const char* s, int n);
  XDRostream&	write(const u_char* s, int n);

  // tellp is for debugging only
  int tellp() { return o.tellp(); }
  
private:  
  XDRostream& operator=(const XDRostream&) { return *this;}

protected:
  XDR xdr;
  char *buf;			// mem-XDR buffer
  int size;			// current buffer size
  ostream o;
};


#endif
